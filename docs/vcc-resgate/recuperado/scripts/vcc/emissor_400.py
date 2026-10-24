# scripts/vcc/emissor_400.py
"""
Etapa 6 de docs/VCC-2026-ROTA.md: `emissor_400` — o gerador das celulas da
submissao. Recebe um delta de pseudobulk previsto por (contexto, alvo) e devolve
400 celulas de contagem crua por par, escritas em streaming num h5ad.

MECANISMO. E o do `baseline.py` da Arc (`_emission_scale` + `_emit_scaled_resample`),
com uma unica adicao nossa — o arredondamento:

  1. sorteia 400 posicoes do pool de controle, uniforme com reposicao;
  2. multiplica cada gene da celula sorteada por uma escala de gene;
  3. arredonda para inteiro (a emissao `dispersed` da Arc devolve float32);
  4. elimina os zeros que o arredondamento criou.

O passo 1 e o que torna a emissao estocastica, e isso nao e detalhe de estilo:
400 celulas identicas (o modo `tile` da Arc) zeram a covariancia prevista e
produzem ~9.701 genes "significativos" por perturbacao contra ~2.739 na
referencia real. Afunda `jaccard` e `fidelity`, duas das seis metricas.

A ESCALA DE GENE, que e onde o delta entra. O delta previsto vive no espaco
log1p de CP50k, o mesmo de `b_ctrl`. Em proporcao de gene:

    p_ctrl = expm1(b_ctrl) / 5e4                      (soma 1,0 — medido)
    p_alvo = expm1(b_ctrl + delta) / 5e4,  renormalizado para somar 1
    escala = p_alvo / p_ctrl,  e 0 onde p_ctrl == 0

A renormalizacao nao e cosmetica: sem ela um delta que move muitos genes para o
mesmo lado desloca o tamanho de biblioteca inteiro, e `nmae` / `mse` pagam por
isso. Com a escala assim, o fator de biblioteca media da celula se cancela — e
exatamente o `profile_g / ctrl_pb_g` da referencia, escrito em proporcao em vez
de em contagem media, porque e nessa unidade que o delta chega.

ARREDONDAMENTO ESTOCASTICO, nao `rint`. `floor(v) + Bernoulli(v - floor(v))`
preserva a media por gene exatamente; `rint` nao preserva (uma escala de 0,8
sobre contagem 3 da 2,4, que `rint` manda para 2 — vies de -0,4 por entrada, no
mesmo sentido em todo gene reprimido). O custo e uma variancia de Bernoulli de
no maximo 0,25 por entrada, desprezivel contra a variancia de contagem que a
celula ja carrega.

GATE. `runs/vcc/eixo/gate_alvos.tsv` marca (contexto, alvo) onde o alvo nao esta
transcrito. CRISPRi reprime transcricao: sem transcrito nao ha o que derrubar, e
a predicao ali e o proprio controle. O delta desses pares e zerado.

STREAMING. O h5ad e aberto por anndata sem X (obs e var sao pequenas), e o grupo
`/X` e criado a mao com `data`/`indices` redimensionaveis. Cada bloco de 400
celulas e apendado. A matriz inteira nunca existe: um bloco de 400 x 18.533 em
CSR e ~27 MB; o pico de RAM e o proprio pool de controle, ~0,9 GB.

ENTRADA DE DELTA (`--delta`, um .npz). Chaves esperadas:

    genes      (18533,)  str   — tem de ser gene_names.csv na ordem
    alvos      (300,)    str   — o conjunto de pert_counts.csv
    contextos  (k,)      str   — subconjunto de {D, E, F}
    delta_D    (300, 18533) float — delta em log1p-CP50k, linha i = alvos[i]
    delta_E, delta_F   idem

Alternativa aceita no lugar de `delta_<ctx>`: `lfc_<ctx>`, log2 da razao de
proporcao do gene. E convertida aqui para o mesmo caminho
(p_alvo proporcional a p_ctrl * 2^lfc, renormalizado), entao as duas rotas
coincidem.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/emissor_400.py --stub --contextos D
    ./.venv/Scripts/python.exe -I scripts/vcc/emissor_400.py --delta runs/vcc/transfer_fit/delta.npz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import cast

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
import scipy.sparse as sp

VCC = Path("data/raw/vcc/controls-final")
CTRL_BULK = Path("runs/vcc/ctrl_bulk/ctrl_bulk.npz")
GATE_ALVOS = Path("runs/vcc/eixo/gate_alvos.tsv")
SAIDA = Path("runs/vcc/submissao")

# O cell-eval2 normaliza o pseudobulk para 50.000 contagens antes do log1p; e o
# espaco em que b_ctrl e o delta vivem.
ESCALA_PSEUDOBULK = 5e4
CELULAS_POR_ALVO = 400
COLUNA_ALVO = "target_gene"
COLUNA_CONTEXTO = "context"
ROTULO_CONTROLE = "non-targeting"
CONTEXTOS_FINAIS = ("D", "E", "F")

# Limites medidos do contrato, de vcc/prep.py e vcc/sizing.py da CLI 0.3.0.
MAX_CONTAGEM_POR_CELULA = 1_000_000
MAX_CELULAS = 400_000
TETO_NNZ = 4_750_000_000

# Chunk dos datasets HDF5 redimensionaveis: 1 Mi entradas = 4 MB por bloco de
# compressao, grande o bastante para o gzip render e pequeno o bastante para o
# cache de chunk do h5py nao crescer.
CHUNK_HDF5 = 1 << 20
NIVEL_GZIP = 4


@dataclass
class Medicao:
    """Numeros que a etapa 6 tem de reportar, acumulados durante a emissao."""

    celulas: int = 0
    nnz: int = 0
    max_contagem_celula: float = 0.0
    min_contagem_celula: float = float("inf")
    nnz_por_bloco: list[int] = field(default_factory=list)
    massa_inalcancavel: list[float] = field(default_factory=list)
    escala_maxima: float = 0.0
    pares_gateados: list[str] = field(default_factory=list)


def _por_gene(valores: object, n_genes: int, nome: str) -> np.ndarray:
    """Converte uma entrada do .npz num vetor float64 de n_genes, ou falha."""
    vetor = np.asarray(valores, dtype=np.float64).ravel()
    if vetor.size != n_genes:
        raise ValueError(f"{nome}: {vetor.size} genes, esperado {n_genes}")
    return vetor


def ler_eixos() -> tuple[list[str], list[str]]:
    """Eixo de genes da saida e lista oficial de alvos, das fontes de verdade."""
    genes = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
    alvos = pd.read_csv(VCC / "pert_counts.csv")["target_gene"].astype(str).tolist()
    if len(set(genes)) != len(genes):
        raise ValueError("gene_names.csv tem simbolo repetido; o eixo de saida nao e valido")
    if ROTULO_CONTROLE in alvos:
        raise ValueError(f"pert_counts.csv lista '{ROTULO_CONTROLE}', que nao se submete")
    return genes, alvos


def ler_gate(alvos: list[str], contextos: tuple[str, ...]) -> dict[tuple[str, str], bool]:
    """Mapa (contexto, alvo) -> alvo transcrito ali. Ausencia do arquivo e erro.

    Falhar e melhor que assumir `ativo`: o gate existe porque predizer efeito num
    alvo sem transcrito gasta as metricas de significancia em troca de nada, e um
    default silencioso de `ativo` desfaria a etapa 1 sem aviso.
    """
    if not GATE_ALVOS.exists():
        raise FileNotFoundError(f"{GATE_ALVOS} nao existe; rode a etapa 1 (eixo_18262) antes")
    tabela = pd.read_csv(GATE_ALVOS, sep="\t")
    mapa: dict[tuple[str, str], bool] = {}
    for contexto, alvo, gate in zip(
        tabela["contexto"].astype(str),
        tabela["alvo"].astype(str),
        tabela["gate"].astype(str),
        strict=True,
    ):
        mapa[(contexto, alvo)] = gate == "ativo"
    faltando = [
        (c, a) for c in contextos for a in alvos if (c, a) not in mapa
    ]
    if faltando:
        raise ValueError(
            f"gate_alvos.tsv nao cobre {len(faltando)} pares, p.ex. {faltando[:5]}"
        )
    return mapa


def carregar_delta(
    caminho: Path, genes: list[str], alvos: list[str], contextos: tuple[str, ...]
) -> tuple[dict[str, np.ndarray], bool]:
    """Le o .npz do `transfer_fit`: matriz (n_alvos, n_genes) por contexto.

    Devolve tambem se o conteudo e log-razao de proporcao (`lfc_<ctx>`) em vez de
    delta aditivo (`delta_<ctx>`). As linhas saem na ordem de `alvos`, qualquer que
    seja a ordem do arquivo — o reindex e feito aqui para que o emissor nao dependa
    de o modelo ter escolhido a mesma ordem.

    Misturar as duas familias de chave no mesmo arquivo e erro, nao conveniencia:
    um unico sinalizador decide a conversao de todos os contextos, e adivinhar por
    contexto esconderia uma troca de unidade num deles.
    """
    with np.load(caminho, allow_pickle=True) as arquivo:
        chaves = set(arquivo.files)
        for obrigatoria in ("genes", "alvos"):
            if obrigatoria not in chaves:
                raise ValueError(f"{caminho}: falta a chave '{obrigatoria}'")
        genes_npz = np.asarray(arquivo["genes"]).astype(str).ravel().tolist()
        if genes_npz != genes:
            raise ValueError(
                f"{caminho}: 'genes' tem {len(genes_npz)} entradas e nao reproduz "
                f"gene_names.csv ({len(genes)}) na ordem; o eixo de saida e o do painel"
            )
        alvos_npz = np.asarray(arquivo["alvos"]).astype(str).ravel().tolist()
        if set(alvos_npz) != set(alvos):
            faltando = sorted(set(alvos) - set(alvos_npz))[:5]
            extras = sorted(set(alvos_npz) - set(alvos))[:5]
            raise ValueError(
                f"{caminho}: 'alvos' difere de pert_counts.csv (faltando {faltando}, "
                f"extras {extras})"
            )
        posicao = {alvo: i for i, alvo in enumerate(alvos_npz)}
        reindex = np.array([posicao[alvo] for alvo in alvos], dtype=np.intp)

        com_delta = [c for c in contextos if f"delta_{c}" in chaves]
        com_lfc = [c for c in contextos if f"lfc_{c}" in chaves]
        if com_delta and com_lfc:
            raise ValueError(
                f"{caminho}: mistura delta_{com_delta} com lfc_{com_lfc}. Um arquivo "
                "carrega uma unidade so"
            )
        if not com_delta and not com_lfc:
            raise ValueError(
                f"{caminho}: falta delta_<ctx> (ou lfc_<ctx>) para {list(contextos)}"
            )
        e_log_razao = bool(com_lfc)
        prefixo = "lfc" if e_log_razao else "delta"
        faltando = [c for c in contextos if f"{prefixo}_{c}" not in chaves]
        if faltando:
            raise ValueError(f"{caminho}: falta {prefixo}_<ctx> para {faltando}")

        saida: dict[str, np.ndarray] = {}
        for contexto in contextos:
            chave = f"{prefixo}_{contexto}"
            bruto = np.asarray(arquivo[chave], dtype=np.float64)
            if bruto.shape != (len(alvos_npz), len(genes)):
                raise ValueError(
                    f"{caminho}: '{chave}' tem forma {bruto.shape}, esperado "
                    f"{(len(alvos_npz), len(genes))}"
                )
            if not np.isfinite(bruto).all():
                ruins = int((~np.isfinite(bruto)).sum())
                raise ValueError(f"{caminho}: '{chave}' tem {ruins} valores nao finitos")
            # `lfc` nao e convertido aqui: a razao de proporcao so fecha contra o
            # `p_ctrl` do contexto, que `escala_de_gene` tem e este leitor nao.
            # O que viaja e o log natural do fator; a composicao e renormalizacao
            # acontecem la, num unico lugar.
            saida[contexto] = bruto[reindex] * (np.log(2.0) if e_log_razao else 1.0)
            print(f"  {chave}: forma {saida[contexto].shape}")
    return saida, e_log_razao


def escala_de_gene(
    b_ctrl: np.ndarray, delta: np.ndarray, *, delta_e_log_razao: bool
) -> tuple[np.ndarray, float, float]:
    """Escala multiplicativa por gene, mais a massa inalcancavel e a escala maxima.

    `delta_e_log_razao=True` significa que `delta` e o log natural de uma razao de
    proporcao (o que `_lfc_para_delta` devolve); False significa que ja e o delta
    aditivo em log1p-CP50k.
    """
    p_ctrl = np.expm1(b_ctrl) / ESCALA_PSEUDOBULK
    if delta_e_log_razao:
        p_alvo = p_ctrl * np.exp(delta)
    else:
        p_alvo = np.expm1(b_ctrl + delta) / ESCALA_PSEUDOBULK
    p_alvo = np.maximum(p_alvo, 0.0)
    soma = float(p_alvo.sum())
    if not np.isfinite(soma) or soma <= 0.0:
        raise ValueError(f"perfil alvo soma {soma}; delta impossivel de renormalizar")
    p_alvo = p_alvo / soma

    ctrl_zero = p_ctrl <= 0.0
    # Divisor seguro antes do `where`: numpy avalia os dois lados da divisao.
    denominador = np.where(ctrl_zero, 1.0, p_ctrl)
    escala = np.where(ctrl_zero, 0.0, p_alvo / denominador)
    # Massa que o perfil alvo pede em gene que o controle nunca expressa: nenhuma
    # celula sorteada tem entrada ali, entao essa fracao da predicao e perdida por
    # construcao. Vale medir, nao corrigir — inventar contagem em gene mudo seria
    # predizer expressao onde o pool nao tem nenhuma evidencia.
    inalcancavel = float(p_alvo[ctrl_zero].sum())
    return escala, inalcancavel, float(escala.max())


def arredondar_estocastico(valores: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """floor(v) + Bernoulli(v - floor(v)), preservando a media exatamente."""
    piso = np.floor(valores)
    fracao = valores - piso
    return piso + (rng.random(valores.shape) < fracao)


def n_linhas(matriz: sp.csr_matrix) -> int:
    """Linhas de um CSR lidas do indptr: `.shape` do scipy e tipado como opcional."""
    return int(np.asarray(matriz.indptr).size) - 1


def emitir_bloco(
    controle: sp.csr_matrix,
    escala: np.ndarray,
    rng: np.random.Generator,
    n_celulas: int,
) -> sp.csr_matrix:
    """400 celulas de contagem crua: reamostra o pool, escala por gene, arredonda."""
    fonte = rng.integers(0, n_linhas(controle), size=n_celulas)
    sorteadas = controle[fonte]
    # A escala fica em float64 e multiplica uma vez: casta-la para float32 antes
    # arredondaria duas vezes (medido ~1,7x o erro no baseline.py da Arc).
    produto = np.asarray(sorteadas.data, dtype=np.float64) * escala[sorteadas.indices]
    contagens = arredondar_estocastico(produto, rng).astype(np.float32)
    # CSR reconstruido em vez de `.data` mutado: o array sorteado e uma view do
    # gather e mutar no lugar amarraria o bloco emitido ao pool de controle.
    bloco = sp.csr_matrix(
        (contagens, np.asarray(sorteadas.indices), np.asarray(sorteadas.indptr)),
        shape=(n_celulas, escala.size),
    )
    bloco.eliminate_zeros()
    bloco.sort_indices()
    return bloco


@dataclass
class DestinoCSR:
    """Os tres datasets do grupo `/X`, abertos uma vez e passados adiante.

    Existe porque `h5py.File.__getitem__` devolve `Group | Dataset | Datatype`:
    reindexar `arquivo["X"]["data"]` a cada bloco significaria reestreitar o tipo
    300 vezes, e a alternativa (ignorar o tipo) e o que o lint recusa.
    """

    dados: h5py.Dataset
    indices: h5py.Dataset
    indptr: h5py.Dataset


def abrir_destino(arquivo: h5py.File) -> DestinoCSR:
    grupo = cast(h5py.Group, arquivo["X"])
    return DestinoCSR(
        dados=cast(h5py.Dataset, grupo["data"]),
        indices=cast(h5py.Dataset, grupo["indices"]),
        indptr=cast(h5py.Dataset, grupo["indptr"]),
    )


def preparar_arquivo(
    caminho: Path, obs: pd.DataFrame, genes: list[str], comprimir: bool
) -> None:
    """Escreve o h5ad com obs e var completos e um `/X` vazio e redimensionavel.

    obs e var saem pelo proprio anndata — assim a codificacao (categorical, indice
    de string) e exatamente a que a versao instalada le de volta, e nao uma imitacao
    a mao que o `vcc prep` recusaria. Com `X=None` o anndata nao escreve `/X`
    nenhum, e o grupo CSR e criado aqui.
    """
    var = pd.DataFrame(index=pd.Index(genes, dtype=object))
    vazio = ad.AnnData(obs=obs, var=var)
    vazio.write_h5ad(caminho, compression="gzip" if comprimir else None)

    opcoes: dict[str, object] = {"chunks": (CHUNK_HDF5,)}
    if comprimir:
        opcoes["compression"] = "gzip"
        opcoes["compression_opts"] = NIVEL_GZIP
    with h5py.File(caminho, "r+") as arquivo:
        if "X" in arquivo:
            del arquivo["X"]
        grupo = arquivo.create_group("X")
        grupo.attrs["encoding-type"] = "csr_matrix"
        grupo.attrs["encoding-version"] = "0.1.0"
        grupo.attrs["shape"] = np.array([len(obs), len(genes)], dtype=np.int64)
        grupo.create_dataset(
            "data", shape=(0,), maxshape=(None,), dtype=np.float32, **opcoes
        )
        grupo.create_dataset(
            "indices", shape=(0,), maxshape=(None,), dtype=np.int32, **opcoes
        )
        # indptr em int64: a submissao completa chega a ~2,0e9 entradas, perto do
        # teto de int32 (2,147e9). O `validate_matrix_structure` da CLI casta para
        # int64 de todo modo.
        indptr = grupo.create_dataset("indptr", shape=(len(obs) + 1,), dtype=np.int64)
        indptr[0] = 0


def apendar_bloco(
    destino: DestinoCSR, bloco: sp.csr_matrix, linha0: int, nnz0: int
) -> int:
    """Apenda um bloco CSR ao `/X` aberto e devolve o novo total de nao nulos."""
    n_novos = int(bloco.nnz)
    destino.dados.resize((nnz0 + n_novos,))
    destino.indices.resize((nnz0 + n_novos,))
    if n_novos:
        destino.dados[nnz0:] = bloco.data
        destino.indices[nnz0:] = np.asarray(bloco.indices, dtype=np.int32)
    altura = n_linhas(bloco)
    # indptr[0] ja vale 0; cada bloco grava so as suas fronteiras de fim de linha.
    destino.indptr[linha0 + 1 : linha0 + altura + 1] = (
        np.asarray(bloco.indptr, dtype=np.int64)[1:] + nnz0
    )
    return nnz0 + n_novos


def carregar_controle(contexto: str, n_genes: int) -> sp.csr_matrix:
    """Pool de controle de um contexto, em CSR canonico de contagem inteira."""
    adata = ad.read_h5ad(VCC / f"context_{contexto}.h5ad")
    if int(adata.n_vars) != n_genes:
        raise ValueError(f"contexto {contexto}: {adata.n_vars} genes, esperado {n_genes}")
    matriz = sp.csr_matrix(adata.X)
    matriz.eliminate_zeros()
    matriz.sum_duplicates()
    if matriz.data.size and not np.array_equal(matriz.data, np.round(matriz.data)):
        raise ValueError(f"contexto {contexto}: o pool de controle nao e contagem inteira")
    return matriz


def montar_obs(contextos: tuple[str, ...], alvos: list[str]) -> pd.DataFrame:
    """obs da submissao: 400 linhas por (contexto, alvo), na ordem de emissao."""
    nomes: list[str] = []
    coluna_alvo: list[str] = []
    coluna_contexto: list[str] = []
    for contexto in contextos:
        for alvo in alvos:
            nomes.extend(f"{contexto}_{alvo}_{i:03d}" for i in range(CELULAS_POR_ALVO))
            coluna_alvo.extend([alvo] * CELULAS_POR_ALVO)
            coluna_contexto.extend([contexto] * CELULAS_POR_ALVO)
    return pd.DataFrame(
        {
            COLUNA_ALVO: pd.Categorical(coluna_alvo, categories=alvos),
            COLUNA_CONTEXTO: pd.Categorical(coluna_contexto, categories=list(contextos)),
        },
        index=pd.Index(nomes, dtype=object),
    )


def emitir(
    contextos: tuple[str, ...],
    alvos: list[str],
    genes: list[str],
    deltas: dict[str, np.ndarray] | None,
    delta_e_log_razao: bool,
    gate: dict[tuple[str, str], bool] | None,
    semente: int,
    caminho: Path,
    comprimir: bool,
) -> Medicao:
    """Loop principal: um bloco de 400 celulas por par, apendado em streaming."""
    n_genes = len(genes)
    n_celulas = len(contextos) * len(alvos) * CELULAS_POR_ALVO
    if n_celulas > MAX_CELULAS:
        raise ValueError(
            f"{n_celulas} celulas excede o teto de {MAX_CELULAS} do `vcc prep`"
        )

    with np.load(CTRL_BULK, allow_pickle=True) as bulk:
        b_ctrl = {c: _por_gene(bulk[f"b_ctrl_{c}"], n_genes, f"b_ctrl_{c}") for c in contextos}

    medicao = Medicao(celulas=n_celulas)
    preparar_arquivo(caminho, montar_obs(contextos, alvos), genes, comprimir)

    linha = 0
    nnz = 0
    with h5py.File(caminho, "r+") as arquivo:
        destino = abrir_destino(arquivo)
        for idx_ctx, contexto in enumerate(contextos):
            inicio_ctx = time.perf_counter()
            controle = carregar_controle(contexto, n_genes)
            pool = n_linhas(controle)
            print(
                f"{contexto}: pool com {pool} celulas, "
                f"{controle.nnz} nao nulos ({controle.nnz / pool:.0f}/celula)"
            )
            for idx_alvo, alvo in enumerate(alvos):
                if deltas is None:
                    delta = np.zeros(n_genes, dtype=np.float64)
                else:
                    delta = deltas[contexto][idx_alvo]
                if gate is not None and not gate[(contexto, alvo)]:
                    # Alvo sem transcrito: a predicao e o controle, nao um efeito.
                    delta = np.zeros(n_genes, dtype=np.float64)
                    medicao.pares_gateados.append(f"{contexto}/{alvo}")

                escala, inalcancavel, escala_max = escala_de_gene(
                    b_ctrl[contexto], delta, delta_e_log_razao=delta_e_log_razao
                )
                medicao.massa_inalcancavel.append(inalcancavel)
                medicao.escala_maxima = max(medicao.escala_maxima, escala_max)

                # Semente derivada do par, nao do contador do laco: o bloco de
                # (D, ABCC1) e o mesmo com --contextos D ou --contextos D,E,F, o
                # que torna a emissao parcial comparavel com a completa.
                rng = np.random.default_rng([semente, idx_ctx, idx_alvo])
                bloco = emitir_bloco(controle, escala, rng, CELULAS_POR_ALVO)

                totais = np.asarray(bloco.sum(axis=1), dtype=np.float64).ravel()
                pior = float(totais.max()) if totais.size else 0.0
                if pior > MAX_CONTAGEM_POR_CELULA:
                    raise ValueError(
                        f"{contexto}/{alvo}: celula com {pior:.0f} contagens, acima do "
                        f"teto de {MAX_CONTAGEM_POR_CELULA} do scorer"
                    )
                vazias = int((totais <= 0.0).sum())
                if vazias:
                    raise ValueError(
                        f"{contexto}/{alvo}: {vazias} celula(s) zeradas pelo "
                        f"arredondamento (escala maxima {escala_max:.3g}). O scorer "
                        "divide pelo tamanho de biblioteca; uma celula nula o quebra"
                    )
                medicao.max_contagem_celula = max(medicao.max_contagem_celula, pior)
                medicao.min_contagem_celula = min(
                    medicao.min_contagem_celula, float(totais.min())
                )
                medicao.nnz_por_bloco.append(int(bloco.nnz))

                nnz = apendar_bloco(destino, bloco, linha, nnz)
                linha += CELULAS_POR_ALVO
                if (idx_alvo + 1) % 50 == 0:
                    print(
                        f"  {contexto}: {idx_alvo + 1}/{len(alvos)} alvos, "
                        f"{nnz:,} nao nulos, {time.perf_counter() - inicio_ctx:.1f} s"
                    )
            del controle
            print(f"{contexto}: pronto em {time.perf_counter() - inicio_ctx:.1f} s")
        destino.indptr[-1] = nnz

    medicao.nnz = nnz
    if nnz > TETO_NNZ:
        raise ValueError(f"{nnz} nao nulos excede o teto de densidade {TETO_NNZ}")
    if linha != n_celulas:
        raise ValueError(f"escritas {linha} linhas, esperado {n_celulas}")
    return medicao


def _sha256(caminho: Path) -> str:
    """Impressao digital de um arquivo, lida em blocos para nao carregar 54 MB."""
    digestor = hashlib.sha256()
    with caminho.open("rb") as fonte:
        for bloco in iter(lambda: fonte.read(1 << 20), b""):
            digestor.update(bloco)
    return digestor.hexdigest()


def main() -> None:
    analisador = argparse.ArgumentParser(description=__doc__ or "", add_help=True)
    analisador.add_argument(
        "--delta",
        type=Path,
        default=None,
        help="npz do transfer_fit com delta_<ctx> (ou lfc_<ctx>) por (contexto, alvo)",
    )
    analisador.add_argument(
        "--stub",
        action="store_true",
        help="emite com delta zero, para testar escritor e validador sem o modelo",
    )
    analisador.add_argument(
        "--contextos",
        default=",".join(CONTEXTOS_FINAIS),
        help="contextos a emitir, separados por virgula (default D,E,F)",
    )
    analisador.add_argument("--semente", type=int, default=0, help="semente raiz (default 0)")
    analisador.add_argument(
        "--nome", default=None, help="nome base da saida (default 'stub' ou 'pred')"
    )
    analisador.add_argument(
        "--sem-gate",
        action="store_true",
        help="nao zera o delta nos pares com alvo sem transcrito",
    )
    analisador.add_argument(
        "--sem-compressao",
        action="store_true",
        help="grava o h5ad sem gzip: ~4x maior, mas bem mais rapido",
    )
    args = analisador.parse_args()

    if args.stub == (args.delta is not None):
        analisador.error("passe exatamente um de --stub ou --delta")
    contextos = tuple(c.strip() for c in str(args.contextos).split(",") if c.strip())
    desconhecidos = [c for c in contextos if c not in CONTEXTOS_FINAIS]
    if desconhecidos:
        analisador.error(
            f"contexto invalido {desconhecidos}; a fase final pontua {CONTEXTOS_FINAIS}"
        )

    SAIDA.mkdir(parents=True, exist_ok=True)
    nome = str(args.nome) if args.nome else ("stub" if args.stub else "pred")
    caminho = SAIDA / f"{nome}.h5ad"
    genes, alvos = ler_eixos()
    gate = None if args.sem_gate else ler_gate(alvos, contextos)

    deltas: dict[str, np.ndarray] | None = None
    delta_e_log_razao = False
    if args.delta is not None:
        deltas, delta_e_log_razao = carregar_delta(
            Path(args.delta), genes, alvos, contextos
        )

    print(
        f"emitindo {len(contextos)} contexto(s) x {len(alvos)} alvos x "
        f"{CELULAS_POR_ALVO} celulas = {len(contextos) * len(alvos) * CELULAS_POR_ALVO} "
        f"celulas em {caminho}"
    )
    processo = psutil.Process()
    inicio = time.perf_counter()
    medicao = emitir(
        contextos=contextos,
        alvos=alvos,
        genes=genes,
        deltas=deltas,
        delta_e_log_razao=delta_e_log_razao,
        gate=gate,
        semente=int(args.semente),
        caminho=caminho,
        comprimir=not args.sem_compressao,
    )
    segundos = time.perf_counter() - inicio
    memoria = processo.memory_info()
    pico_bytes = int(getattr(memoria, "peak_wset", memoria.rss))
    tamanho = caminho.stat().st_size
    nnz_bloco = np.asarray(medicao.nnz_por_bloco, dtype=np.float64)

    relatorio = {
        "arquivo": str(caminho),
        "modo": "stub" if args.stub else "delta",
        "delta": None if args.delta is None else str(args.delta),
        # Impressao digital do delta, para o `validar_numerico` checar IDENTIDADE
        # em vez de so regressao. Medido: duas submissoes que diferem apenas na
        # parcela da media tem correlacao +0,996, e a inclinacao do validador
        # entre elas da 0,996 — dentro da janela [0,85; 1,15], logo ele NAO as
        # distingue. Nenhuma janela estreita o suficiente para separar as duas
        # sobreviveria ao ruido de amostragem de 400 celulas, entao o caminho e
        # identidade de conteudo, nao aperto de limite.
        "delta_sha256": None if args.delta is None else _sha256(Path(args.delta)),
        "espaco_do_delta": "log_razao_de_proporcao" if delta_e_log_razao else "log1p_cp50k",
        "contextos": list(contextos),
        "alvos": len(alvos),
        "genes": len(genes),
        "celulas_por_alvo": CELULAS_POR_ALVO,
        "celulas": medicao.celulas,
        "semente_raiz": int(args.semente),
        "semente_por_bloco": "np.random.default_rng([semente, idx_contexto, idx_alvo])",
        "gate_aplicado": not args.sem_gate,
        "pares_gateados": medicao.pares_gateados,
        "compressao": "nenhuma" if args.sem_compressao else f"gzip-{NIVEL_GZIP}",
        "tempo_s": round(segundos, 2),
        "tempo_por_bloco_ms": round(1000.0 * segundos / max(len(nnz_bloco), 1), 2),
        "pico_ram_mib": round(pico_bytes / 2**20, 1),
        "tamanho_bytes": tamanho,
        "tamanho_mib": round(tamanho / 2**20, 1),
        "nnz_total": medicao.nnz,
        "nnz_por_celula_medio": round(medicao.nnz / medicao.celulas, 1),
        "nnz_por_bloco_min": int(nnz_bloco.min()) if nnz_bloco.size else 0,
        "nnz_por_bloco_max": int(nnz_bloco.max()) if nnz_bloco.size else 0,
        "teto_nnz": TETO_NNZ,
        "fracao_do_teto_nnz": round(medicao.nnz / TETO_NNZ, 4),
        "nnz_projetado_3_contextos": int(
            round(medicao.nnz * len(CONTEXTOS_FINAIS) / len(contextos))
        ),
        "fracao_do_teto_projetada": round(
            medicao.nnz * len(CONTEXTOS_FINAIS) / len(contextos) / TETO_NNZ, 4
        ),
        "contagem_por_celula_max": medicao.max_contagem_celula,
        "contagem_por_celula_min": medicao.min_contagem_celula,
        "teto_contagem_por_celula": MAX_CONTAGEM_POR_CELULA,
        "escala_de_gene_maxima": medicao.escala_maxima,
        "massa_inalcancavel_media": float(np.mean(medicao.massa_inalcancavel)),
        "massa_inalcancavel_max": float(np.max(medicao.massa_inalcancavel)),
        "arredondamento": "estocastico (floor + Bernoulli da fracao)",
    }
    (SAIDA / f"medicao_{nome}.json").write_text(
        json.dumps(relatorio, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(
        f"\n{medicao.celulas} celulas, {medicao.nnz:,} nao nulos "
        f"({relatorio['nnz_por_celula_medio']}/celula, "
        f"{relatorio['fracao_do_teto_nnz']:.2%} do teto), "
        f"{relatorio['tamanho_mib']} MiB, {segundos:.1f} s, "
        f"pico {relatorio['pico_ram_mib']} MiB"
    )
    print(f"medicao em {SAIDA / f'medicao_{nome}.json'}")


if __name__ == "__main__":
    main()
