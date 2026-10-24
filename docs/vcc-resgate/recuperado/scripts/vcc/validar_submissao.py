# scripts/vcc/validar_submissao.py
"""
Etapa 6 de docs/VCC-2026-ROTA.md: o portao entre nos e uma submissao queimada.
Sao 2 por dia e so a ultima conta, entao cada item do contrato e testado aqui,
contra o arquivo gerado, antes de qualquer envio.

OS ITENS SAO OS DO SCORER, NAO OS NOSSOS. Cada checagem abaixo espelha uma
rejeicao concreta de `vcc/prep.py` da CLI 0.3.0 (ou de `cell_eval2` atraves
dela), com o numero que ela usa:

  estrutura     `validate_matrix_structure` — indptr monotonico, do tamanho certo,
                comecando em zero, com fim dentro de `data`, e indices em faixa.
                Um CSR incoerente nao e pego pelas checagens de valor: `.data`
                varre bem e a corrupcao so morde quando o scipy anda os ponteiros
                de linha. No servidor isso e SIGSEGV no meio do load, e a submissao
                fica em `scoring` sem erro nenhum.
  canonico      `canonicalize_matrix` — sem folga no fim de `data`/`indices` e sem
                indice de coluna repetido na linha. O escritor de arquivo do
                scorer (`require_unique_indices=True`) levanta nos repetidos.
  valores       `validate_matrix_values` + `_require_counts` — finito, nao
                negativo, inteiro, soma por celula <= 1e6, nenhuma perturbacao
                inteiramente nula.
  rotulos       `validate_targets` — por contexto, exatamente os 300 alvos
                oficiais, exatamente 400 celulas cada; nenhuma linha de controle.
  eixos         `run_prep` — var_names identico a gene_names.csv na ordem, sem
                repetido, 18.533 colunas.
  tamanho       `MAX_CELL_DIM` 400.000 celulas e `sizing.HARD_CAP_NNZ` 4,75e9
                entradas armazenadas.

Mais um item que o contrato nao escreve e o scorer precisa: nenhuma celula com
soma zero. O `prep` so recusa perturbacao inteiramente nula, mas a normalizacao
por CPM divide pelo tamanho de biblioteca da celula, e uma celula zerada pelo
arredondamento entra ali como divisao por zero. E reportado como item proprio
para que a distincao fique visivel.

TUDO EM STREAMING. O arquivo completo tem ~2e9 entradas; ler `X` inteiro seria
~16 GB. As varreduras de valor andam em blocos de linha, com as fronteiras lidas
do proprio `indptr`, entao cada item e medido no arquivo todo a custo de RAM
constante.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/validar_submissao.py runs/vcc/submissao/stub.h5ad --contextos D
    ./.venv/Scripts/python.exe -I scripts/vcc/validar_submissao.py runs/vcc/submissao/pred.h5ad

SEGUNDA OPINIAO, que e a autoridade. O proprio `vcc prep --dry-run` valida sem
escrever nada, sem rede (com `--contexts` explicito nao consulta a fase) e sem
gastar cota de submissao:

    vcc prep runs/vcc/submissao/pred.h5ad \
      -g data/raw/vcc/controls-final/gene_names.csv \
      --perts data/raw/vcc/controls-final/pert_counts.csv \
      --contexts D,E,F --dry-run

Rode os dois. A diferenca e de proposito: o `prep` para no primeiro erro, entao
cada defeito custa uma ida e volta; este script mede todos os itens de uma vez e
nomeia cada um que quebrou. Os dois foram conferidos um contra o outro em
2026-10-22, nos dois sentidos — aprovacao e recusa.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import h5py
import numpy as np
import pandas as pd

VCC = Path("data/raw/vcc/controls-final")

CELULAS_POR_ALVO = 400
N_GENES = 18533
MAX_CONTAGEM_POR_CELULA = 1_000_000
MAX_CELULAS = 400_000
TETO_NNZ = 4_750_000_000
ROTULO_CONTROLE = "non-targeting"
COLUNA_ALVO = "target_gene"
COLUNA_CONTEXTO = "context"
CONTEXTOS_FINAIS = ("D", "E", "F")

# Linhas por varredura. 2048 x ~5.700 nao nulos = ~12 Mi entradas, ~46 MB entre
# `data` (float32) e `indices` (int32) — grande o bastante para o custo por bloco
# desaparecer e pequeno o bastante para o pico nao depender do tamanho do arquivo.
BLOCO_LINHAS = 2048


@dataclass
class Item:
    """Um item do contrato, com o veredito e o numero que o sustenta."""

    nome: str
    passou: bool
    detalhe: str


class Laudo:
    """Acumula os itens e decide o veredito final."""

    def __init__(self) -> None:
        self.itens: list[Item] = []

    def registrar(self, nome: str, passou: bool, detalhe: str) -> bool:
        self.itens.append(Item(nome, passou, detalhe))
        marca = "ok  " if passou else "FALHA"
        print(f"  [{marca}] {nome}: {detalhe}")
        return passou

    @property
    def falhas(self) -> list[Item]:
        return [i for i in self.itens if not i.passou]


def _texto(valores: np.ndarray) -> np.ndarray:
    """Array de string vindo do HDF5, que pode chegar como bytes ou como objeto."""
    plano = np.asarray(valores).ravel()
    if plano.size and isinstance(plano[0], bytes):
        return np.array([v.decode("utf-8") for v in plano], dtype=object)
    return plano.astype(object)


def coluna_de_texto(grupo: h5py.Group, nome: str) -> np.ndarray:
    """Le uma coluna de obs/var como texto, qualquer que seja a codificacao.

    O anndata escreve a mesma informacao de quatro formas (categorical, string-array,
    nullable-string-array, array simples) e qual delas sai depende da versao que
    gravou. O validador tem de ler o arquivo que recebeu, nao o que teria escrito.
    """
    objeto = grupo[nome]
    if isinstance(objeto, h5py.Group):
        tipo = objeto.attrs.get("encoding-type", "")
        if tipo == "categorical":
            categorias = _texto(cast(h5py.Dataset, objeto["categories"])[:])
            codigos = np.asarray(cast(h5py.Dataset, objeto["codes"])[:], dtype=np.int64)
            if codigos.size and (codigos.min() < 0 or codigos.max() >= categorias.size):
                raise ValueError(f"obs['{nome}']: codigo categorico fora das categorias")
            return categorias[codigos]
        if "values" in objeto:
            return _texto(cast(h5py.Dataset, objeto["values"])[:])
        raise ValueError(f"obs['{nome}']: grupo de codificacao '{tipo}' nao reconhecido")
    return _texto(cast(h5py.Dataset, objeto)[:])


def ler_eixos_oficiais() -> tuple[list[str], list[str]]:
    genes = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
    alvos = pd.read_csv(VCC / "pert_counts.csv")["target_gene"].astype(str).tolist()
    return genes, alvos


@dataclass
class Varredura:
    """O que uma passada em streaming por `X` mediu."""

    nnz: int
    n_nao_finitos: int
    n_negativos: int
    n_fracionarios: int
    n_zeros_explicitos: int
    n_indices_fora: int
    n_pares_desordenados: int
    totais_por_celula: np.ndarray
    nnz_por_celula: np.ndarray
    valor_maximo: float
    primeiro_defeito: str


def varrer(
    dados: h5py.Dataset,
    indices: h5py.Dataset,
    indptr: np.ndarray,
    n_obs: int,
    n_vars: int,
) -> Varredura:
    """Uma passada por `X` em blocos de linha, medindo todos os itens de valor."""
    totais = np.zeros(n_obs, dtype=np.float64)
    por_celula = np.zeros(n_obs, dtype=np.int64)
    nao_finitos = negativos = fracionarios = zeros = fora = desordenados = 0
    maximo = 0.0
    defeito = ""

    for inicio in range(0, n_obs, BLOCO_LINHAS):
        fim = min(inicio + BLOCO_LINHAS, n_obs)
        fronteiras = indptr[inicio : fim + 1]
        a, b = int(fronteiras[0]), int(fronteiras[-1])
        comprimentos = np.diff(fronteiras)
        por_celula[inicio:fim] = comprimentos
        if b <= a:
            continue
        valores = np.asarray(dados[a:b], dtype=np.float64)
        colunas = np.asarray(indices[a:b], dtype=np.int64)

        finito = np.isfinite(valores)
        n_bloco_nao_finitos = int((~finito).sum())
        nao_finitos += n_bloco_nao_finitos
        # As comparacoes abaixo usam so os valores finitos: NaN responde False a
        # tudo, e contar um NaN como "nao negativo" diria que o item passou.
        bons = valores[finito]
        negativos += int((bons < 0.0).sum())
        fracionarios += int((bons != np.trunc(bons)).sum())
        zeros += int((bons == 0.0).sum())
        if bons.size:
            maximo = max(maximo, float(bons.max()))

        fora += int(((colunas < 0) | (colunas >= n_vars)).sum())

        # Ordem estrita dentro da linha: pega desordem e indice repetido de uma
        # vez. A comparacao so vale entre entradas da mesma linha, por isso a
        # mascara de fronteira.
        if colunas.size > 1:
            linha_por_entrada = np.repeat(np.arange(fim - inicio), comprimentos)
            mesma_linha = linha_por_entrada[1:] == linha_por_entrada[:-1]
            crescente = colunas[1:] > colunas[:-1]
            desordenados += int((mesma_linha & ~crescente).sum())
        else:
            linha_por_entrada = np.repeat(np.arange(fim - inicio), comprimentos)

        seguros = np.where(finito, valores, 0.0)
        totais[inicio:fim] = np.bincount(
            linha_por_entrada, weights=seguros, minlength=fim - inicio
        )
        if not defeito and (n_bloco_nao_finitos or negativos or fracionarios):
            # A primeira linha com defeito, para que a mensagem de falha diga onde.
            ruim = ~finito | (valores < 0.0) | (valores != np.trunc(valores))
            if ruim.any():
                primeira = int(linha_por_entrada[int(np.argmax(ruim))]) + inicio
                defeito = f"primeira celula afetada: linha {primeira}"

    return Varredura(
        nnz=int(indptr[-1]),
        n_nao_finitos=nao_finitos,
        n_negativos=negativos,
        n_fracionarios=fracionarios,
        n_zeros_explicitos=zeros,
        n_indices_fora=fora,
        n_pares_desordenados=desordenados,
        totais_por_celula=totais,
        nnz_por_celula=por_celula,
        valor_maximo=maximo,
        primeiro_defeito=defeito,
    )


def checar_estrutura(
    laudo: Laudo, indptr: np.ndarray, n_dados: int, n_indices: int, n_obs: int
) -> bool:
    """Os itens de `validate_matrix_structure` e `canonicalize_matrix`."""
    ok = laudo.registrar(
        "estrutura_indptr_tamanho",
        indptr.size == n_obs + 1,
        f"{indptr.size} entradas, esperado {n_obs + 1}",
    )
    if not ok:
        # Sem indptr do tamanho certo nada a jusante e interpretavel.
        return False
    diferencas = np.diff(indptr)
    ok &= laudo.registrar(
        "estrutura_indptr_inicio", int(indptr[0]) == 0, f"indptr[0] = {int(indptr[0])}"
    )
    n_negativas = int((diferencas < 0).sum())
    ok &= laudo.registrar(
        "estrutura_indptr_monotonico",
        n_negativas == 0,
        f"{n_negativas} linha(s) com comprimento negativo",
    )
    ok &= laudo.registrar(
        "estrutura_pares",
        n_dados == n_indices,
        f"data {n_dados:,} vs indices {n_indices:,}",
    )
    fim = int(indptr[-1])
    ok &= laudo.registrar(
        "estrutura_fim_dentro",
        fim <= n_dados,
        f"indptr[-1] = {fim:,}, data = {n_dados:,}",
    )
    # Folga no fim e legal para o scipy e proibida aqui: o `prep` a corta, mas o
    # servidor dimensiona a submissao por `data.shape[0]` enquanto o `prep` julgou
    # por `nnz`, e os dois portoes de tamanho passam a discordar.
    ok &= laudo.registrar(
        "canonico_sem_folga",
        fim == n_dados,
        f"{n_dados - fim:,} entrada(s) de folga depois de indptr[-1]",
    )
    return ok


def checar_rotulos(
    laudo: Laudo,
    alvo_por_celula: np.ndarray,
    contexto_por_celula: np.ndarray,
    alvos_oficiais: list[str],
    contextos_esperados: tuple[str, ...],
) -> None:
    """Os itens de `validate_targets` e a recusa de controle de `run_prep`."""
    n_controle = int((alvo_por_celula == ROTULO_CONTROLE).sum())
    laudo.registrar(
        "sem_linha_de_controle",
        n_controle == 0,
        f"{n_controle} celula(s) rotuladas '{ROTULO_CONTROLE}'",
    )

    presentes = sorted(set(contexto_por_celula.tolist()))
    laudo.registrar(
        "contextos",
        presentes == sorted(contextos_esperados),
        f"presentes {presentes}, esperados {sorted(contextos_esperados)}",
    )

    esperados = set(alvos_oficiais)
    erros_alvo: list[str] = []
    erros_contagem: list[str] = []
    for contexto in contextos_esperados:
        no_contexto = contexto_por_celula == contexto
        if not no_contexto.any():
            erros_alvo.append(f"{contexto}: nenhuma celula")
            continue
        rotulos = alvo_por_celula[no_contexto]
        observados = set(rotulos.tolist()) - {ROTULO_CONTROLE}
        faltando = sorted(esperados - observados)
        extras = sorted(observados - esperados)
        if faltando or extras:
            erros_alvo.append(
                f"{contexto}: {len(faltando)} faltando ({faltando[:3]}), "
                f"{len(extras)} extras ({extras[:3]})"
            )
        valores, contagens = np.unique(rotulos, return_counts=True)
        erradas = [
            f"{v}={c}"
            for v, c in zip(valores.tolist(), contagens.tolist(), strict=True)
            if v in esperados and c != CELULAS_POR_ALVO
        ]
        if erradas:
            erros_contagem.append(f"{contexto}: {len(erradas)} alvo(s) {erradas[:3]}")

    laudo.registrar(
        "alvos_por_contexto",
        not erros_alvo,
        "; ".join(erros_alvo) if erros_alvo else f"{len(esperados)} alvos oficiais em cada contexto",
    )
    laudo.registrar(
        "celulas_por_alvo",
        not erros_contagem,
        "; ".join(erros_contagem)
        if erros_contagem
        else f"exatamente {CELULAS_POR_ALVO} por (contexto, alvo)",
    )


def validar(caminho: Path, contextos_esperados: tuple[str, ...]) -> Laudo:
    laudo = Laudo()
    genes_oficiais, alvos_oficiais = ler_eixos_oficiais()

    with h5py.File(caminho, "r") as arquivo:
        if "X" not in arquivo or not isinstance(arquivo["X"], h5py.Group):
            laudo.registrar("matriz_csr", False, "/X ausente ou nao e grupo CSR/CSC")
            return laudo
        grupo = cast(h5py.Group, arquivo["X"])
        tipo = str(grupo.attrs.get("encoding-type", ""))
        laudo.registrar("matriz_csr", tipo == "csr_matrix", f"encoding-type = '{tipo}'")
        forma = np.asarray(grupo.attrs["shape"], dtype=np.int64)
        n_obs, n_vars = int(forma[0]), int(forma[1])

        obs = cast(h5py.Group, arquivo["obs"])
        var = cast(h5py.Group, arquivo["var"])
        indice_var = str(var.attrs.get("_index", "_index"))
        genes = coluna_de_texto(var, indice_var)

        colunas_obs = [str(c) for c in np.asarray(obs.attrs.get("column-order", [])).ravel()]
        tem_colunas = COLUNA_ALVO in colunas_obs and COLUNA_CONTEXTO in colunas_obs
        laudo.registrar(
            "colunas_obs",
            tem_colunas,
            f"colunas {colunas_obs}; exigidas '{COLUNA_ALVO}' e '{COLUNA_CONTEXTO}'",
        )
        if not tem_colunas:
            return laudo
        alvo_por_celula = coluna_de_texto(obs, COLUNA_ALVO)
        contexto_por_celula = coluna_de_texto(obs, COLUNA_CONTEXTO)

        laudo.registrar(
            "forma_coerente",
            alvo_por_celula.size == n_obs and genes.size == n_vars,
            f"X {n_obs} x {n_vars}, obs {alvo_por_celula.size}, var {genes.size}",
        )
        laudo.registrar(
            "dimensao_de_genes", n_vars == N_GENES, f"{n_vars} colunas, exigidas {N_GENES}"
        )
        laudo.registrar(
            "eixo_de_genes",
            genes.tolist() == genes_oficiais,
            "identico a gene_names.csv na ordem"
            if genes.tolist() == genes_oficiais
            else f"divergente (primeira diferenca em {_primeira_diferenca(genes, genes_oficiais)})",
        )
        repetidos = genes.size - len(set(genes.tolist()))
        laudo.registrar(
            "genes_sem_repetido", repetidos == 0, f"{repetidos} simbolo(s) repetido(s)"
        )

        n_celulas_esperadas = len(contextos_esperados) * len(alvos_oficiais) * CELULAS_POR_ALVO
        laudo.registrar(
            "total_de_celulas",
            n_obs == n_celulas_esperadas,
            f"{n_obs:,} celulas, esperado {n_celulas_esperadas:,}",
        )
        laudo.registrar(
            "teto_de_celulas",
            n_obs <= MAX_CELULAS,
            f"{n_obs:,} celulas, teto {MAX_CELULAS:,}",
        )

        checar_rotulos(
            laudo, alvo_por_celula, contexto_por_celula, alvos_oficiais, contextos_esperados
        )

        dados = cast(h5py.Dataset, grupo["data"])
        indices = cast(h5py.Dataset, grupo["indices"])
        indptr = np.asarray(cast(h5py.Dataset, grupo["indptr"])[:], dtype=np.int64)
        if not checar_estrutura(
            laudo, indptr, int(dados.shape[0]), int(indices.shape[0]), n_obs
        ):
            return laudo

        inicio = time.perf_counter()
        medida = varrer(dados, indices, indptr, n_obs, n_vars)
        print(f"  (varredura de {medida.nnz:,} entradas em {time.perf_counter() - inicio:.1f} s)")

    laudo.registrar(
        "estrutura_indices_em_faixa",
        medida.n_indices_fora == 0,
        f"{medida.n_indices_fora} indice(s) fora de 0..{n_vars - 1}",
    )
    laudo.registrar(
        "canonico_ordem_estrita",
        medida.n_pares_desordenados == 0,
        f"{medida.n_pares_desordenados} par(es) fora de ordem ou com indice repetido",
    )
    laudo.registrar(
        "sem_zeros_explicitos",
        medida.n_zeros_explicitos == 0,
        f"{medida.n_zeros_explicitos} zero(s) armazenado(s)",
    )
    laudo.registrar(
        "valores_finitos",
        medida.n_nao_finitos == 0,
        f"{medida.n_nao_finitos} NaN/inf. {medida.primeiro_defeito}".strip(),
    )
    laudo.registrar(
        "valores_nao_negativos",
        medida.n_negativos == 0,
        f"{medida.n_negativos} valor(es) negativo(s). {medida.primeiro_defeito}".strip(),
    )
    laudo.registrar(
        "valores_inteiros",
        medida.n_fracionarios == 0,
        f"{medida.n_fracionarios} valor(es) fracionario(s), maximo {medida.valor_maximo:g}. "
        f"{medida.primeiro_defeito}".strip(),
    )

    totais = medida.totais_por_celula
    pior = float(totais.max()) if totais.size else 0.0
    acima = int((totais > MAX_CONTAGEM_POR_CELULA).sum())
    laudo.registrar(
        "contagem_por_celula",
        acima == 0,
        f"{acima} celula(s) acima de {MAX_CONTAGEM_POR_CELULA:,}; maior soma {pior:,.0f}",
    )
    vazias = int((totais <= 0.0).sum())
    laudo.registrar(
        "celula_nao_nula",
        vazias == 0,
        f"{vazias} celula(s) com soma zero (a normalizacao por CPM divide por isso)",
    )

    nulas: list[str] = []
    for contexto in contextos_esperados:
        no_contexto = contexto_por_celula == contexto
        for alvo in alvos_oficiais:
            mascara = no_contexto & (alvo_por_celula == alvo)
            if mascara.any() and not totais[mascara].any():
                nulas.append(f"{contexto}/{alvo}")
    laudo.registrar(
        "perturbacao_nao_nula",
        not nulas,
        f"{len(nulas)} perturbacao(oes) inteiramente nula(s): {nulas[:5]}"
        if nulas
        else "nenhuma perturbacao inteiramente nula",
    )

    laudo.registrar(
        "teto_de_densidade",
        medida.nnz <= TETO_NNZ,
        f"{medida.nnz:,} entradas armazenadas, teto {TETO_NNZ:,} "
        f"({medida.nnz / TETO_NNZ:.2%}); {medida.nnz / max(n_obs, 1):.0f}/celula "
        f"(min {int(medida.nnz_por_celula.min())}, max {int(medida.nnz_por_celula.max())})",
    )
    return laudo


def _primeira_diferenca(observado: np.ndarray, esperado: list[str]) -> str:
    """Posicao e par da primeira divergencia entre dois eixos de genes."""
    lista = observado.tolist()
    for i, (a, b) in enumerate(zip(lista, esperado, strict=False)):
        if a != b:
            return f"posicao {i}: '{a}' vs '{b}'"
    return f"comprimento {len(lista)} vs {len(esperado)}"


def main() -> None:
    analisador = argparse.ArgumentParser(description=__doc__ or "", add_help=True)
    analisador.add_argument("arquivo", type=Path, help="h5ad gerado pelo emissor_400")
    analisador.add_argument(
        "--contextos",
        default=",".join(CONTEXTOS_FINAIS),
        help="contextos que o arquivo deve cobrir (default D,E,F)",
    )
    analisador.add_argument(
        "--laudo", type=Path, default=None, help="caminho do laudo JSON (default ao lado do h5ad)"
    )
    args = analisador.parse_args()

    caminho = Path(args.arquivo)
    if not caminho.exists():
        analisador.error(f"{caminho} nao existe")
    contextos = tuple(c.strip() for c in str(args.contextos).split(",") if c.strip())
    completo = sorted(contextos) == sorted(CONTEXTOS_FINAIS)

    print(f"validando {caminho} ({caminho.stat().st_size / 2**20:.1f} MiB)")
    print(f"contextos exigidos: {list(contextos)}" + ("" if completo else "  [PARCIAL]"))
    laudo = validar(caminho, contextos)

    destino = Path(args.laudo) if args.laudo else caminho.with_name(f"validacao_{caminho.stem}.json")
    destino.write_text(
        json.dumps(
            {
                "arquivo": str(caminho),
                "contextos_exigidos": list(contextos),
                "painel_completo": completo,
                "itens": [
                    {"nome": i.nome, "passou": i.passou, "detalhe": i.detalhe}
                    for i in laudo.itens
                ],
                "itens_testados": len(laudo.itens),
                "falhas": [i.nome for i in laudo.falhas],
                "veredito": "aprovado" if not laudo.falhas else "reprovado",
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    if laudo.falhas:
        print(f"\nREPROVADO: {len(laudo.falhas)} de {len(laudo.itens)} itens falharam")
        for item in laudo.falhas:
            print(f"  - {item.nome}: {item.detalhe}")
        print(f"laudo em {destino}")
        sys.exit(1)

    print(f"\nAPROVADO: {len(laudo.itens)} itens do contrato, todos passaram")
    if not completo:
        print(
            "Painel PARCIAL: a cobertura de contexto foi testada contra "
            f"{list(contextos)}, nao contra {list(CONTEXTOS_FINAIS)}. "
            "Uma submissao e um arquivo so, com os tres."
        )
    print(f"laudo em {destino}")


if __name__ == "__main__":
    main()
