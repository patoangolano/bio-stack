# scripts/vcc/avaliador_proxy.py
"""Monta o par (real, predito) que exercita a sequencia de pontuacao do cell-eval2.

Etapa 5 de docs/VCC-2026-ROTA.md. O desafio esconde os bundles de D/E/F, entao o
`avaliador_proxy` precisa de um dataset com a MESMA FORMA do oficial — contagens
cruas, 18 533 genes, 400 celulas por perturbacao, rotulo de controle
`non-targeting`, coluna `target_gene` — e com estrutura de DE de verdade, senao o
gate de CPM do preset nao admite perturbacao nenhuma e as seis metricas saem
degeneradas.

O que e real aqui:

* as celulas. Todas saem das 18 400 celulas de controle entregues em
  `data/raw/vcc/controls-final/context_{D,E,F}.h5ad` — contagens medidas, nao
  simuladas. Distribuicao de biblioteca, esparsidade (5 678 nnz/celula medido em
  D) e correlacao entre genes vem do dado.
* o efeito de perturbacao. As razoes de dobra vem de `fold_changes.npy` de um dos
  datasets preprocessados do PIE (ver `FONTES`), medidas pela Arc. Com
  `--fonte-efeitos orion --contexto-pie hct116` ha 148 dos 300 alvos do painel
  VCC com efeito medido; com `arc_vcc_25 / ARC_H1`, 12.

O que e proxy, declarado:

* o emparelhamento contexto-efeito. O efeito medido numa linhagem (H1
  embrionaria, HCT116 colorretal) e aplicado sobre celulas de D/E/F, que sao
  outro tipo celular. Isso nao e uma predicao biologica — e um gerador de dado
  com a estatistica certa para o avaliador morder. A consequencia para a leitura
  de qualquer pontuacao obtida contra este lado real esta em
  `scripts/vcc/avaliador_recorte.py`.
* o eixo de genes. Nenhuma fonte cobre os 18 533: `arc_vcc_25` cruza 12 286 e
  `orion` 14 690. O que nao cruza entra com razao de dobra 1 (sem efeito), assim
  como o que a Arc marcou `tested = False`.

O lado predito e o mesmo mecanismo com o efeito ENCOLHIDO por um expoente e com
jitter por gene: direcao majoritariamente certa, magnitude errada. E o regime em
que um modelo real cai, e o unico que produz as seis metricas fora dos extremos.
Dois casos de borda sao referencias uteis e saem dos mesmos parametros:
`--encolhimento 1.0 --jitter 0.0` da o QUASE-ORACULO (o efeito real, so com ruido
de amostragem independente) e `--encolhimento 0.0 --jitter 0.0` da o CONTROLE
PURO (razao de dobra 1 em todo gene, isto e, delta zero).
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

RAIZ = Path(__file__).resolve().parents[2]
CONTROLES = RAIZ / "data" / "raw" / "vcc" / "controls-final"
PIE = RAIZ / "data" / "raw" / "pie" / "arrays"

#: O rotulo e a coluna vem do `manifest.json` do painel, nao de convencao nossa.
COLUNA_PERT = "target_gene"
ROTULO_CONTROLE = "non-targeting"

#: De onde vem a razao de dobra medida, e qual contexto de cada fonte usar.
#: A escolha decide quantas perturbacoes do PAINEL o proxy consegue rotular, e a
#: diferenca e grande — medido em 2026-10-22:
#:
#:   arc_vcc_25 / ARC_H1   146 nocautes,  12 deles entre os 300 alvos do VCC
#:   orion      / hct116  8 987 linhas, 148 deles entre os 300 alvos do VCC
#:
#: Os 12 do arc_vcc_25 sao o mesmo 12 que o memorando da rota tirou do
#: `train.json` do PIE; os 148 do orion sao o mesmo 148. Para um proxy com
#: rotulos REAIS do painel e mais de uma duzia de perturbacoes, orion e a unica
#: fonte que serve.
FONTES: dict[str, tuple[str, str]] = {
    # `ARC_H1_VAL` existe mas e o fold de validacao, e nao e o default.
    "arc_vcc_25": ("PIE_arc_vcc_25", "ARC_H1"),
    "orion": ("PIE_x_atlas_orion", "hct116"),
    "jiang": ("PIE_jiang", ""),
    "replogle": ("PIE_replogle_nadig_essential", ""),
}


def _forma(m: sp.csr_matrix) -> tuple[int, int]:
    """Forma de uma CSR como par de inteiros.

    O scipy tipa `.shape` como opcional porque a classe base admite matriz sem
    forma definida; uma CSR construida sempre tem duas dimensoes.
    """
    return cast(tuple[int, int], m.shape)


@dataclass(frozen=True)
class Efeitos:
    """Razoes de dobra por (perturbacao, gene) no eixo de 18 533 genes."""

    alvos: list[str]
    no_painel_vcc: list[bool]
    fold_changes: np.ndarray  # (n_alvos, n_genes) float64, 1.0 = sem efeito
    n_genes_cobertos: int
    n_pares_efetivos: int
    fonte: str
    contexto_fonte: str


def _carregar_efeitos_pie(eixo_genes: np.ndarray, alvos_vcc: set[str], n_alvos: int,
                          *, min_genes_moventes: int, limiar_movente: float,
                          so_alvos_vcc: bool, fonte_efeitos: str,
                          contexto_pie: str | None) -> Efeitos:
    """Razoes de dobra medidas por alvo, projetadas no nosso eixo de genes.

    Escolhe os alvos com MAIS genes moventes primeiro, e entre eles da preferencia
    aos que estao no painel VCC. O gate do preset exige que a referencia tenha
    >= 10 genes significativos por perturbacao (`min_gate_size`); um alvo com
    efeito fraco atravessa o avaliador como NaN e sai do agregado, o que esvazia a
    verificacao em vez de testa-la.

    ⚠️ QUANTOS ALVOS DO PAINEL CADA FONTE SUSTENTA, medido em 2026-10-22:
    `arc_vcc_25 / ARC_H1` tem 146 nocautes e exatamente **12** entre os 300 alvos
    do VCC (17 contando `ARC_H1_VAL`, que e o fold de validacao);
    `orion / hct116` tem 8 987 linhas e **148**. Os dois numeros sao os mesmos que
    o memorando da rota tirou do `train.json` do PIE. Entao um proxy com rotulos
    REAIS do painel e mais de uma duzia de perturbacoes exige
    `--fonte-efeitos orion`. Sem `--so-alvos-vcc` a selecao completa com nocautes
    que nao sao alvos do painel, o que nao prejudica o avaliador — o que ele mede
    e a FORMA do dado (numero de perturbacoes, profundidade, estrutura de DE), nao
    a identidade dos genes — mas impede recortar uma submissao real contra este
    lado, porque os rotulos nao existem na submissao.
    """
    diretorio, ctx_default = FONTES[fonte_efeitos]
    base = PIE / diretorio / "preprocessed"
    meta = json.loads((base / "meta.json").read_text(encoding="utf-8"))
    genes_pie = np.asarray(meta["genes"], dtype=str)
    pert_para_id = meta["pert_to_id"]
    nome_ctx = contexto_pie or ctx_default
    if nome_ctx not in meta["context_to_id"]:
        raise SystemExit(
            f"contexto {nome_ctx!r} nao existe em {fonte_efeitos}; disponiveis: "
            f"{sorted(meta['context_to_id'])}"
        )
    ctx_id = meta["context_to_id"][nome_ctx]

    # mmap, nao leitura inteira: o `fold_changes.npy` do orion tem 1,18 GB e o
    # `tested.npy` 295 MB, e so as linhas dos candidatos sao lidas (148 de 8 987
    # no caso que importa). Ler os tres inteiros custaria 2,7 GB para usar 13 MB.
    fc = np.load(base / "fold_changes.npy", mmap_mode="r")
    tested = np.load(base / "tested.npy", mmap_mode="r")
    pert_ids = np.load(base / "pert_ids.npy")
    ctx_ids = np.load(base / "ctx_ids.npy")

    # Indice do nosso eixo por simbolo; genes do PIE fora dele sao descartados.
    posicao = {g: i for i, g in enumerate(eixo_genes)}
    col_pie = np.array([i for i, g in enumerate(genes_pie) if g in posicao], dtype=np.int64)
    col_nossa = np.array([posicao[genes_pie[i]] for i in col_pie], dtype=np.int64)

    id_para_pert = {v: k for k, v in pert_para_id.items()}
    candidatos: list[tuple[bool, int, str, int]] = []
    for linha in range(fc.shape[0]):
        if int(ctx_ids[linha]) != ctx_id:
            continue
        nome = id_para_pert.get(int(pert_ids[linha]))
        if nome is None:
            continue
        no_painel = nome in alvos_vcc
        if so_alvos_vcc and not no_painel:
            continue
        # "Movente" = testado, finito, positivo e a mais de `limiar_movente` em
        # |log2fc| de 1. O default 0,1 (~7%) nao e generoso: a referencia tem
        # 18 400 celulas contra 400 perturbadas, e nesse desequilibrio o Wilcoxon
        # detecta deslocamento bem abaixo de 20%. O que o limiar protege e contra
        # contar como efeito um gene cujo deslocamento o arredondamento para
        # inteiro apagaria.
        v = fc[linha, col_pie]
        t = tested[linha, col_pie]
        seguro = np.where(v > 0.0, v, 1.0)
        movente = (t & np.isfinite(v) & (v > 0.0)
                   & (np.abs(np.log2(seguro)) > limiar_movente))
        candidatos.append((no_painel, int(movente.sum()), nome, linha))

    if not candidatos:
        raise SystemExit(
            f"nenhum nocaute de {fonte_efeitos} em {nome_ctx} atende o filtro — confira "
            f"pert_counts.csv e {diretorio}/preprocessed/meta.json"
        )
    # Alvo do painel primeiro, depois o efeito mais forte: se o proxy couber nos
    # 12 do painel, usa so eles, e o relatorio diz quantos sao.
    candidatos.sort(key=lambda c: (not c[0], -c[1]))
    escolhidos = [c for c in candidatos if c[1] >= min_genes_moventes][:n_alvos]
    if len(escolhidos) < n_alvos:
        raise SystemExit(
            f"so {len(escolhidos)} nocautes tem >= {min_genes_moventes} genes moventes a "
            f"|log2fc| > {limiar_movente}; pedidos {n_alvos}. Baixe --min-genes-moventes, "
            "baixe --n-perts, ou tire --so-alvos-vcc (que limita a 12)."
        )

    out = np.ones((len(escolhidos), eixo_genes.size), dtype=np.float64)
    efetivos = 0
    for k, (_painel, _n, _nome, linha) in enumerate(escolhidos):
        v = fc[linha, col_pie].astype(np.float64)
        t = tested[linha, col_pie]
        valido = t & np.isfinite(v) & (v > 0.0)
        out[k, col_nossa[valido]] = v[valido]
        efetivos += int(valido.sum())
    return Efeitos(
        alvos=[nome for _p, _n, nome, _l in escolhidos],
        no_painel_vcc=[bool(c[0]) for c in escolhidos],
        fold_changes=out,
        n_genes_cobertos=int(col_nossa.size),
        n_pares_efetivos=efetivos,
        fonte=fonte_efeitos,
        contexto_fonte=nome_ctx,
    )


def _emitir_bloco(ctrl: sp.csr_matrix, indices_fonte: np.ndarray,
                  fold_change: np.ndarray) -> sp.csr_matrix:
    """400 celulas perturbadas: controle reamostrado, escalado por gene, arredondado.

    O mecanismo e o do `baseline.py` da Arc — sorteia controles e escala por gene —
    e o contrato de submissao exige contagem INTEIRA, nao negativa e finita. O
    arredondamento e nosso: a emissao `dispersed` da Arc devolve float32.

    Multiplica `data` por `fold_change[indices]` no lugar, sem densificar: um bloco
    denso de 400 x 18 533 em float64 seria 59 MB por perturbacao, e o ganho em
    velocidade nao paga o pico de memoria quando o bloco e 3% denso.
    """
    bloco = ctrl[indices_fonte].tocsr()
    bloco.data *= fold_change[bloco.indices]
    np.rint(bloco.data, out=bloco.data)
    bloco.eliminate_zeros()  # zero explicito e recusado pelo contrato de submissao
    return bloco.astype(np.float32)


def _montar_lado(ctrl: sp.csr_matrix, eixo_genes: np.ndarray, alvos: list[str],
                 fold_changes: np.ndarray, particoes: list[np.ndarray],
                 celulas_por_pert: int, rng: np.random.Generator,
                 *, controles: np.ndarray | None) -> ad.AnnData:
    """Um lado (real ou predito) como AnnData de contagens cruas.

    `controles` nao-nulo acrescenta o pool de controle com rotulo `non-targeting`.

    ⚠️ OS DOIS LADOS PRECISAM DO POOL, e isso contraria o que o contrato de
    submissao do VCC diz. Medido em 2026-10-22: `cell_eval2.io.validate_pair`
    exige conjunto de perturbacoes IDENTICO nos dois lados e aborta com
    "perturbation sets differ: 4 pred vs 5 real; real-only=['non-targeting']".
    O arquivo `.vcc` que se envia ao desafio nao pode ter linha `non-targeting`,
    mas o h5ad que chega ao `run` tem — o passo de preparacao da Arc reinjeta os
    controles reais na predicao, e e esse artefato que o avaliador le.

    O conteudo dessas linhas e inerte no caminho de CPU: com
    `control_source: real` (o que o preset `vcc2026` fixa), `run._pred_de_input`
    DESCARTA as linhas de controle da predicao e dobra o pool de controle do lado
    REAL como grupo de referencia. Elas existem para satisfazer o validador.
    """
    blocos: list[sp.csr_matrix] = []
    rotulos: list[np.ndarray] = []
    if controles is not None:
        blocos.append(ctrl[controles].tocsr().astype(np.float32))
        rotulos.append(np.full(controles.size, ROTULO_CONTROLE, dtype=object))
    for k, alvo in enumerate(alvos):
        fonte = rng.choice(particoes[k], size=celulas_por_pert, replace=True)
        blocos.append(_emitir_bloco(ctrl, fonte, fold_changes[k]))
        rotulos.append(np.full(celulas_por_pert, alvo, dtype=object))

    X = sp.vstack(blocos, format="csr")
    del blocos
    pert = np.concatenate(rotulos)
    obs = pd.DataFrame(
        {COLUNA_PERT: pd.Categorical(pert)},
        index=pd.Index([f"c{i:07d}" for i in range(pert.size)], name=None),
    )
    var = pd.DataFrame(index=pd.Index(eixo_genes, name=None))
    return ad.AnnData(X=X, obs=obs, var=var)


def _particionar(n_ctrl_total: int, pool_controle: np.ndarray, n_alvos: int,
                 rng: np.random.Generator, *, fonte: str) -> list[np.ndarray]:
    """Pools de celulas-fonte, um por alvo, sempre disjuntos ENTRE SI.

    Disjunto entre si importa por dois motivos medidos no codigo do cell-eval2: a
    ancora de replicata divide o lado real em metades e pontua metade contra
    metade, e o `pds_cosine` ordena cada perfil previsto contra o de TODAS as
    perturbacoes. Duas perturbacoes que reamostram as mesmas celulas-fonte
    compartilham ruido de amostragem e inflam as duas coisas.

    ⚠️ **VAZAMENTO MEDIDO, e e o oposto do que o paragrafo acima protege.** As
    particoes sao disjuntas entre PERTURBACOES, mas o lado real e o lado predito
    de uma MESMA perturbacao sorteiam da mesma particao. Com 40 perturbacoes cada
    particao tem ~460 celulas, e duas amostras com reposicao do mesmo pool de 460
    se parecem mais entre si do que amostras de outro pool — rastro de identidade
    de celula que o `pds_cosine` le como acerto. Medido em 2026-10-22, com 40
    perturbacoes, para uma predicao de CONTROLE PURO (delta zero, que nao deveria
    discriminar nada):

        particao por perturbacao   pds bruto 0,99936   avg_score +0,188
        pool unico (`--pred-pool-unico`) pds bruto 0,41859   avg_score -0,033

    Ou seja o rastro valia +1,106 de `pds_cosine` bruto e +0,221 de `avg_score`, e
    era ele — nao biologia — que punha o controle puro ACIMA de zero. Com o rastro
    removido o controle cai abaixo de zero, que e onde o desenho do modelo assume
    que ele esta.

    Consequencia pratica: **para qualquer braco de referencia emitido por este
    script, use `--pred-pool-unico`.** Uma submissao real nao precisa da flag — o
    `emissor_400` sorteia do pool de controle inteiro, logo ja esta no regime sem
    rastro. A flag existe para que as referencias (oraculo, controle) sejam
    comparaveis com ela.

    `fonte` decide a relacao com o POOL DE CONTROLE, e e escolha, nao detalhe:

    * `disjunta` — as fontes saem das celulas que o pool de controle nao usou.
      Nenhuma celula aparece nos dois lados, nem escalada. Exige
      `n_ctrl < 18400`, porque o painel so entrega 18 400 celulas de controle e
      nao entrega as 138 400 de verdade-fundamental.
    * `compartilhada` — as fontes particionam as 18 400 inteiras, de modo que o
      pool de controle pode ser o painel completo. Cada celula perturbada e
      entao uma copia escalada de uma celula que tambem esta na referencia. E o
      unico jeito de ter os 18 400 controles oficiais com dado nosso, e o preco e
      pequeno: 400 celulas de 18 400 sao 2,2% da media da referencia.
    """
    if fonte == "disjunta":
        livres = np.setdiff1d(np.arange(n_ctrl_total), pool_controle, assume_unique=False)
    else:
        livres = np.arange(n_ctrl_total)
    rng.shuffle(livres)
    partes = np.array_split(livres, n_alvos)
    magras = [i for i, p in enumerate(partes) if p.size < 50]
    if magras:
        raise SystemExit(
            f"particao {magras[0]} tem {partes[magras[0]].size} celulas-fonte (< 50). "
            f"{livres.size} celulas livres nao sustentam {n_alvos} pools disjuntos com "
            f"--fonte {fonte} — baixe --n-ctrl, baixe --n-perts, ou use "
            "--fonte compartilhada."
        )
    return [np.asarray(p) for p in partes]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__ or "", formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--contexto", default="D", choices=["D", "E", "F"])
    p.add_argument("--fonte-efeitos", default="arc_vcc_25", choices=sorted(FONTES),
                   help="de onde vem a razao de dobra medida. `orion` cobre 148 dos 300 "
                        "alvos do painel contra 12 de `arc_vcc_25` [default: %(default)s]")
    p.add_argument("--contexto-pie", default=None,
                   help="contexto dentro da fonte de efeitos; o default de cada fonte esta "
                        "em FONTES")
    p.add_argument("--n-perts", type=int, default=8,
                   help="perturbacoes no par (real, predito) [default: %(default)s]")
    p.add_argument("--n-ctrl", type=int, default=4000,
                   help="celulas no pool de controle do lado real; 18400 e o painel "
                        "inteiro [default: %(default)s]")
    p.add_argument("--fonte", default="compartilhada", choices=["disjunta", "compartilhada"],
                   help="de onde saem as celulas-fonte perturbadas em relacao ao pool de "
                        "controle; ver _particionar [default: %(default)s]")
    p.add_argument("--celulas-por-pert", type=int, default=400,
                   help="fixado em 400 pelo painel oficial [default: %(default)s]")
    p.add_argument("--encolhimento", type=float, default=0.5,
                   help="expoente aplicado a razao de dobra no lado predito: 1.0 copia o "
                        "efeito real, 0.0 prediz ausencia de efeito [default: %(default)s]")
    p.add_argument("--jitter", type=float, default=0.08,
                   help="desvio lognormal por gene somado ao efeito predito, para que a "
                        "predicao nao seja o real reescalado [default: %(default)s]")
    p.add_argument("--pred-pool-unico", action="store_true",
                   help="o lado PREDITO sorteia toda perturbacao do pool de controle inteiro, "
                        "em vez da particao por perturbacao. Serve para medir quanto do "
                        "`pds_cosine` vem de identidade de celula compartilhada entre os dois "
                        "lados e nao de acerto biologico. Nao altera o lado real, logo o "
                        "bundle continua valendo")
    p.add_argument("--min-genes-moventes", type=int, default=200,
                   help="piso de genes com efeito por alvo [default: %(default)s]")
    p.add_argument("--limiar-movente", type=float, default=0.1,
                   help="|log2fc| acima do qual um gene conta como movente [default: %(default)s]")
    p.add_argument("--so-alvos-vcc", action="store_true",
                   help="usa somente nocautes que sao alvos do painel VCC 2026. Teto "
                        "medido: 12 com arc_vcc_25/ARC_H1, 148 com orion/hct116. "
                        "Obrigatorio para recortar uma submissao real contra este lado")
    p.add_argument("--semente", type=int, default=0)
    p.add_argument("--saida", required=True, help="diretorio onde gravar real.h5ad e pred.h5ad")
    args = p.parse_args()

    saida = Path(args.saida)
    saida.mkdir(parents=True, exist_ok=True)

    eixo_genes = pd.read_csv(CONTROLES / "gene_names.csv")["gene_name"].to_numpy(dtype=str)
    alvos_vcc = set(pd.read_csv(CONTROLES / "pert_counts.csv")["target_gene"].astype(str))

    efeitos = _carregar_efeitos_pie(eixo_genes, alvos_vcc, args.n_perts,
                                    min_genes_moventes=args.min_genes_moventes,
                                    limiar_movente=args.limiar_movente,
                                    so_alvos_vcc=args.so_alvos_vcc,
                                    fonte_efeitos=args.fonte_efeitos,
                                    contexto_pie=args.contexto_pie)

    caminho = CONTROLES / f"context_{args.contexto}.h5ad"
    print(f"lendo {caminho.name} ...", file=sys.stderr)
    ctrl_ad = ad.read_h5ad(caminho)
    if list(ctrl_ad.var_names) != list(eixo_genes):
        raise SystemExit("o eixo de genes do h5ad nao bate com gene_names.csv")
    ctrl = sp.csr_matrix(ctrl_ad.X)
    del ctrl_ad

    # Um gerador por papel, todos derivados da mesma semente: `np.random.seed` e
    # estado global de processo e qualquer biblioteca importada depois consome dele.
    semeador = np.random.SeedSequence(args.semente)
    rng_ctrl, rng_real, rng_pred, rng_jitter = (
        np.random.default_rng(s) for s in semeador.spawn(4)
    )

    n_total = _forma(ctrl)[0]
    pool_controle = rng_ctrl.choice(n_total, size=args.n_ctrl, replace=False)
    particoes = _particionar(n_total, pool_controle, args.n_perts, rng_ctrl, fonte=args.fonte)

    # O lado predito: efeito encolhido mais jitter lognormal por (alvo, gene). O
    # encolhimento sozinho seria o real reescalado, e `delta_pearson`/`pds_cosine`
    # nao distinguiriam isso de um modelo perfeito com escala errada.
    fc_pred = efeitos.fold_changes ** args.encolhimento
    if args.jitter > 0.0:
        fc_pred = fc_pred * np.exp(rng_jitter.normal(0.0, args.jitter, size=fc_pred.shape))

    print(f"montando lado real: {args.n_ctrl} controles + "
          f"{args.n_perts}x{args.celulas_por_pert} perturbadas", file=sys.stderr)
    real = _montar_lado(ctrl, eixo_genes, efeitos.alvos, efeitos.fold_changes, particoes,
                        args.celulas_por_pert, rng_real, controles=pool_controle)
    real.write_h5ad(saida / "real.h5ad", compression=None)
    forma_real, nnz_real = real.shape, int(real.X.nnz)  # type: ignore[union-attr]
    del real

    print("montando lado predito", file=sys.stderr)
    # ⚠️ O lado real e o predito sorteiam, por padrao, da MESMA particao por
    # perturbacao. Isso deixa um rastro de identidade de celula que o
    # `pds_cosine` consegue ler: ele pergunta se o perfil previsto de p esta mais
    # perto do perfil real de p do que do de qualquer outra perturbacao, e duas
    # amostras com reposicao do mesmo pool de ~460 celulas se parecem mais entre
    # si do que com amostras de outro pool, independente de biologia.
    # `--pred-pool-unico` tira o rastro, para medir quanto dele havia.
    particoes_pred = ([np.arange(n_total)] * args.n_perts if args.pred_pool_unico
                      else particoes)
    # MESMAS celulas de controle do lado real, nao uma reamostragem: e o que o
    # passo de preparacao da Arc faz ao reinjetar os controles na submissao, e
    # mantem `pds_cosine` e as metricas de expressao lendo um unico pool.
    pred = _montar_lado(ctrl, eixo_genes, efeitos.alvos, fc_pred, particoes_pred,
                        args.celulas_por_pert, rng_pred, controles=pool_controle)
    pred.write_h5ad(saida / "pred.h5ad", compression=None)
    forma_pred, nnz_pred = pred.shape, int(pred.X.nnz)  # type: ignore[union-attr]
    del pred

    resumo = {
        "contexto": args.contexto,
        "fonte_celulas": str(caminho.relative_to(RAIZ)).replace("\\", "/"),
        "fonte_efeitos": efeitos.fonte,
        "contexto_fonte_efeitos": efeitos.contexto_fonte,
        "alvos": efeitos.alvos,
        "n_alvos_no_painel_vcc": sum(efeitos.no_painel_vcc),
        "so_alvos_vcc": args.so_alvos_vcc,
        "limiar_movente": args.limiar_movente,
        "min_genes_moventes": args.min_genes_moventes,
        "n_perts": args.n_perts,
        "n_ctrl": args.n_ctrl,
        "fonte": args.fonte,
        "pred_pool_unico": args.pred_pool_unico,
        "celulas_por_pert": args.celulas_por_pert,
        "encolhimento_pred": args.encolhimento,
        "jitter_pred": args.jitter,
        "semente": args.semente,
        "genes_do_eixo_cobertos_pelo_pie": efeitos.n_genes_cobertos,
        "pares_alvo_gene_com_fold_change": efeitos.n_pares_efetivos,
        "real": {"forma": list(forma_real), "nnz": nnz_real,
                 "nnz_por_celula": round(nnz_real / forma_real[0], 1)},
        "pred": {"forma": list(forma_pred), "nnz": nnz_pred,
                 "nnz_por_celula": round(nnz_pred / forma_pred[0], 1)},
    }
    (saida / "proxy.json").write_text(json.dumps(resumo, indent=2, ensure_ascii=False) + "\n",
                                      encoding="utf-8")
    print(json.dumps(resumo, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
