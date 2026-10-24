#!/usr/bin/env python
"""Validador de CONTRATO da submissao do VCC 2026: 28 itens, zero numerico.

Confere estrutura, rotulo e eixo de um `.h5ad` antes do `vcc prep`. NAO confere
se os numeros estao certos -- um arquivo de delta zero passa 28/28 aqui. Para
isso existe `validar_numerico.py`, e os DOIS rodam antes de qualquer envio.

O contrato vem do `manifest.json` do pacote `controls-final` (coluna de alvo,
coluna de contexto, rotulo de controle, 18.533 genes, 300 alvos, 400 celulas
por alvo, contextos D/E/F) e das travas do `vcc prep` 0.3.0 (teto de celulas,
de contagem por celula e de nao nulos).

Leitura em blocos de linhas direto do HDF5: a submissao completa tem ~2,1e9
nao nulos e nunca e materializada. `obs`, `var` e `X` sao lidos crus, sem o
carregador do anndata, porque o arquivo pode estar malformado de proposito.

Se o `indptr` estiver quebrado (itens 13-18), o validador PARA no item 18: o
mapeamento linha -> valor fica indefinido e os itens seguintes mediriam lixo.
Se `X` nao for CSR, para no item 12, pelo mesmo motivo.

Uso:
    uv run python scripts/vcc/validar_submissao.py submissao.h5ad
    uv run python scripts/vcc/validar_submissao.py forja.h5ad --contextos D --json laudo.json
Sai com 0 se aprovado e 1 se reprovado.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[2]
REF_PADRAO = RAIZ / "data" / "raw" / "vcc" / "controls-final"

# Travas do `vcc prep` 0.3.0 (`vcc prep --help`).
TETO_DE_CELULAS = 400_000  # --max-cell-dim
TETO_DE_CONTAGEM = 1_000_000  # --max-counts-per-cell
TETO_DE_NNZ = 4_750_000_000  # --max-nnz, "the scoring hardware's limit"

LINHAS_POR_BLOCO = 4096
MAX_EXEMPLOS = 10


@dataclass(frozen=True)
class Item:
    nome: str
    passou: bool
    detalhe: str


@dataclass
class Laudo:
    arquivo: str
    contextos_exigidos: list[str]
    painel_completo: bool
    itens: list[Item] = field(default_factory=list)

    def anotar(self, nome: str, passou: bool, detalhe: str) -> bool:
        self.itens.append(Item(nome, bool(passou), detalhe))
        return bool(passou)

    @property
    def falhas(self) -> list[str]:
        return [i.nome for i in self.itens if not i.passou]

    @property
    def itens_testados(self) -> int:
        return len(self.itens)

    @property
    def veredito(self) -> str:
        return "reprovado" if self.falhas else "aprovado"

    def para_dict(self) -> dict:
        return {
            "arquivo": self.arquivo,
            "contextos_exigidos": self.contextos_exigidos,
            "painel_completo": self.painel_completo,
            "itens": [asdict(i) for i in self.itens],
            "itens_testados": self.itens_testados,
            "falhas": self.falhas,
            "veredito": self.veredito,
        }


@dataclass(frozen=True)
class Referencia:
    genes: list[str]
    alvos: list[str]
    celulas_esperadas: dict[str, int]
    contextos: list[str]
    col_alvo: str
    col_contexto: str
    controle: str


def ler_referencia(pasta: Path) -> Referencia:
    manifest = json.loads((pasta / "manifest.json").read_text(encoding="utf-8"))
    genes = pd.read_csv(pasta / "gene_names.csv")["gene_name"].astype(str).tolist()
    perts = pd.read_csv(pasta / "pert_counts.csv")
    alvos = perts["target_gene"].astype(str).tolist()
    # Uma coluna `n_cells` no pert_counts vence o valor do manifest, como no `vcc prep`.
    if "n_cells" in perts.columns:
        esperadas = dict(zip(alvos, perts["n_cells"].astype(int), strict=True))
    else:
        esperadas = dict.fromkeys(alvos, int(manifest["cells_per_pert"]))
    return Referencia(
        genes=genes,
        alvos=alvos,
        celulas_esperadas=esperadas,
        contextos=[str(c) for c in manifest["contexts"]],
        col_alvo=str(manifest["pert_col"]),
        col_contexto=str(manifest["context_col"]),
        controle=str(manifest["control_label"]),
    )


def _lista(valores: Sequence[str]) -> str:
    v = list(valores)
    return repr(v) if len(v) <= MAX_EXEMPLOS else repr(v[:MAX_EXEMPLOS])[:-1] + ", ...]"


def _numero(x: float) -> str:
    return str(int(x)) if np.isfinite(x) and float(x).is_integer() else f"{x:g}"


def _ler_tabela(f: h5py.File, chave: str) -> pd.DataFrame:
    g = f.get(chave)
    if not isinstance(g, h5py.Group):
        return pd.DataFrame()
    tabela = ad.io.read_elem(g)
    return tabela if isinstance(tabela, pd.DataFrame) else pd.DataFrame()


def _encoding(no: h5py.Group | h5py.Dataset) -> str:
    v = no.attrs.get("encoding-type", "")
    return v.decode() if isinstance(v, bytes) else str(v)


def _itens_de_rotulo(laudo: Laudo, ref: Referencia, obs: pd.DataFrame, exigidos: list[str]):
    """Itens 9-12. Devolve (alvo, contexto) como arrays de str, ou None se a coluna falta."""
    alvo = obs[ref.col_alvo].astype(str).to_numpy() if ref.col_alvo in obs else None
    ctx = obs[ref.col_contexto].astype(str).to_numpy() if ref.col_contexto in obs else None
    ausente = f"nao avaliado: coluna {ref.col_contexto!r} ausente em obs"
    if alvo is None:
        ausente = f"nao avaliado: coluna {ref.col_alvo!r} ausente em obs"
        laudo.anotar("sem_linha_de_controle", False, ausente)
    else:
        n_ctrl = int((alvo == ref.controle).sum())
        laudo.anotar("sem_linha_de_controle", n_ctrl == 0, f"{n_ctrl} celula(s) rotuladas {ref.controle!r}")

    if alvo is None or ctx is None:
        for nome in ("contextos", "alvos_por_contexto", "celulas_por_alvo"):
            laudo.anotar(nome, False, ausente)
        return None

    presentes = sorted(set(ctx.tolist()))
    laudo.anotar(
        "contextos",
        set(presentes) == set(exigidos),
        f"presentes {presentes}, esperados {exigidos}",
    )

    oficiais = set(ref.alvos)
    problemas = []
    for c in exigidos:
        achados = set(alvo[(ctx == c) & (alvo != ref.controle)].tolist())
        faltam, sobram = oficiais - achados, achados - oficiais
        if faltam or sobram:
            problemas.append(
                f"{c}: faltam {len(faltam)} {_lista(sorted(faltam))}, "
                f"sobram {len(sobram)} {_lista(sorted(sobram))}"
            )
    laudo.anotar(
        "alvos_por_contexto",
        not problemas,
        "; ".join(problemas) or f"{len(oficiais)} alvos oficiais em cada contexto",
    )

    contagem = Counter(zip(ctx.tolist(), alvo.tolist(), strict=True))
    errados = []
    for c in exigidos:
        fora = [
            f"{a}={contagem.get((c, a), 0)}"
            for a in ref.alvos
            if contagem.get((c, a), 0) != ref.celulas_esperadas[a]
        ]
        if fora:
            errados.append(f"{c}: {len(fora)} alvo(s) {_lista(fora)}")
    unicos = sorted(set(ref.celulas_esperadas.values()))
    esperado = f"exatamente {unicos[0]}" if len(unicos) == 1 else "o n_cells oficial"
    laudo.anotar(
        "celulas_por_alvo", not errados, "; ".join(errados) or f"{esperado} por (contexto, alvo)"
    )
    return alvo, ctx


def _itens_de_estrutura(laudo: Laudo, indptr: np.ndarray, n_linhas: int, nd: int, ni: int) -> bool:
    """Itens 13-18, so sobre o indptr e os comprimentos. Devolve se o indptr e confiavel."""
    ok = laudo.anotar(
        "estrutura_indptr_tamanho",
        indptr.size == n_linhas + 1,
        f"{indptr.size} entradas, esperado {n_linhas + 1}",
    )
    inicio = int(indptr[0]) if indptr.size else None
    ok &= laudo.anotar(
        "estrutura_indptr_inicio",
        inicio == 0,
        f"indptr[0] = {inicio}" if inicio is not None else "indptr vazio",
    )
    negativas = int((np.diff(indptr) < 0).sum()) if indptr.size > 1 else 0
    ok &= laudo.anotar(
        "estrutura_indptr_monotonico", negativas == 0, f"{negativas} linha(s) com comprimento negativo"
    )
    ok &= laudo.anotar("estrutura_pares", nd == ni, f"data {nd:,} vs indices {ni:,}")
    fim = int(indptr[-1]) if indptr.size else 0
    ok &= laudo.anotar(
        "estrutura_fim_dentro", 0 <= fim <= min(nd, ni), f"indptr[-1] = {fim:,}, data = {nd:,}"
    )
    folga = max(0, max(nd, ni) - fim)
    ok &= laudo.anotar(
        "canonico_sem_folga", folga == 0, f"{folga} entrada(s) de folga depois de indptr[-1]"
    )
    return ok


@dataclass
class _Acumulador:
    """O que se soma entre blocos de linhas para os itens 19-27."""

    fora: int = 0
    desordem: int = 0
    zeros: int = 0
    nao_finitos: int = 0
    negativos: int = 0
    fracionarios: int = 0
    maximo: float = -np.inf
    acima_do_teto: int = 0
    maior_soma: float = 0.0
    celulas_nulas: int = 0
    primeira: dict[str, int] = field(default_factory=dict)

    def marcar_primeira(self, chave: str, mascara: np.ndarray, rel: np.ndarray, r0: int) -> None:
        if chave in self.primeira or not mascara.any():
            return
        p = int(np.flatnonzero(mascara)[0])
        # Com linhas vazias, varias linhas tem o mesmo inicio; side="right" - 1
        # cai na ultima delas, que e a nao vazia que contem a posicao p.
        self.primeira[chave] = r0 + int(np.searchsorted(rel, p, side="right")) - 1

    def sufixo(self, chave: str) -> str:
        linha = self.primeira.get(chave)
        return f" primeira celula afetada: linha {linha}" if linha is not None else ""


def _varrer_valores(
    x: h5py.Group,
    indptr: np.ndarray,
    ng: int,
    grupos: np.ndarray | None,
    n_grupos: int,
    linhas_por_bloco: int,
) -> tuple[_Acumulador, np.ndarray | None]:
    dados, indices = x["data"], x["indices"]
    assert isinstance(dados, h5py.Dataset) and isinstance(indices, h5py.Dataset)
    acc = _Acumulador()
    nao_nulas_por_grupo = np.zeros(n_grupos, dtype=np.int64) if grupos is not None else None
    n_linhas = indptr.size - 1
    for r0 in range(0, n_linhas, linhas_por_bloco):
        r1 = min(r0 + linhas_por_bloco, n_linhas)
        a, b = int(indptr[r0]), int(indptr[r1])
        dat = np.asarray(dados[a:b])
        ind = np.asarray(indices[a:b]).astype(np.int64)
        rel = indptr[r0 : r1 + 1] - a

        acc.fora += int(((ind < 0) | (ind >= ng)).sum())
        if ind.size > 1:
            # O par (ultimo de uma linha, primeiro da seguinte) nao e uma violacao.
            d = np.diff(ind)
            dentro = np.ones(d.size, dtype=bool)
            fronteiras = rel[1:-1] - 1
            dentro[fronteiras[(fronteiras >= 0) & (fronteiras < d.size)]] = False
            acc.desordem += int(((d <= 0) & dentro).sum())

        acc.zeros += int((dat == 0).sum())
        finito = np.isfinite(dat)
        negativo = finito & (dat < 0)
        fracionario = finito & (dat != np.floor(dat))
        acc.nao_finitos += int((~finito).sum())
        acc.negativos += int(negativo.sum())
        acc.fracionarios += int(fracionario.sum())
        acc.marcar_primeira("nao_finitos", ~finito, rel, r0)
        acc.marcar_primeira("negativos", negativo, rel, r0)
        acc.marcar_primeira("fracionarios", fracionario, rel, r0)
        if finito.any():
            acc.maximo = max(acc.maximo, float(dat[finito].max()))

        # Soma por celula por soma acumulada: reduceat erra em linha vazia, e o
        # NaN/inf e zerado antes para nao contaminar as linhas seguintes do bloco.
        acumulada = np.concatenate([[0.0], np.cumsum(np.where(finito, dat, 0), dtype=np.float64)])
        somas = acumulada[rel[1:]] - acumulada[rel[:-1]]
        acc.acima_do_teto += int((somas > TETO_DE_CONTAGEM).sum())
        if somas.size:
            acc.maior_soma = max(acc.maior_soma, float(somas.max()))
        acc.celulas_nulas += int((somas == 0).sum())
        if grupos is not None and nao_nulas_por_grupo is not None:
            nao_nulas_por_grupo += np.bincount(
                grupos[r0:r1], weights=(somas != 0), minlength=n_grupos
            ).astype(np.int64)
    return acc, nao_nulas_por_grupo


def validar(
    caminho: Path,
    referencia: Path = REF_PADRAO,
    contextos: Sequence[str] | None = None,
    linhas_por_bloco: int = LINHAS_POR_BLOCO,
) -> Laudo:
    ref = ler_referencia(Path(referencia))
    exigidos = list(contextos) if contextos is not None else ref.contextos
    laudo = Laudo(str(caminho), exigidos, sorted(exigidos) == sorted(ref.contextos))

    with h5py.File(caminho, "r") as f:
        obs, var = _ler_tabela(f, "obs"), _ler_tabela(f, "var")
        x = f.get("X")

        csr = isinstance(x, h5py.Group) and _encoding(x) == "csr_matrix"
        if csr:
            detalhe_x = "encoding-type = 'csr_matrix'"
        elif isinstance(x, h5py.Dataset):
            detalhe_x = f"X e denso, forma {tuple(x.shape)}; exigido csr_matrix"
        elif isinstance(x, h5py.Group):
            detalhe_x = f"encoding-type = {_encoding(x)!r}; exigido 'csr_matrix'"
        else:
            detalhe_x = "X ausente"
        laudo.anotar("matriz_csr", csr, detalhe_x)

        colunas = [str(c) for c in obs.columns]
        laudo.anotar(
            "colunas_obs",
            ref.col_alvo in colunas and ref.col_contexto in colunas,
            f"colunas {colunas}; exigidas {ref.col_alvo!r} e {ref.col_contexto!r}",
        )

        n_obs, n_var = len(obs), len(var)
        if isinstance(x, h5py.Group) and "shape" in x.attrs:
            forma = tuple(int(v) for v in np.asarray(x.attrs["shape"]).tolist())
        elif isinstance(x, h5py.Dataset):
            forma = tuple(int(v) for v in x.shape)
        else:
            forma = (-1, -1)
        coerente = laudo.anotar(
            "forma_coerente",
            forma == (n_obs, n_var),
            f"X {forma[0]} x {forma[1]}, obs {n_obs}, var {n_var}",
        )

        ng = len(ref.genes)
        laudo.anotar("dimensao_de_genes", n_var == ng, f"{n_var} colunas, exigidas {ng}")
        eixo = var.index.astype(str).tolist()
        if eixo == ref.genes:
            detalhe_eixo = "identico a gene_names.csv na ordem"
        else:
            # strict=False de proposito: comprimentos diferentes sao um dos casos tratados.
            pares_eixo = zip(eixo, ref.genes, strict=False)
            pos = next((i for i, (u, v) in enumerate(pares_eixo) if u != v), None)
            detalhe_eixo = (
                f"divergente (primeira diferenca em posicao {pos}: {eixo[pos]!r} vs {ref.genes[pos]!r})"
                if pos is not None
                else f"divergente (comprimento {len(eixo)} vs {ng})"
            )
        laudo.anotar("eixo_de_genes", eixo == ref.genes, detalhe_eixo)
        repetidos = len(eixo) - len(set(eixo))
        laudo.anotar("genes_sem_repetido", repetidos == 0, f"{repetidos} simbolo(s) repetido(s)")

        esperado = len(exigidos) * sum(ref.celulas_esperadas.values())
        laudo.anotar(
            "total_de_celulas", n_obs == esperado, f"{n_obs:,} celulas, esperado {esperado:,}"
        )
        laudo.anotar(
            "teto_de_celulas",
            n_obs <= TETO_DE_CELULAS,
            f"{n_obs:,} celulas, teto {TETO_DE_CELULAS:,}",
        )

        rotulos = _itens_de_rotulo(laudo, ref, obs, exigidos)

        if not (csr and coerente) or not isinstance(x, h5py.Group):
            return laudo

        indptr = np.asarray(x["indptr"]).astype(np.int64)
        dados, indices = x["data"], x["indices"]
        assert isinstance(dados, h5py.Dataset) and isinstance(indices, h5py.Dataset)
        if not _itens_de_estrutura(laudo, indptr, n_obs, dados.shape[0], indices.shape[0]):
            return laudo

        grupos, rotulos_de_grupo = None, np.array([], dtype=object)
        if rotulos is not None:
            alvo, ctx = rotulos
            grupos, rotulos_de_grupo = pd.factorize(pd.Series(ctx) + "/" + pd.Series(alvo))
        acc, nao_nulas = _varrer_valores(
            x, indptr, ng, grupos, len(rotulos_de_grupo), linhas_por_bloco
        )

    laudo.anotar("estrutura_indices_em_faixa", acc.fora == 0, f"{acc.fora} indice(s) fora de 0..{ng - 1}")
    laudo.anotar(
        "canonico_ordem_estrita",
        acc.desordem == 0,
        f"{acc.desordem} par(es) fora de ordem ou com indice repetido",
    )
    laudo.anotar("sem_zeros_explicitos", acc.zeros == 0, f"{acc.zeros} zero(s) armazenado(s)")
    laudo.anotar(
        "valores_finitos", acc.nao_finitos == 0, f"{acc.nao_finitos} NaN/inf.{acc.sufixo('nao_finitos')}"
    )
    laudo.anotar(
        "valores_nao_negativos",
        acc.negativos == 0,
        f"{acc.negativos} valor(es) negativo(s).{acc.sufixo('negativos')}",
    )
    laudo.anotar(
        "valores_inteiros",
        acc.fracionarios == 0,
        f"{acc.fracionarios} valor(es) fracionario(s), maximo {_numero(acc.maximo)}."
        f"{acc.sufixo('fracionarios')}",
    )
    laudo.anotar(
        "contagem_por_celula",
        acc.acima_do_teto == 0,
        f"{acc.acima_do_teto} celula(s) acima de {TETO_DE_CONTAGEM:,}; maior soma {_numero(acc.maior_soma)}",
    )
    laudo.anotar(
        "celula_nao_nula",
        acc.celulas_nulas == 0,
        f"{acc.celulas_nulas} celula(s) com soma zero (a normalizacao por CPM divide por isso)",
    )
    if nao_nulas is None:
        laudo.anotar("perturbacao_nao_nula", False, "nao avaliado: colunas de rotulo ausentes em obs")
    else:
        nulas = [str(rotulos_de_grupo[i]) for i in np.flatnonzero(nao_nulas == 0)]
        laudo.anotar(
            "perturbacao_nao_nula",
            not nulas,
            f"{len(nulas)} perturbacao(oes) inteiramente nula(s): {_lista(nulas)}"
            if nulas
            else "nenhuma perturbacao inteiramente nula",
        )

    nnz = int(indptr[-1])
    por_celula = np.diff(indptr)
    laudo.anotar(
        "teto_de_densidade",
        nnz <= TETO_DE_NNZ,
        f"{nnz:,} entradas armazenadas, teto {TETO_DE_NNZ:,} ({100 * nnz / TETO_DE_NNZ:.2f}%); "
        f"{nnz / max(n_obs, 1):.0f}/celula (min {int(por_celula.min())}, max {int(por_celula.max())})",
    )
    return laudo


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("arquivo", type=Path, help="o .h5ad a validar")
    ap.add_argument("--ref", type=Path, default=REF_PADRAO, help="pasta do controls-final extraido")
    ap.add_argument(
        "--contextos",
        default=None,
        help="contextos exigidos, separados por virgula (padrao: os do manifest, D,E,F)",
    )
    ap.add_argument("--json", type=Path, default=None, help="grava o laudo completo neste arquivo")
    ap.add_argument("--linhas-por-bloco", type=int, default=LINHAS_POR_BLOCO)
    args = ap.parse_args(argv)

    contextos = [c for c in args.contextos.split(",") if c] if args.contextos is not None else None
    laudo = validar(args.arquivo, args.ref, contextos, args.linhas_por_bloco)
    for item in laudo.itens:
        print(f"{'ok   ' if item.passou else 'FALHA'}  {item.nome:<28} {item.detalhe}")
    print(
        f"\n{laudo.veredito.upper()}: {laudo.itens_testados - len(laudo.falhas)}/{laudo.itens_testados}"
        f" itens; contextos {laudo.contextos_exigidos}"
        f"{'' if laudo.painel_completo else ' (painel PARCIAL)'}"
    )
    if args.json is not None:
        args.json.write_text(json.dumps(laudo.para_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if laudo.veredito == "aprovado" else 1


if __name__ == "__main__":
    sys.exit(main())
