"""
Figura de arvore filogenetica com suporte nos ramos, a partir de Newick do
IQ-TREE (rotulos "SH-aLRT/aBayes/UFBoot", so os campos pedidos) ou do
wASTRAL (probabilidade local).

    uv run python scripts/filogenia/desenhar_arvore.py arvore.treefile saida --outgroup Homo_sapiens

Regras de figura da stack: eixo rotulado com unidade, n declarado no titulo.
Ramo forte (todos os criterios presentes passam: SH-aLRT >= 80, UFBoot >= 95,
aBayes ou probabilidade local >= 0,95) fica em preto; o resto, em cinza com o
valor escrito.

Suporte e propriedade de RAMO (biparticao), mas o Newick o guarda no no.
Reenraizar move o no de lugar, e a Biopython so leva junto o comprimento —
o rotulo ficaria no ramo vizinho (medido). Por isso o suporte e lido por
biparticao antes de enraizar e reatribuido por biparticao depois.
"""

from __future__ import annotations

import argparse
import math
import re
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from Bio import Phylo  # noqa: E402
from Bio.Phylo.BaseTree import Clade, Tree  # noqa: E402

#: Ordem dos campos no rotulo do IQ-TREE com os padroes do pipeline.
ROTULOS_PADRAO = ("alrt", "abayes", "ufboot")

#: Limiar por campo; o ramo e forte quando todos os campos presentes passam.
LIMIAR = {"alrt": 80.0, "ufboot": 95.0, "abayes": 0.95, "prob": 0.95}

LEGENDA = {
    "genes": "Ramos em preto: SH-aLRT ≥ 80 e UFBoot ≥ 95 (e aBayes ≥ 0,95), critério do IQ-TREE "
    "para árvore de um loco; em cinza, os valores.",
    "supermatriz": "Supermatriz: o suporte satura com muitos sítios — julgue os ramos pelo gCF/sCFL "
    "(06_concordancia). Preto/cinza só repete o critério SH-aLRT ≥ 80 e UFBoot ≥ 95.",
    "especies": "Ramos em preto: probabilidade local ≥ 0,95; em cinza, o valor. Comprimentos do "
    "wASTRAL híbrido não são unidades coalescentes: figura sem escala.",
}

Suporte = dict[str, float]


def _campos(rotulo: str, rotulos: tuple[str, ...]) -> Suporte | None:
    """
    Le um rotulo de suporte. Usa `rotulos` quando o numero de campos bate; senao
    infere pela contagem (2 = SH-aLRT/UFBoot, o padrao do IQ-TREE com -B e
    --alrt; 3 = SH-aLRT/aBayes/UFBoot). Campo vazio ("/96", "/0.97/96") conta
    na posicao e fica de fora do resultado: descartar antes de casar deslocaria
    os nomes e leria o aBayes como SH-aLRT. Devolve None se nao reconhecer — e
    o desenho mostra o rotulo cru, em cinza.
    """
    try:
        valores = [float(p) if p.strip() else None for p in rotulo.split("/")]
    except ValueError:
        return None
    n = len(valores)
    if n == len(rotulos):
        nomes = rotulos
    elif n == 3:
        nomes = ("alrt", "abayes", "ufboot")
    elif n == 2:
        nomes = ("alrt", "ufboot")
    elif n == 1 and valores[0] is not None:
        nomes = ("prob",) if valores[0] <= 1.0 else ("ufboot",)
    else:
        return None
    lido = {k: v for k, v in zip(nomes, valores, strict=True) if v is not None}
    return lido or None


def _forte(sup: Suporte) -> bool:
    return all(v >= LIMIAR[k] for k, v in sup.items() if k in LIMIAR)


def _texto(sup: Suporte) -> str:
    ordem = [k for k in ("alrt", "abayes", "ufboot", "prob") if k in sup]
    return "/".join(f"{sup[k]:.3g}" if sup[k] <= 1.0 else f"{sup[k]:g}" for k in ordem)


def _chave(folhas: frozenset[str], todas: frozenset[str], ref: str) -> frozenset[str]:
    """Biparticao canonica: o lado que nao tem o taxon de referencia."""
    return todas - folhas if ref in folhas else folhas


def _ler_suportes(
    arv: Tree, rotulos: tuple[str, ...]
) -> dict[frozenset[str], tuple[Suporte | None, str]]:
    todas = frozenset(t.name for t in arv.get_terminals())
    ref = min(todas)
    mapa: dict[frozenset[str], tuple[Suporte | None, str]] = {}
    for c in arv.find_clades():
        if c.is_terminal() or c is arv.root:
            continue
        cru = c.name if c.name else (f"{c.confidence:g}" if c.confidence is not None else "")
        if not cru:
            continue
        folhas = frozenset(t.name for t in c.get_terminals())
        mapa[_chave(folhas, todas, ref)] = (_campos(cru, rotulos), cru)
        c.name, c.confidence = None, None
    return mapa


def enraizar(arv: Tree, outgroup: str) -> str:
    """Enraiza no outgroup presente; devolve o que foi usado (ou "" se nenhum)."""
    nomes = {t.name for t in arv.get_terminals()}
    presentes = [o.strip() for o in outgroup.split(",") if o.strip() in nomes]
    if not presentes:
        arv.root_at_midpoint()
        return ""
    # Primeiro num taxon so (sempre valido); depois, com a arvore ja enraizada
    # ali, o ancestral comum dos outgroups passa a ser um clado bem definido.
    arv.root_with_outgroup(presentes[0])
    if len(presentes) > 1:
        arv.root_with_outgroup(arv.common_ancestor(*presentes))
    return ",".join(presentes)


def desenhar_arvore(
    entrada: Path,
    saida_base: Path,
    titulo: str = "",
    outgroup: str = "",
    unidade: str = "substituições por sítio",
    rotulos: tuple[str, ...] = ROTULOS_PADRAO,
    tipo: str = "genes",
) -> list[Path]:
    # O wASTRAL escreve ":nan" em ramo interno de comprimento indefinido. Lido
    # cru, o parser Newick da Biopython joga o NaN no suporte do no ("0.94:nan"
    # vira confidence=nan) e o valor real some — medido. Tirar o ":nan" antes.
    texto = re.sub(r":nan\b", "", entrada.read_text(encoding="utf-8"), flags=re.IGNORECASE)
    # A Biopython nao declara __all__ em Bio.Phylo; o basedpyright trata read/draw
    # como nao exportados, mas sao a API publica documentada.
    arv = Phylo.read(StringIO(texto), "newick")  # pyright: ignore[reportPrivateImportUsage]
    sem_comprimento = any(
        c.branch_length is None
        or (isinstance(c.branch_length, float) and math.isnan(c.branch_length))
        for c in arv.find_clades()
        if c is not arv.root
    )
    if sem_comprimento:
        for c in arv.find_clades():
            c.branch_length = 1.0

    mapa = _ler_suportes(arv, rotulos)
    raiz = enraizar(arv, outgroup)
    arv.ladderize(reverse=True)

    # Reatribui por biparticao. Os dois filhos da raiz definem a MESMA
    # biparticao (o ramo que a raiz partiu): o suporte vai so para um deles.
    todas = frozenset(t.name for t in arv.get_terminals())
    ref = min(todas)
    sup: dict[int, tuple[Suporte | None, str]] = {}
    usadas: set[frozenset[str]] = set()
    for c in arv.find_clades(order="level"):
        if c.is_terminal() or c is arv.root:
            continue
        chave = _chave(frozenset(t.name for t in c.get_terminals()), todas, ref)
        if chave in mapa and chave not in usadas:
            sup[id(c)] = mapa[chave]
            usadas.add(chave)

    n = arv.count_terminals()
    altura = max(3.0, 0.32 * n + 1.2)
    fig, ax = plt.subplots(figsize=(8.0, altura))

    def rotulo_no(c: Clade) -> str:
        return c.name.replace("_", " ") if c.is_terminal() and c.name else ""

    def rotulo_ramo(c: Clade) -> str | None:
        if c.is_terminal() or id(c) not in sup:
            return None
        valores, cru = sup[id(c)]
        if valores is None:
            return cru  # nao reconhecido: mostra cru, nunca esconde como forte
        return None if _forte(valores) else _texto(valores)

    def cor_ramo(c: Clade) -> str:
        if c.is_terminal() or id(c) not in sup:
            return "#000000"  # raiz criada no reenraizamento: sem suporte, nao fraco
        valores, _ = sup[id(c)]
        return "#000000" if valores is not None and _forte(valores) else "#8c8c8c"

    for c in arv.find_clades():
        c.color = cor_ramo(c)

    Phylo.draw(  # pyright: ignore[reportPrivateImportUsage]
        arv,
        axes=ax,
        do_show=False,
        label_func=rotulo_no,  # pyright: ignore[reportArgumentType]
        branch_labels=rotulo_ramo,
    )
    nomes = {t.name.replace("_", " ") for t in arv.get_terminals() if t.name}
    for txt in ax.texts:
        if txt.get_text().strip() in nomes:
            txt.set_fontstyle("italic")
        else:
            txt.set_fontsize(7)
            txt.set_color("0.35")

    ax.set_ylabel("")
    ax.set_yticks([])
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.set_xlabel("sem escala (cladograma)" if sem_comprimento else unidade)
    if sem_comprimento:
        ax.set_xticks([])
    detalhe = f"n = {n} táxons" + (
        f"; raiz: {raiz.replace('_', ' ')}" if raiz else "; raiz no ponto médio (sem outgroup)"
    )
    ax.set_title(f"{titulo}\n{detalhe}" if titulo else detalhe, fontsize=10, loc="left")
    fig.text(0.01, 0.005, LEGENDA.get(tipo, LEGENDA["genes"]), fontsize=7, color="0.35")
    fig.tight_layout(rect=(0, 0.02, 1, 1))

    saida_base.parent.mkdir(parents=True, exist_ok=True)
    feitos = []
    for ext in ("png", "svg"):
        destino = saida_base.with_suffix(f".{ext}")
        fig.savefig(destino, dpi=200 if ext == "png" else None)
        feitos.append(destino)
    plt.close(fig)
    return feitos


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("arvore", type=Path)
    ap.add_argument("saida", type=Path, help="caminho base, sem extensao (gera .png e .svg)")
    ap.add_argument("--outgroup", default="")
    ap.add_argument("--titulo", default="")
    ap.add_argument("--unidade", default="substituições por sítio")
    ap.add_argument("--tipo", default="genes", choices=sorted(LEGENDA))
    ap.add_argument(
        "--rotulos",
        default=",".join(ROTULOS_PADRAO),
        help="campos do rotulo interno, na ordem do IQ-TREE (resumo.json -> rotulos_suporte)",
    )
    a = ap.parse_args()
    rotulos = tuple(r for r in a.rotulos.split(",") if r)
    for f in desenhar_arvore(a.arvore, a.saida, a.titulo, a.outgroup, a.unidade, rotulos, a.tipo):
        print(f)


if __name__ == "__main__":
    main()
