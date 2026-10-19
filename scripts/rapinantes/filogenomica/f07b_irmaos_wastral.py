# f07b_irmaos_wastral.py — sister group (and its local PP) of a few focal clades in each wASTRAL tree, rooted on Gallus
# Uso: uv run python f07b_irmaos_wastral.py runs/filogenia/<execucao>/07_arvore_especies/especies.tre [...]
import sys
from pathlib import Path

from Bio import Phylo

FOCOS = {
    "Strisores (Apus+Calypte)": {"Apus_apus", "Calypte_anna"},
    "Leptosomus": {"Leptosomus_discolor"},
    "Falconiformes+Cariama": {
        "Herpetotheres_cachinnans",
        "Ibycter_americanus",
        "Falco_peregrinus",
        "Falco_naumanni",
        "Cariama_cristata",
    },
    "Psittacopasserae": {
        "Strigops_habroptila",
        "Nestor_notabilis",
        "Melopsittacus_undulatus",
        "Acanthisitta_chloris",
        "Taeniopygia_guttata",
        "Corvus_cornix",
        "Chiroxiphia_lanceolata",
    },
}


def resumo(nomes: set[str]) -> str:
    amostra = ", ".join(sorted(nomes)[:5])
    return f"{len(nomes)} terminais: {amostra}{', ...' if len(nomes) > 5 else ''}"


for caminho in map(Path, sys.argv[1:]):
    arvore = Phylo.read(caminho, "newick")  # pyright: ignore[reportPrivateImportUsage]
    arvore.root_with_outgroup("Gallus_gallus")
    print(f"## {caminho.parent.parent.name}")
    for nome, membros in FOCOS.items():
        alvo = (
            arvore.common_ancestor(*membros)
            if len(membros) > 1
            else next(arvore.find_clades(name=next(iter(membros))))
        )
        caminho_raiz = arvore.get_path(alvo)
        pai = caminho_raiz[-2] if len(caminho_raiz) > 1 else arvore.root
        irmaos = [c for c in pai.clades if c is not alvo]
        folhas_irmas = {t.name for c in irmaos for t in c.get_terminals()}
        pp = pai.confidence if pai.confidence is not None else float("nan")
        print(f"  {nome:28s} irmão de [{resumo(folhas_irmas)}]  PP do nó pai={pp:.3f}")
