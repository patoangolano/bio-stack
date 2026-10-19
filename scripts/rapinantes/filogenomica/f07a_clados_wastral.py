# f07a_clados_wastral.py — monophyly and local PP of the raptor-relevant clades in each wASTRAL tree
# Uso: uv run python f07a_clados_wastral.py runs/filogenia/<execucao>/07_arvore_especies/especies.tre [...]
import sys
from pathlib import Path

from Bio import Phylo

ACCIPITRIFORMES = {
    "Sagittarius_serpentarius",
    "Pandion_haliaetus",
    "Elanus_axillaris",
    "Gypaetus_barbatus",
    "Leptodon_cayanensis",
    "Aegypius_monachus",
    "Pithecophaga_jefferyi",
    "Harpia_harpyja",
    "Aquila_chrysaetos",
    "Astur_gentilis",
    "Haliaeetus_albicilla",
    "Buteo_buteo",
}
CATHARTIDAE = {"Gymnogyps_californianus", "Sarcoramphus_papa"}
STRIGIFORMES = {"Tyto_alba", "Strix_aluco", "Athene_noctua", "Ninox_strenua"}
FALCONIFORMES = {
    "Herpetotheres_cachinnans",
    "Ibycter_americanus",
    "Falco_peregrinus",
    "Falco_naumanni",
}
CARIAMA = {"Cariama_cristata"}
LEPTOSOMUS = {"Leptosomus_discolor"}
CORACIIMORPHAE = {
    "Colius_striatus",
    "Apaloderma_vittatum",
    "Buceros_rhinoceros",
    "Merops_nubicus",
    "Dryobates_pubescens",
    "Indicator_indicator",
    "Pogoniulus_pusillus",
}
PSITTACOPASSERAE = {
    "Strigops_habroptila",
    "Nestor_notabilis",
    "Melopsittacus_undulatus",
    "Acanthisitta_chloris",
    "Taeniopygia_guttata",
    "Corvus_cornix",
    "Chiroxiphia_lanceolata",
}

ACCIPITRIMORPHAE = ACCIPITRIFORMES | CATHARTIDAE
CLADOS = {
    "Accipitriformes": ACCIPITRIFORMES,
    "Accipitrimorphae": ACCIPITRIMORPHAE,
    "Strigiformes": STRIGIFORMES,
    "Falconiformes": FALCONIFORMES,
    "Hieraves (Strig+Accipitrimorphae)": STRIGIFORMES | ACCIPITRIMORPHAE,
    "Hieraves+Leptosomus": STRIGIFORMES | ACCIPITRIMORPHAE | LEPTOSOMUS,
    "Afroaves": ACCIPITRIMORPHAE | STRIGIFORMES | LEPTOSOMUS | CORACIIMORPHAE,
    "Coraciimorphae+Leptosomus": CORACIIMORPHAE | LEPTOSOMUS,
    "Strig+Coraciimorphae(+Lepto)": STRIGIFORMES | CORACIIMORPHAE | LEPTOSOMUS,
    "Falconiformes+Cariama": FALCONIFORMES | CARIAMA,
    "Eufalconimorphae": FALCONIFORMES | PSITTACOPASSERAE,
    "Australaves": FALCONIFORMES | PSITTACOPASSERAE | CARIAMA,
    "Psittacopasserae": PSITTACOPASSERAE,
    "Psittacopasserae+Coraciimorphae": PSITTACOPASSERAE | CORACIIMORPHAE,
    "todos rapinantes+Cariama+Lepto": ACCIPITRIMORPHAE
    | STRIGIFORMES
    | FALCONIFORMES
    | CARIAMA
    | LEPTOSOMUS,
    "Telluraves": ACCIPITRIMORPHAE
    | STRIGIFORMES
    | FALCONIFORMES
    | CARIAMA
    | LEPTOSOMUS
    | CORACIIMORPHAE
    | PSITTACOPASSERAE,
}


def relatorio(caminho: Path) -> None:
    arvore = Phylo.read(caminho, "newick")  # pyright: ignore[reportPrivateImportUsage]
    arvore.root_with_outgroup("Gallus_gallus")
    folhas = {t.name for t in arvore.get_terminals()}
    print(f"## {caminho.parent.parent.name}  ({len(folhas)} terminais)")
    for nome, membros in CLADOS.items():
        presentes = membros & folhas
        if len(presentes) < 2:
            print(f"  {nome:36s} n/a")
            continue
        mrca = arvore.common_ancestor(*presentes)
        extras = {t.name for t in mrca.get_terminals()} - presentes
        if extras:
            print(
                f"  {nome:36s} NAO monofilético (+{len(extras)}: {', '.join(sorted(extras)[:4])})"
            )
        else:
            pp = mrca.confidence if mrca.confidence is not None else float("nan")
            bl = mrca.branch_length if mrca.branch_length is not None else float("nan")
            print(f"  {nome:36s} sim  PP={pp:.3f}  ramo={bl:.3f} CU")


for p in sys.argv[1:]:
    relatorio(Path(p))
