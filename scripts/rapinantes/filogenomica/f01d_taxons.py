"""
f01d — tabela de taxons do braco de filogenomica dos rapinantes.

Le as tabelas do `datasets summary genome` (f01*.sh, no WSL) e grava
data/proc/rapinantes/filogenomica/taxons.tsv com: nome usado nos alinhamentos,
especie, acessao, fonte (RefSeq | miniprot), tax_id, clado, papel, nivel da
montagem, N50 de contig e de scaffold, e — para miniprot — o taxon RefSeq cuja
proteina serve de consulta.

    uv run python scripts/rapinantes/filogenomica/f01d_taxons.py
"""

from __future__ import annotations

import csv
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]
AMOSTRAGEM = Path(r"\\wsl.localhost\Ubuntu-24.04\home\manec\rapinantes\amostragem")
SAIDA = RAIZ / "data" / "proc" / "rapinantes" / "filogenomica" / "taxons.tsv"

# (acessao, clado, papel, consulta-miniprot). Ordem: grupo externo -> Neoaves -> Telluraves.
SELECAO: list[tuple[str, str, str, str]] = [
    # Galloanserae: raiz
    ("GCF_016699485.2", "Galliformes: Phasianidae", "raiz", ""),
    ("GCF_047663525.1", "Anseriformes: Anatidae", "raiz (2o Galloanserae)", ""),
    # Neoaves fora de Telluraves
    ("GCF_036013475.1", "Columbiformes: Columbidae", "Columbaves", ""),
    ("GCF_009769465.1", "Musophagiformes: Musophagidae", "Columbaves (Otidimorphae)", ""),
    ("GCF_036417535.1", "Otidiformes: Otididae", "Columbaves (Otidimorphae)", ""),
    ("GCF_017976375.1", "Cuculiformes: Cuculidae", "Columbaves (Otidimorphae)", ""),
    ("GCF_003957555.1", "Apodiformes: Trochilidae", "Strisores (H5)", ""),
    ("GCF_020740795.1", "Apodiformes: Apodidae", "Strisores", ""),
    ("GCF_000700745.2", "Caprimulgiformes: Caprimulgidae", "Strisores (H5: corujas + bacuraus)", ""),
    ("GCF_030867145.1", "Opisthocomiformes: Opisthocomidae", "Neoaves (hoatzin)", ""),
    ("GCF_028858705.1", "Gruiformes: Gruidae", "Gruiformes", ""),
    ("GCF_028500815.1", "Charadriiformes: Laridae", "Charadriiformes", ""),
    ("GCF_000690775.1", "Eurypygiformes: Eurypygidae", "Phaethoquornithes", ""),
    ("GCF_000687285.1", "Phaethontiformes: Phaethontidae", "Phaethoquornithes (P. aethereus GCF_964289735.1 sem dados no NCBI Orthologs)", ""),
    ("GCF_030936135.1", "Gaviiformes: Gaviidae", "Aequornithes", ""),
    ("GCF_965638725.1", "Sphenisciformes: Spheniscidae", "Aequornithes", ""),
    ("GCF_000690835.1", "Procellariiformes: Procellariidae", "Aequornithes", ""),
    ("GCF_034638445.1", "Ciconiiformes: Ciconiidae", "Aequornithes (H6: cathartideos + cegonhas)", ""),
    ("GCF_030463565.1", "Pelecaniformes: Pelecanidae", "Aequornithes", ""),
    # Telluraves — Accipitriformes
    ("GCF_018139145.2", "Accipitriformes: Cathartidae", "rapinante", ""),
    ("GCA_037962945.1", "Accipitriformes: Cathartidae", "rapinante (2o Cathartidae)", "Gymnogyps_californianus"),
    ("GCA_034782735.1", "Accipitriformes: Sagittariidae", "rapinante", "Aquila_chrysaetos"),
    ("GCA_026930845.1", "Accipitriformes: Pandionidae", "rapinante", "Aquila_chrysaetos"),
    ("GCA_026212715.1", "Accipitriformes: Accipitridae: Elaninae", "rapinante", "Aquila_chrysaetos"),
    ("GCA_028022735.1", "Accipitriformes: Accipitridae: Gypaetinae", "rapinante", "Aquila_chrysaetos"),
    ("GCA_055772235.1", "Accipitriformes: Accipitridae: Perninae", "rapinante", "Aquila_chrysaetos"),
    ("GCA_982185325.1", "Accipitriformes: Accipitridae: Aegypiinae", "rapinante", "Aquila_chrysaetos"),
    ("GCA_054553165.1", "Accipitriformes: Accipitridae: Circaetinae", "rapinante", "Aquila_chrysaetos"),
    ("GCF_026419915.1", "Accipitriformes: Accipitridae: Harpiinae", "rapinante", ""),
    ("GCF_900496995.4", "Accipitriformes: Accipitridae: Aquilinae", "rapinante", ""),
    ("GCF_929443795.2", "Accipitriformes: Accipitridae: Accipitrinae", "rapinante", ""),
    ("GCF_947461875.1", "Accipitriformes: Accipitridae: Haliaeetinae", "rapinante", ""),
    ("GCF_964188355.1", "Accipitriformes: Accipitridae: Buteoninae", "rapinante", ""),
    # Falconiformes
    ("GCA_013399355.1", "Falconiformes: Falconidae: Herpetotherinae", "rapinante", "Falco_peregrinus"),
    ("GCA_035590275.1", "Falconiformes: Falconidae: Polyborinae", "rapinante", "Falco_peregrinus"),
    ("GCF_023634155.1", "Falconiformes: Falconidae: Falconinae", "rapinante", ""),
    ("GCF_017639655.2", "Falconiformes: Falconidae: Falconinae", "rapinante", ""),
    # Strigiformes
    ("GCF_018691265.1", "Strigiformes: Tytonidae", "rapinante", ""),
    ("GCF_031877795.1", "Strigiformes: Strigidae: Striginae", "rapinante", ""),
    ("GCF_965140245.1", "Strigiformes: Strigidae: Surniinae", "rapinante", ""),
    ("GCA_982290615.1", "Strigiformes: Strigidae: Ninoxinae", "rapinante", "Strix_aluco"),
    # demais Telluraves
    ("GCF_000690535.1", "Cariamiformes: Cariamidae", "Australaves (H7)", ""),
    ("GCF_004027225.2", "Psittaciformes: Strigopidae", "Australaves", ""),
    ("GCF_000696875.1", "Psittaciformes: Strigopidae (Nestor)", "Australaves", ""),
    ("GCF_012275295.1", "Psittaciformes: Psittaculidae", "Australaves", ""),
    ("GCF_000695815.1", "Passeriformes: Acanthisittidae", "Australaves", ""),
    ("GCF_009829145.1", "Passeriformes: Pipridae (suboscine)", "Australaves", ""),
    ("GCF_048771995.1", "Passeriformes: Estrildidae (oscine)", "Australaves", ""),
    ("GCF_000738735.6", "Passeriformes: Corvidae (oscine)", "Australaves", ""),
    ("GCF_028858725.1", "Coliiformes: Coliidae", "Afroaves (Coraciimorphae)", ""),
    ("GCF_034619465.1", "Leptosomiformes: Leptosomidae", "Afroaves (Coraciimorphae)", ""),
    ("GCF_000703405.1", "Trogoniformes: Trogonidae", "Afroaves (Coraciimorphae)", ""),
    ("GCF_000710305.1", "Bucerotiformes: Bucerotidae", "Afroaves (Coraciimorphae)", ""),
    ("GCF_000691845.1", "Coraciiformes: Meropidae", "Afroaves (Coraciimorphae)", ""),
    ("GCF_014839835.1", "Piciformes: Picidae", "Afroaves (Coraciimorphae)", ""),
    ("GCF_027791375.1", "Piciformes: Indicatoridae", "Afroaves (Coraciimorphae)", ""),
    ("GCF_015220805.1", "Piciformes: Lybiidae", "Afroaves (Coraciimorphae)", ""),
]


def nome_taxon(especie: str) -> str:
    """Genero_especie (sem subespecie): o mesmo nome em todos os locos."""
    return "_".join(especie.split()[:2])


def main() -> None:
    meta: dict[str, dict[str, str]] = {}
    for tsv in sorted(AMOSTRAGEM.glob("*.tsv")):
        with tsv.open(encoding="utf-8") as fh:
            for linha in csv.DictReader(fh, delimiter="\t"):
                meta.setdefault(linha["Assembly Accession"], linha)
    SAIDA.parent.mkdir(parents=True, exist_ok=True)
    campos = [
        "taxon", "especie", "acessao", "fonte", "tax_id", "clado", "papel", "nivel",
        "contig_n50", "scaffold_n50", "anotacao", "consulta_miniprot",
    ]
    nomes: set[str] = set()
    with SAIDA.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for acc, clado, papel, consulta in SELECAO:
            m = meta.get(acc)
            if m is None:
                raise SystemExit(f"{acc} nao esta nas tabelas de {AMOSTRAGEM}")
            taxon = nome_taxon(m["Organism Name"])
            if taxon in nomes:
                raise SystemExit(f"nome repetido: {taxon}")
            nomes.add(taxon)
            fonte = "RefSeq" if acc.startswith("GCF_") else "miniprot"
            w.writerow({
                "taxon": taxon,
                "especie": m["Organism Name"],
                "acessao": acc,
                "fonte": fonte,
                "tax_id": m["Organism Taxonomic ID"],
                "clado": clado,
                "papel": papel,
                "nivel": m["Assembly Level"],
                "contig_n50": m.get("Assembly Stats Contig N50", ""),
                "scaffold_n50": m.get("Assembly Stats Scaffold N50", ""),
                "anotacao": m.get("Annotation Name", "") if fonte == "RefSeq" else "",
                "consulta_miniprot": consulta,
            })
    n_rs = sum(1 for a, *_ in SELECAO if a.startswith("GCF_"))
    print(f"{len(SELECAO)} taxons ({n_rs} RefSeq, {len(SELECAO) - n_rs} miniprot) -> {SAIDA}")


if __name__ == "__main__":
    main()
