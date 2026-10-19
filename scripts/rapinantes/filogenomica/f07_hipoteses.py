"""
f07 — arvores de restricao e arvores iniciais para os testes de topologia (AU/KH/SH).

Para cada hipotese H: (1) a arvore de restricao do IQ-TREE (-g), multifurcada, com
os clados que H exige; (2) uma arvore inicial compativel, feita por SPR sobre a
arvore ML da supermatriz: o clado X e podado e reenxertado como irmao de Y (a
busca restrita parte dela, -t). Se a arvore ML ja satisfaz H, a melhor arvore
restrita e a propria ML e a busca e pulada (coluna `ja_na_ml`).

Saida em <dir_saida>: H*.restricao.nwk, H*.inicial.nwk, ml.nwk, astral.nwk e
hipoteses.tsv.

    uv run python scripts/rapinantes/filogenomica/f07_hipoteses.py <run_id> <dir_saida>
"""

from __future__ import annotations

import copy
import csv
import io
import re
import sys
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

RAIZ = Path(__file__).resolve().parents[3]

G: dict[str, set[str]] = {
    "CATH": {"Gymnogyps_californianus", "Sarcoramphus_papa"},
    "SAG": {"Sagittarius_serpentarius"},
    "PAN": {"Pandion_haliaetus"},
    "ACCIPITRIDAE": {"Elanus_axillaris", "Gypaetus_barbatus", "Leptodon_cayanensis", "Aegypius_monachus",
                     "Pithecophaga_jefferyi", "Harpia_harpyja", "Aquila_chrysaetos", "Astur_gentilis",
                     "Haliaeetus_albicilla", "Buteo_buteo"},
    "STRIG": {"Tyto_alba", "Strix_aluco", "Athene_noctua", "Ninox_strenua"},
    "CORAC": {"Colius_striatus", "Leptosomus_discolor", "Apaloderma_vittatum", "Buceros_rhinoceros",
              "Merops_nubicus", "Dryobates_pubescens", "Indicator_indicator", "Pogoniulus_pusillus"},
    "FALC": {"Herpetotheres_cachinnans", "Ibycter_americanus", "Falco_peregrinus", "Falco_naumanni"},
    "PSIT": {"Strigops_habroptila", "Nestor_notabilis", "Melopsittacus_undulatus"},
    "PASS": {"Acanthisitta_chloris", "Chiroxiphia_lanceolata", "Taeniopygia_guttata", "Corvus_cornix"},
    "CARI": {"Cariama_cristata"},
    "CAPRI": {"Antrostomus_carolinensis"},
    "CICO": {"Ciconia_boyciana"},
    "RAIZ": {"Gallus_gallus", "Anas_platyrhynchos"},
}
G["ACC"] = G["SAG"] | G["PAN"] | G["ACCIPITRIDAE"]
G["ACCMORPH"] = G["CATH"] | G["ACC"]
G["AFRO"] = G["ACCMORPH"] | G["STRIG"] | G["CORAC"]
G["PSPA"] = G["PSIT"] | G["PASS"]
G["EUFAL"] = G["FALC"] | G["PSPA"]
G["AUST"] = G["CARI"] | G["EUFAL"]

# id, descricao, clados exigidos (restricao), X (podado), Y (novo irmao de X)
HIPOTESES: list[tuple[str, str, list[str], str, str]] = [
    ("H1", "Accipitrimorphae irmao do resto de Afroaves ((Strigiformes,Coraciimorphae),Accipitrimorphae)",
     ["AFRO", "ACCMORPH", "STRIG+CORAC"], "STRIG", "CORAC"),
    ("H1b", "Strigiformes irmao do resto de Afroaves ((Accipitrimorphae,Coraciimorphae),Strigiformes)",
     ["AFRO", "STRIG", "ACCMORPH+CORAC"], "STRIG", "ACCMORPH+CORAC"),
    ("H2", "Strigiformes irmao de Accipitrimorphae", ["AFRO", "STRIG+ACCMORPH"], "STRIG", "ACCMORPH"),
    ("H3", "Strigiformes irmao de Coraciimorphae (Accipitrimorphae livre)", ["STRIG+CORAC"], "STRIG", "CORAC"),
    ("H4", "Falconiformes s.l. (Falconidae + Accipitrimorphae)", ["FALC+ACCMORPH"], "FALC", "ACCMORPH"),
    ("H4b", "Falconidae + Accipitriformes s.s. (sem Cathartidae)", ["FALC+ACC"], "FALC", "ACC"),
    ("H5", "Strigiformes + Caprimulgidae", ["STRIG+CAPRI"], "STRIG", "CAPRI"),
    ("H6", "Cathartidae + Ciconiiformes", ["CATH+CICO"], "CATH", "CICO"),
    ("H7a", "Cariamiformes irmao de Eufalconimorphae", ["AUST", "EUFAL"], "CARI", "EUFAL"),
    ("H7b", "Cariamiformes irmao so de Falconiformes", ["CARI+FALC"], "CARI", "FALC"),
    ("S1", "Sagittarius irmao de (Pandion + Accipitridae)", ["ACC", "PAN+ACCIPITRIDAE"], "SAG", "PAN+ACCIPITRIDAE"),
    ("S2", "Sagittarius irmao de Accipitridae (Pandion fora)", ["ACC", "SAG+ACCIPITRIDAE"], "SAG", "ACCIPITRIDAE"),
    ("S3", "Sagittarius + Pandion", ["ACC", "SAG+PAN"], "SAG", "PAN"),
]


def conjunto(expr: str, taxa: set[str]) -> set[str]:
    s: set[str] = set()
    for parte in expr.split("+"):
        s |= G[parte]
    return s & taxa


def ler_arvore(caminho: Path) -> Tree:
    txt = caminho.read_text(encoding="utf-8").strip()
    txt = re.sub(r":nan", "", txt)
    txt = re.sub(r"\)[^:,;()]+", ")", txt)  # rotulos internos (suporte) fora
    return Phylo.read(io.StringIO(txt), "newick")  # pyright: ignore[reportPrivateImportUsage]


def folhas(c: Clade) -> set[str]:
    return {t.name for t in c.get_terminals()}


def pai(arv: Tree, alvo: Clade) -> Clade:
    for c in arv.find_clades(order="level"):
        if any(f is alvo for f in c.clades):
            return c
    raise ValueError("sem pai")


def mrca_monofiletico(arv: Tree, s: set[str]) -> Clade | None:
    c = arv.common_ancestor(sorted(s))
    return c if folhas(c) == s else None


def enraizar(arv: Tree, taxa: set[str]) -> None:
    fora = sorted(G["RAIZ"] & taxa)
    arv.root_with_outgroup(fora)


def spr(arv: Tree, x: set[str], y: set[str]) -> Tree | None:
    t = copy.deepcopy(arv)
    cx = mrca_monofiletico(t, x)
    if cx is None:
        return None
    p = pai(t, cx)
    p.clades = [c for c in p.clades if c is not cx]
    if len(p.clades) == 1 and p is not t.root:
        so = p.clades[0]
        avo = pai(t, p)
        so.branch_length = (so.branch_length or 0) + (p.branch_length or 0)
        avo.clades = [so if c is p else c for c in avo.clades]
    cy = mrca_monofiletico(t, y)
    if cy is None:
        return None
    if cy is t.root:
        return None
    py = pai(t, cy)
    bl = cy.branch_length or 0.0
    novo = Clade(branch_length=bl / 2, clades=[cy, cx])
    cy.branch_length = bl / 2
    py.clades = [novo if c is cy else c for c in py.clades]
    return t


def satisfaz(arv: Tree, clados: list[set[str]], taxa: set[str]) -> bool:
    bip = set()
    for c in arv.find_clades():
        f = frozenset(folhas(c))
        bip.add(f)
        bip.add(frozenset(taxa - f))
    return all(frozenset(s) in bip or frozenset(taxa - s) in bip for s in clados)


def restricao(clados: list[set[str]], taxa: set[str]) -> str:
    def monta(s: set[str]) -> str:
        filhos = [c for c in clados if c < s]
        maximos = [c for c in filhos if not any(c < d for d in filhos)]
        usados = set().union(*maximos) if maximos else set()
        partes = [monta(c) for c in sorted(maximos, key=lambda c: sorted(c))]
        partes += sorted(s - usados)
        return "(" + ",".join(partes) + ")"
    return monta(taxa) + ";"


def main() -> None:
    run = RAIZ / "runs" / "filogenia" / sys.argv[1]
    saida = Path(sys.argv[2])
    saida.mkdir(parents=True, exist_ok=True)
    ml = ler_arvore(run / "05_concatenado" / "concatenado.treefile")
    taxa = folhas(ml.root)
    enraizar(ml, taxa)
    Phylo.write(ml, saida / "ml.nwk", "newick", plain=False)  # pyright: ignore[reportPrivateImportUsage]
    astral = ler_arvore(run / "07_arvore_especies" / "especies.tre")
    for c in astral.find_clades():
        c.branch_length = None
    Phylo.write(astral, saida / "astral.nwk", "newick", plain=True)  # pyright: ignore[reportPrivateImportUsage]
    linhas = []
    for hid, desc, exig, x, y in HIPOTESES:
        clados = [conjunto(e, taxa) for e in exig]
        clados = [c for c in clados if 1 < len(c) < len(taxa) - 1]
        ja = satisfaz(ml, clados, taxa)
        (saida / f"{hid}.restricao.nwk").write_text(restricao(clados, taxa) + "\n", encoding="utf-8")
        ini = ml if ja else spr(ml, conjunto(x, taxa), conjunto(y, taxa))
        ok_ini = ini is not None and satisfaz(ini, clados, taxa)
        if ini is not None and ok_ini:
            Phylo.write(ini, saida / f"{hid}.inicial.nwk", "newick", plain=False)  # pyright: ignore[reportPrivateImportUsage]
        linhas.append({"hipotese": hid, "descricao": desc, "clados": " | ".join(exig),
                       "ja_na_ml": "sim" if ja else "nao",
                       "inicial_spr": "ok" if ok_ini else "falhou (X ou Y nao monofiletico na ML)"})
    with (saida / "hipoteses.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    for linha in linhas:
        print("\t".join(linha.values()))


if __name__ == "__main__":
    main()
