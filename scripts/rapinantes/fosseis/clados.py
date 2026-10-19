"""Clados nomeados (conjuntos de terminais atuais) por matriz, a partir de taxa.tsv."""

from __future__ import annotations

import csv
from pathlib import Path

BASE = Path(__file__).resolve().parents[3] / "runs" / "rapinantes" / "fosseis"

CAPRIMULGIFORMES_TRAD = {
    "Steatornithidae",
    "Podargidae",
    "Caprimulgidae",
    "Nyctibiidae",
    "Aegothelidae",
}
FALCONIFORMES_SL = {"Accipitridae", "Pandionidae", "Sagittariidae", "Cathartidae", "Falconidae"}


def ler_taxa(id_: str) -> list[dict]:
    with open(BASE / "matrizes" / id_ / "taxa.tsv", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def familia_de(r: dict) -> str:
    """Familia efetiva: NCBI, ou o proprio nome quando o terminal ja e uma familia."""
    t = r["terminal"]
    if r["familia"]:
        return r["familia"]
    if t.endswith("idae"):
        return t
    return {"Strigidae_slash_Tytonidae": "Strigidae", "Strigiformes": "Strigidae"}.get(t, "")


def clados(id_: str) -> dict[str, frozenset[str]]:
    tx = [r for r in ler_taxa(id_) if r["tipo"] == "atual"]
    por: dict[str, set[str]] = {}
    for r in tx:
        por.setdefault(r["clado"], set()).add(r["terminal"])
    g = {k: set(v) for k, v in por.items() if k not in ("outros", "fossil")}

    def u(*nomes: str) -> set[str]:
        out: set[str] = set()
        for n in nomes:
            out |= g.get(n, set())
        return out

    def tem(*nomes: str) -> bool:
        return all(g.get(n) for n in nomes)

    if tem("Psittaciformes", "Passeriformes"):
        g["Psittacopasserae"] = u("Psittaciformes", "Passeriformes")
    if tem("Falconiformes", "Psittaciformes", "Passeriformes"):
        g["Eufalconimorphae"] = u("Falconiformes", "Psittaciformes", "Passeriformes")
    if tem("Cariamiformes", "Falconiformes", "Psittaciformes", "Passeriformes"):
        g["Australaves"] = u("Cariamiformes", "Falconiformes", "Psittaciformes", "Passeriformes")
    if tem("Accipitriformes", "Strigiformes"):
        g["Hieraves"] = u("Accipitriformes", "Strigiformes")
    if tem("Accipitriformes", "Strigiformes", "Coraciimorphae"):
        g["Afroaves"] = u("Accipitriformes", "Strigiformes", "Coraciimorphae")
    tell = u(
        "Accipitriformes",
        "Strigiformes",
        "Coraciimorphae",
        "Cariamiformes",
        "Falconiformes",
        "Psittaciformes",
        "Passeriformes",
    )
    if len(tell) >= 2:
        g["Telluraves"] = tell
    fam: dict[str, set[str]] = {}
    for r in tx:
        f = familia_de(r)
        if f:
            fam.setdefault(f, set()).add(r["terminal"])
    for f, v in fam.items():
        g.setdefault(f, set(v))
    # grupos "tradicionais" (morfologia classica)
    fsl = {r["terminal"] for r in tx if familia_de(r) in FALCONIFORMES_SL}
    if len({familia_de(r) for r in tx if r["terminal"] in fsl}) >= 2:
        g["Falconiformes_sl"] = fsl
        sem_cat = {r["terminal"] for r in tx if familia_de(r) in FALCONIFORMES_SL - {"Cathartidae"}}
        if sem_cat != fsl and len(sem_cat) >= 2:
            g["Falconiformes_sl_sem_Cathartidae"] = sem_cat
        acc_fal = {r["terminal"] for r in tx if familia_de(r) in {"Accipitridae", "Falconidae"}}
        if len(acc_fal) >= 2 and acc_fal != fsl:
            g["Accipitridae_Falconidae"] = acc_fal
    owls = {r["terminal"] for r in tx if familia_de(r) in {"Strigidae", "Tytonidae"}}
    capr = {r["terminal"] for r in tx if familia_de(r) in CAPRIMULGIFORMES_TRAD}
    if owls and capr:
        g["Strigiformes_Caprimulgiformes"] = owls | capr
    if owls and fsl:
        g["Raptores_trad"] = owls | fsl
    cic = {r["terminal"] for r in tx if familia_de(r) == "Ciconiidae"}
    cat = {r["terminal"] for r in tx if familia_de(r) == "Cathartidae"}
    if cic and cat:
        g["Cathartidae_Ciconiidae"] = cic | cat
    for r in tx:
        g.setdefault(r["terminal"], {r["terminal"]})
    return {k: frozenset(v) for k, v in g.items() if v}


TRADICIONAIS = [
    "Falconiformes_sl",
    "Falconiformes_sl_sem_Cathartidae",
    "Accipitridae_Falconidae",
    "Raptores_trad",
    "Strigiformes_Caprimulgiformes",
    "Cathartidae_Ciconiidae",
]
MOLECULARES = [
    "Accipitriformes",
    "Strigiformes",
    "Hieraves",
    "Afroaves",
    "Coraciimorphae",
    "Psittacopasserae",
    "Eufalconimorphae",
    "Australaves",
    "Telluraves",
]
