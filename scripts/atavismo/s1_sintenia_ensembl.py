"""Sintenia do cluster SCPP de esmalte (ENAM/AMBN): vizinhos em Homo e busca
dos ortologos desses vizinhos em Gallus gallus. Ensembl REST."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import requests

BASE = "https://rest.ensembl.org"
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S = requests.Session()
S.headers.update({"Content-Type": "application/json"})


def get(path: str, **params) -> dict | list:
    for attempt in range(5):
        r = S.get(BASE + path, params=params, timeout=60)
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 2)))
            continue
        if r.status_code >= 500:
            time.sleep(2 * (attempt + 1))
            continue
        r.raise_for_status()
        return r.json()
    raise RuntimeError(f"falhou: {path}")


def lookup_symbol(species: str, symbol: str) -> dict | None:
    try:
        return get(f"/lookup/symbol/{species}/{symbol}", expand=0)
    except requests.HTTPError:
        return None


def genes_in_region(species: str, chrom: str, start: int, end: int) -> list[dict]:
    feats = get(f"/overlap/region/{species}/{chrom}:{start}-{end}", feature="gene")
    out = []
    for f in feats:
        out.append(
            {
                "id": f["id"],
                "name": f.get("external_name") or f["id"],
                "biotype": f.get("biotype"),
                "start": f["start"],
                "end": f["end"],
                "strand": f["strand"],
            }
        )
    return sorted(out, key=lambda g: g["start"])


res: dict = {}

# --- 1. ancoras humanas -------------------------------------------------
human = {}
for sym in ["ENAM", "AMBN", "AMTN", "MMP20", "ODAM", "KLK4", "DSPP", "FAM83H", "HES7"]:
    g = lookup_symbol("homo_sapiens", sym)
    human[sym] = (
        None
        if g is None
        else {
            "id": g["id"],
            "chr": g["seq_region_name"],
            "start": g["start"],
            "end": g["end"],
            "strand": g["strand"],
            "biotype": g["biotype"],
        }
    )
    print(f"human {sym}: {human[sym]}", file=sys.stderr)
res["human_anchors"] = human

# --- 2. vizinhanca humana do cluster ENAM/AMBN --------------------------
e, a = human["ENAM"], human["AMBN"]
lo = min(e["start"], a["start"]) - 400_000
hi = max(e["end"], a["end"]) + 400_000
res["human_region"] = {"chr": e["chr"], "start": lo, "end": hi}
res["human_neighbors"] = [
    g for g in genes_in_region("homo_sapiens", e["chr"], lo, hi) if g["biotype"] == "protein_coding"
]
for g in res["human_neighbors"]:
    print(f"  HS {g['start']:>12}-{g['end']:<12} {g['strand']:+d} {g['name']}", file=sys.stderr)

OUT.joinpath("s1_sintenia_humana.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
print("\nOK -> s1_sintenia_humana.json")
