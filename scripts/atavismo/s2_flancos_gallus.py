"""Ortologos em Gallus gallus dos genes que flanqueiam o cluster SCPP de esmalte."""
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


def get(path: str, **params):
    for attempt in range(6):
        r = S.get(BASE + path, params=params, timeout=90)
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 2)));  continue
        if r.status_code >= 500:
            time.sleep(2 * (attempt + 1));  continue
        if r.status_code == 400:
            return None
        r.raise_for_status()
        return r.json()
    return None


# 1) estender a janela humana pro lado proximal (antes de ODAM)
ext = get("/overlap/region/homo_sapiens/4:69300000-70300000", feature="gene")
prox = sorted(
    [g for g in ext if g.get("biotype") == "protein_coding"], key=lambda g: g["start"]
)
print("=== Homo sapiens chr4, lado proximal ao cluster ===", file=sys.stderr)
for g in prox:
    print(f"  {g['start']:>12}-{g['end']:<12} {g['strand']:+d} {g.get('external_name')}", file=sys.stderr)

# 2) ortologos em Gallus dos flancos candidatos
FLANCOS = ["JCHAIN", "UTP3", "RUFY3", "GRSF1", "MOB1B", "DCK", "SULT1B1", "SULT1E1",
           "CSN1S1", "STATH", "ODAM", "AMTN", "AMBN", "ENAM", "MUC7", "CSN3", "FDCSP"]
gal = {}
for sym in FLANCOS:
    h = get(f"/homology/symbol/homo_sapiens/{sym}", target_species="gallus_gallus",
            type="orthologues", sequence="none", format="condensed")
    hits = []
    if h and h.get("data"):
        for hom in h["data"][0].get("homologies", []):
            tid = hom.get("id") or hom["target"]["id"]
            lk = get(f"/lookup/id/{tid}")
            if lk:
                hits.append({
                    "id": tid, "name": lk.get("display_name"),
                    "chr": lk.get("seq_region_name"), "start": lk.get("start"),
                    "end": lk.get("end"), "strand": lk.get("strand"),
                    "type": hom.get("type"), "conf": hom.get("is_high_confidence"),
                })
    gal[sym] = hits
    print(f"GG ortologo de {sym}: {hits if hits else 'NENHUM'}", file=sys.stderr)

OUT.joinpath("s2_flancos_gallus.json").write_text(
    json.dumps({"homo_proximal": prox, "gallus_orthologs": gal}, indent=2), encoding="utf-8")
print("\nOK -> s2_flancos_gallus.json")
