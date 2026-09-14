"""Conteudo genico do intervalo de Gallus gallus chr4 onde o cluster SCPP deveria estar."""
from __future__ import annotations
import json, time, sys
from pathlib import Path
import requests

BASE = "https://rest.ensembl.org"
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S = requests.Session(); S.headers.update({"Content-Type": "application/json"})

def get(path, **params):
    last = None
    for attempt in range(6):
        try:
            r = S.get(BASE + path, params=params, timeout=120)
        except requests.RequestException as exc:
            last = exc; time.sleep(3 * (attempt + 1)); continue
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(float(r.headers.get("Retry-After", 3 * (attempt + 1)))); continue
        r.raise_for_status(); return r.json()
    raise RuntimeError(f"{path}: {last}")

def region(species, chrom, start, end, label):
    feats = get(f"/overlap/region/{species}/{chrom}:{start}-{end}", feature="gene")
    pc = sorted([g for g in feats if g.get("biotype") == "protein_coding"], key=lambda g: g["start"])
    allg = sorted(feats, key=lambda g: g["start"])
    print(f"\n=== {label}  {chrom}:{start}-{end}  ({len(allg)} genes, {len(pc)} codificantes) ===")
    for g in allg:
        print(f"  {g['start']:>12}-{g['end']:<12} {g['strand']:+d} {g.get('external_name') or g['id']:<22} {g.get('biotype')}")
    return allg

res = {}
# intervalo de Gallus entre DCK/JCHAIN e SULT1B1 (chr4)
res["gallus_chr4_48.6_49.7"] = region("gallus_gallus", "4", 48_600_000, 49_700_000, "Gallus gallus chr4 proximal a JCHAIN")
res["gallus_chr4_49.65_51.4"] = region("gallus_gallus", "4", 49_650_000, 51_400_000, "Gallus gallus chr4 entre DCK e SULT1B1")
OUT.joinpath("s3_intervalo_gallus.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
print("\nOK -> s3_intervalo_gallus.json")
