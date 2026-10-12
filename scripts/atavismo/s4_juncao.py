from __future__ import annotations

import json
import time
from pathlib import Path

import requests

BASE = "https://rest.ensembl.org"
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S = requests.Session()
S.headers.update({"Content-Type": "application/json"})


def get(path, **p):
    for a in range(6):
        try:
            r = S.get(BASE + path, params=p, timeout=120)
        except requests.RequestException:
            time.sleep(3 * (a + 1))
            continue
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(3 * (a + 1))
            continue
        r.raise_for_status()
        return r.json()
    raise RuntimeError(path)


def region(sp, c, s, e, label):
    f = get(f"/overlap/region/{sp}/{c}:{s}-{e}", feature="gene")
    allg = sorted(f, key=lambda g: g["start"])
    print(f"\n=== {label} {c}:{s}-{e} ===")
    for g in allg:
        print(
            f"  {g['start']:>12}-{g['end']:<12} {g['strand']:+d} {g.get('external_name') or g['id']:<24} {g.get('biotype')}"
        )
    return allg


res = {}
res["gallus_pos_SULT1E1"] = region(
    "gallus_gallus", "4", 51_300_000, 53_200_000, "Gallus chr4 distal a SULT1E1"
)
OUT.joinpath("s4_juncao_gallus.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
