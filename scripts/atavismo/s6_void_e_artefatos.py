from __future__ import annotations

import json
import re
import time
from pathlib import Path

import requests

BASE = "https://rest.ensembl.org"
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S = requests.Session()


def jget(path, **p):
    S.headers.update({"Content-Type": "application/json"})
    for a in range(6):
        try:
            r = S.get(BASE + path, params=p, timeout=180)
        except requests.RequestException:
            time.sleep(3 * (a + 1))
            continue
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(3 * (a + 1))
            continue
        if r.status_code == 400:
            return None
        r.raise_for_status()
        return r.json()
    return None


def seq(species, chrom, s, e):
    r = S.get(
        f"{BASE}/sequence/region/{species}/{chrom}:{s}-{e}",
        headers={"Content-Type": "text/x-fasta"},
        timeout=600,
    )
    r.raise_for_status()
    return r.text


res = {}
# --- A. o vazio de 712 kb e sequencia real ou gap de montagem? ---
fa = seq("gallus_gallus", "4", 51_300_000, 52_030_000)
body = "".join(fa.splitlines()[1:])
n = body.upper().count("N")
runs = [m.span() for m in re.finditer(r"N{100,}", body.upper())]
res["void"] = {
    "chr": "4",
    "start": 51_300_000,
    "end": 52_030_000,
    "len": len(body),
    "n_bases": n,
    "pct_N": round(100 * n / len(body), 3),
    "n_runs_ge100bp": len(runs),
    "maior_run": max((b - a for a, b in runs), default=0),
}
print("VAZIO Gallus chr4:51,300,000-52,030,000 ->", json.dumps(res["void"]))
Path(OUT / "gallus_chr4_void.fa").write_text(fa, encoding="utf-8")

# --- B. artefatos: FAM83H e HES7 ---
for sp, sym in [
    ("homo_sapiens", "FAM83H"),
    ("gallus_gallus", "FAM83H"),
    ("homo_sapiens", "HES7"),
    ("gallus_gallus", "HES7"),
    ("gallus_gallus", "HES5"),
    ("gallus_gallus", "HES6"),
    ("gallus_gallus", "HES4"),
    ("gallus_gallus", "HEY1"),
    ("gallus_gallus", "LFNG"),
    ("gallus_gallus", "DLL1"),
]:
    g = jget(f"/lookup/symbol/{sp}/{sym}")
    res[f"lookup_{sp}_{sym}"] = g
    print(
        f"{sp:16} {sym:8} ->",
        "NENHUM"
        if not g
        else f"{g['id']} chr{g['seq_region_name']}:{g['start']}-{g['end']} {g['biotype']} | {g.get('description')}",
    )

# ortologo de HES7 humano em Gallus, por homologia (nao por simbolo)
for sym in ["HES7", "FAM83H", "ODAM", "KLK4", "AMTN", "ENAM", "AMBN"]:
    h = jget(
        f"/homology/symbol/homo_sapiens/{sym}",
        target_species="gallus_gallus",
        type="orthologues",
        format="condensed",
        sequence="none",
    )
    ids = [
        (x.get("id"), x.get("type"))
        for x in (h["data"][0]["homologies"] if h and h.get("data") else [])
    ]
    res[f"homol_gallus_{sym}"] = ids
    print(f"homologia HS {sym} -> Gallus:", ids or "NENHUM")

OUT.joinpath("s6_void_artefatos.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
