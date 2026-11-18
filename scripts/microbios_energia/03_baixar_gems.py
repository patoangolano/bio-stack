"""Baixa modelos metabólicos em escala genômica (BiGG) para a FBA da pesquisa
de microrganismos energéticos. Saída: data/proc/gem/<id>.json"""

from __future__ import annotations

import urllib.request
from pathlib import Path

MODELOS = {
    "iHN637": "Clostridium ljungdahlii DSM 13528 (acetogênico, Wood-Ljungdahl)",
    "iJN678": "Synechocystis sp. PCC 6803 (fotoautotrófico oxigênico)",
    "iAF692": "Methanosarcina barkeri Fusaro (metanogênico versátil)",
}

dest = Path("data/proc/gem")
dest.mkdir(parents=True, exist_ok=True)
for mid, desc in MODELOS.items():
    out = dest / f"{mid}.json"
    if out.exists() and out.stat().st_size > 0:
        print(f"{mid}: já existe ({out.stat().st_size} B)")
        continue
    url = f"http://bigg.ucsd.edu/static/models/{mid}.json"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            out.write_bytes(r.read())
        print(f"{mid}: OK {out.stat().st_size} B — {desc}")
    except Exception as e:  # noqa: BLE001
        print(f"{mid}: FALHOU {e!r}")
