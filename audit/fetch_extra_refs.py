"""Fetch & verify primary metadata (Crossref + Europe PMC) for candidate v1.0.1 references:
- IQ-TREE 2 (Minh et al. 2020, 10.1093/molbev/msaa015) and IQ-TREE 3 (Wong et al. 2025, if published)
- Sackton et al. 2019 (Palaeognathae / ratite paraphyly, 10.1126/science.aat7244)
- Stiller et al. 2024 (avian phylogenomics, 10.1038/s41586-024-07323-1)
- Clopper & Pearson 1934 (exact binomial CI, 10.1093/biomet/26.4.404)
"""
from __future__ import annotations
import json, subprocess
from pathlib import Path

EV = Path(__file__).resolve().parent / "evidence" / "bibliography"
CANDIDATES = {
    "minh2020": "10.1093/molbev/msaa015",
    "sackton2019": "10.1126/science.aat7244",
    "stiller2024": "10.1038/s41586-024-07323-1",
    "clopper1934": "10.1093/biomet/26.4.404",
    "iqtree3_check": "10.32942/X2P62N",  # check if IQ-TREE 3 EcoEvoRxiv / MBE exists via Europe PMC search
}

for key, doi in CANDIDATES.items():
    p = subprocess.run(["curl.exe", "-sS", "-L", "--max-time", "30", "-w", "\n%{http_code}",
                        f"https://api.crossref.org/works/{doi}?mailto=bio-stack-audit@example.org"],
                       capture_output=True, text=True, encoding="utf-8")
    body, _, code = p.stdout.rpartition("\n")
    print(f"{key:15s} ({doi}): CR HTTP {code}", end=" ")
    if code == "200":
        (EV / "crossref" / f"{key}.json").write_text(body, encoding="utf-8")
        m = json.loads(body)["message"]
        yrs = [m[k]["date-parts"][0][0] for k in ("published-print","published-online","issued") if k in m and m[k].get("date-parts")]
        print(f"| {(m.get('title') or [''])[0][:70]} | {(m.get('container-title') or [''])[0]} | yr={yrs} | vol={m.get('volume')} | pg={m.get('page') or m.get('article-number')} | n_auth={len(m.get('author',[]))}")
    else:
        print()

# Search Europe PMC for IQ-TREE 3 paper
p = subprocess.run(["curl.exe", "-sS", "-L", "--max-time", "30",
                    "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22IQ-TREE%203%22&format=json&resultType=core"],
                   capture_output=True, text=True, encoding="utf-8")
if p.returncode == 0 and p.stdout.strip():
    res = json.loads(p.stdout).get("resultList", {}).get("result", [])
    print(f"EPMC IQ-TREE 3 hits: {len(res)}")
    for r in res[:3]:
        print(" ", r.get("pubYear"), r.get("doi"), r.get("title"), r.get("journalTitle"))
