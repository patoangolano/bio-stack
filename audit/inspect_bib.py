"""Analyze Phase E bibliography comparison.json, retry failed Crossref queries, and audit in-text citations."""
from __future__ import annotations
import json, re, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "audit" / "evidence" / "bibliography"
comp = json.loads((EV / "comparison.json").read_text(encoding="utf-8"))

print(f"Total entries in comparison.json: {len(comp)}")
# Retry any Crossref failures
retry_keys = [r["key"] for r in comp if not r.get("crossref_ok")]
print("Crossref missing/failed:", retry_keys)
for r in comp:
    if not r.get("crossref_ok") and r["bib"].get("doi"):
        doi = r["bib"]["doi"]
        p = subprocess.run(["curl.exe", "-sS", "-L", "--max-time", "60", "-w", "\n%{http_code}",
                            f"https://api.crossref.org/works/{doi}?mailto=bio-stack-audit@example.org"],
                           capture_output=True, text=True, encoding="utf-8")
        body, _, code = p.stdout.rpartition("\n")
        print(f"  Retry Crossref {r['key']} ({doi}): HTTP {code}")
        if code == "200":
            (EV / "crossref" / f"{r['key']}.json").write_text(body, encoding="utf-8")

print("\n=== ENTRIES WITH ISSUES OR NOTES ===")
for r in comp:
    if r["issues"] or r["notes"]:
        print(f"[{r['key']}] CR={r['crossref_ok']} EP={r['epmc_ok']} year={r['bib'].get('year')} doi={r['bib'].get('doi')}")
        for x in r["issues"]:
            print(f"   ISSUE: {x}")
        for x in r["notes"]:
            print(f"   NOTE : {x}")

# Check in-text citations in manuscript.typ vs references.bib
typ = (ROOT / "submission_biorxiv" / "manuscript.typ").read_text(encoding="utf-8")
bib_keys = {r["key"] for r in comp}
# Match @key (not @tab1, @fig1) and #pcite(<key>)
at_cites = set(re.findall(r"@([a-z]+[0-9]{4}[a-z]?)", typ))
pcites = set(re.findall(r"#pcite\(<([a-z]+[0-9]{4}[a-z]?)>\)", typ))
used_keys = at_cites | pcites
print(f"\nBib keys: {len(bib_keys)} | Cited in manuscript.typ: {len(used_keys)}")
print("Cited but NOT in bib:", sorted(used_keys - bib_keys))
print("In bib but NOT cited:", sorted(bib_keys - used_keys))

# Print 2025/2026 entries in full detail
print("\n=== 2025/2026 ENTRIES DETAIL ===")
for r in comp:
    if str(r["bib"].get("year")) in ("2025", "2026"):
        print(r["key"], "| bib:", r["bib"])
        print("  CR title:", r.get("crossref_title"), "years:", r.get("crossref_years"), "vol:", r.get("crossref_volume"), "pg:", r.get("crossref_page"), "art:", r.get("crossref_article_number"))
        print("  EP title:", r.get("epmc_title"), "pubYear:", r.get("epmc_pubYear"), "journal:", r.get("epmc_journal"), "pmid:", r.get("epmc_pmid"))
