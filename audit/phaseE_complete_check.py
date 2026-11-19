"""Complete Phase E audit from cached Crossref + Europe PMC JSONs (including retried camacho2009,
katoh2013, liang2021), plus semantic check of abstracts against in-text claims and lookup of IQ-TREE /
palaeognath phylogeny citations if needed."""
from __future__ import annotations
import json, re, unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "audit" / "evidence" / "bibliography"
BIB = ROOT / "submission_biorxiv" / "references.bib"
TYP = ROOT / "submission_biorxiv" / "manuscript.typ"

def norm(s: str) -> str:
    s = s.replace("\\&", "&")
    s = re.sub(r"<[^>]+>", "", s)
    s = re.sub(r"[{}\\]", "", s)
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.lower().replace("&amp;", "&")
    return re.sub(r"[^a-z0-9]+", " ", s).strip()

def parse_bib(text: str) -> list[dict]:
    entries = []
    for m in re.finditer(r"@(\w+)\{([^,]+),\n(.*?)\n\}\n", text, re.S):
        e = {"type": m.group(1), "key": m.group(2)}
        for fm in re.finditer(r"^\s*(\w+)\s*=\s*\{(.*)\},?\s*$", m.group(3), re.M):
            e[fm.group(1)] = fm.group(2)
        entries.append(e)
    return entries

entries = parse_bib(BIB.read_text(encoding="utf-8"))
results = []
for e in entries:
    k = e["key"]
    cr_p = EV / "crossref" / f"{k}.json"
    ep_p = EV / "epmc" / f"{k}.json"
    cr = None
    if cr_p.exists() and cr_p.stat().st_size > 10:
        try:
            cr = json.loads(cr_p.read_text(encoding="utf-8")).get("message")
        except Exception:
            cr = None
    ep = None
    if ep_p.exists() and ep_p.stat().st_size > 10:
        try:
            rl = json.loads(ep_p.read_text(encoding="utf-8")).get("resultList", {}).get("result", [])
            ep = rl[0] if rl else None
        except Exception:
            ep = None

    issues = []
    # Compare against CR
    cr_title = " ".join(cr.get("title", [])) if cr else ""
    cr_authors = cr.get("author", []) if cr else []
    cr_vol = cr.get("volume", "") if cr else ""
    cr_pg = (cr.get("page") or cr.get("article-number") or "") if cr else ""
    cr_journal = (cr.get("container-title") or [""])[0] if cr else ""
    cr_yrs = []
    if cr:
        for yk in ("published-print", "published-online", "issued", "published"):
            if yk in cr and cr[yk].get("date-parts") and cr[yk]["date-parts"][0][0]:
                cr_yrs.append(str(cr[yk]["date-parts"][0][0]))
    bib_authors = [a.strip() for a in re.split(r"\s+and\s+", e.get("author", "")) if a.strip()]
    if cr:
        if norm(e["title"]) != norm(cr_title):
            issues.append(f"title: bib='{e['title']}' vs CR='{cr_title}'")
        if str(e["year"]) not in cr_yrs:
            issues.append(f"year: bib={e['year']} vs CR={sorted(set(cr_yrs))}")
        if e.get("volume", "") != str(cr_vol):
            issues.append(f"vol: bib='{e.get('volume','')}' vs CR='{cr_vol}'")
        if e.get("pages", "").replace("--", "-") != str(cr_pg):
            issues.append(f"pages: bib='{e.get('pages','')}' vs CR='{cr_pg}'")
        if len(bib_authors) != len(cr_authors):
            issues.append(f"n_authors: bib={len(bib_authors)} vs CR={len(cr_authors)}")
        if cr.get("updated-by"):
            issues.append(f"CR updated-by: {cr['updated-by']}")
    else:
        issues.append("CR missing")
    ep_abs = (ep.get("abstractText") or "") if ep else ""
    ep_vol = (ep.get("journalInfo") or {}).get("volume", "") if ep else ""
    ep_pg = ep.get("pageInfo", "") if ep else ""
    ep_nauth = len((ep.get("authorList") or {}).get("author", [])) if ep else 0
    ep_pt = (ep.get("pubTypeList") or {}).get("pubType", []) if ep else []
    ep_cc = (ep.get("commentCorrectionList") or {}).get("commentCorrection", []) if ep else []
    retract_flags = [x for x in ep_pt if "retract" in x.lower()] + [c for c in ep_cc if "retract" in str(c.get("type","")).lower() or "concern" in str(c.get("type","")).lower()]
    if retract_flags:
        issues.append(f"RETRACTION/CONCERN: {retract_flags}")

    results.append({
        "key": k, "doi": e.get("doi", ""), "year": e.get("year", ""),
        "cr_ok": cr is not None, "ep_ok": ep is not None,
        "bib_nauth": len(bib_authors), "cr_nauth": len(cr_authors), "ep_nauth": ep_nauth,
        "bib_vol": e.get("volume", ""), "cr_vol": cr_vol, "ep_vol": ep_vol,
        "bib_pg": e.get("pages", ""), "cr_pg": cr_pg, "ep_pg": ep_pg,
        "cr_title": cr_title, "cr_journal": cr_journal,
        "ep_pmid": ep.get("pmid", "") if ep else "",
        "ep_abstract_snippet": re.sub(r"<[^>]+>", "", ep_abs)[:260],
        "issues": issues,
    })

(EV / "audited_46_summary.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"Audited {len(results)} entries: CR_ok={sum(r['cr_ok'] for r in results)}/46, EP_ok={sum(r['ep_ok'] for r in results)}/46")
print("\n=== ALL ENTRIES WITH ANY DIFFERENCE ===")
for r in results:
    if r["issues"]:
        print(f"  {r['key']:18s} (CR={r['cr_ok']}, EP={r['ep_ok']}, auth bib/CR/EP={r['bib_nauth']}/{r['cr_nauth']}/{r['ep_nauth']}, vol={r['bib_vol']}/{r['cr_vol']}/{r['ep_vol']}, pg={r['bib_pg']}/{r['cr_pg']}/{r['ep_pg']})")
        for x in r["issues"]:
            print(f"     -> {x}")

print("\n=== CHECK KEY ABSTRACTS (2024-2026 + methodological refs) ===")
for r in results:
    if r["key"] in ("kantaputra2026", "emerling2026", "chun2025", "randall2024", "chow2024", "redlich2024", "kirilenko2023", "weisman2022", "weisman2020", "saputra2021", "prudent2016", "kowalczyk2019", "brocklehurst2021"):
        print(f"  [{r['key']}] ({r['year']}) {r['cr_title']}")
        print(f"     ABS: {r['ep_abstract_snippet']}")
