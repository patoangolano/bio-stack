"""Phase E reference audit: query Crossref + Europe PMC for every bib entry, store raw output, compare fields.

Run: uv run python audit/bibliography_audit.py
Outputs: audit/evidence/bibliography/{crossref,epmc}/<key>.json, audit/evidence/bibliography/comparison.json
Network via curl.exe (documented public APIs only).
"""
from __future__ import annotations

import json
import re
import subprocess
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIB = ROOT / "submission_biorxiv" / "references.bib"
EV = ROOT / "audit" / "evidence" / "bibliography"
(EV / "crossref").mkdir(parents=True, exist_ok=True)
(EV / "epmc").mkdir(parents=True, exist_ok=True)
MAILTO = "mailto=bio-stack-audit@example.org"


def curl(url: str) -> tuple[int, str]:
    p = subprocess.run(["curl.exe", "-sS", "-L", "--max-time", "60", "-w", "\n%{http_code}", url],
                       capture_output=True, text=True, encoding="utf-8")
    out = p.stdout
    body, _, code = out.rpartition("\n")
    try:
        return int(code), body
    except ValueError:
        return 0, p.stderr


def parse_bib(text: str) -> list[dict]:
    entries = []
    for m in re.finditer(r"@(\w+)\{([^,]+),\n(.*?)\n\}\n", text, re.S):
        e = {"type": m.group(1), "key": m.group(2)}
        for fm in re.finditer(r"^\s*(\w+)\s*=\s*\{(.*)\},?\s*$", m.group(3), re.M):
            e[fm.group(1)] = fm.group(2)
        entries.append(e)
    return entries


def norm(s: str) -> str:
    s = s.replace("\\&", "&").replace("\\u0026", "&").replace("\\\u0026", "&")
    s = re.sub(r"<[^>]+>", "", s)
    s = re.sub(r"[{}\\]", "", s)
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.lower().replace("&amp;", "&")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return s.strip()


def bib_families(author: str) -> list[str]:
    out = []
    for a in re.split(r"\s+and\s+", author):
        a = a.strip()
        if a.startswith("{"):
            out.append(norm(a))
        elif "," in a:
            out.append(norm(a.split(",")[0]))
        else:
            out.append(norm(a.split()[-1]))
    return out


def main() -> None:
    entries = parse_bib(BIB.read_text(encoding="utf-8"))
    from concurrent.futures import ThreadPoolExecutor

    def work(e):
        key = e["key"]
        r: dict = {"key": key, "bib": e, "issues": [], "notes": []}
        doi = e.get("doi", "")
        # ---- Crossref
        cr = None
        if doi:
            code, body = curl(f"https://api.crossref.org/works/{doi}?{MAILTO}")
            (EV / "crossref" / f"{key}.json").write_text(body, encoding="utf-8")
            r["crossref_http"] = code
            if code == 200:
                try:
                    cr = json.loads(body)["message"]
                except Exception as ex:  # noqa: BLE001
                    r["issues"].append(f"crossref parse error {ex}")
            else:
                r["issues"].append(f"Crossref HTTP {code}")
        else:
            r["issues"].append("no DOI in bib")
        # ---- Europe PMC
        ep = None
        q = None
        if e.get("pmid"):
            q = f"EXT_ID:{e['pmid']}%20AND%20SRC:MED"
        elif doi:
            q = f"DOI:%22{doi}%22"
        if q:
            code, body = curl(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=core")
            (EV / "epmc" / f"{key}.json").write_text(body, encoding="utf-8")
            r["epmc_http"] = code
            if code == 200:
                try:
                    res = json.loads(body)["resultList"]["result"]
                    ep = res[0] if res else None
                    if ep is None:
                        r["notes"].append("EuropePMC: 0 hits")
                except Exception as ex:  # noqa: BLE001
                    r["issues"].append(f"epmc parse error {ex}")
            else:
                r["issues"].append(f"EuropePMC HTTP {code}")

        r["crossref_ok"] = cr is not None
        r["epmc_ok"] = ep is not None
        # ---- compare vs Crossref
        if cr:
            cr_title = norm(" ".join(cr.get("title", [""])))
            if norm(e["title"]) != cr_title:
                r["issues"].append(f"title differs: crossref='{' '.join(cr.get('title', []))}'")
            fam = [norm(a.get("family", a.get("name", ""))) for a in cr.get("author", [])]
            bf = bib_families(e["author"])
            r["crossref_nauthors"] = len(fam)
            r["bib_nauthors"] = len(bf)
            if fam != bf:
                r["issues"].append(f"authors differ: bib n={len(bf)} (first {bf[0]}, last {bf[-1]}) vs crossref n={len(fam)} (first {fam[0] if fam else None}, last {fam[-1] if fam else None})")
            ct = [norm(x) for x in cr.get("container-title", [])]
            if norm(e["journal"]) not in ct:
                r["issues"].append(f"journal differs: crossref={cr.get('container-title')}")
            yrs = []
            for k in ("published-print", "published-online", "issued", "published"):
                if k in cr and cr[k].get("date-parts") and cr[k]["date-parts"][0][0]:
                    yrs.append((k, cr[k]["date-parts"][0][0]))
            r["crossref_years"] = yrs
            if str(e["year"]) not in [str(y) for _, y in yrs]:
                r["issues"].append(f"year differs: crossref={yrs}")
            elif len({y for _, y in yrs}) > 1:
                r["notes"].append(f"crossref year ambiguity (print/online): {yrs}")
            if "volume" in e and e["volume"] != cr.get("volume", ""):
                r["issues"].append(f"volume differs: crossref={cr.get('volume')}")
            if "volume" not in e and cr.get("volume"):
                r["issues"].append(f"volume missing: crossref={cr.get('volume')}")
            if "number" in e and e["number"] != cr.get("issue", ""):
                r["notes"].append(f"issue: bib={e['number']} crossref={cr.get('issue')}")
            pg_bib = e.get("pages", "").replace("--", "-")
            pg_cr = cr.get("page") or cr.get("article-number") or ""
            if pg_bib != pg_cr:
                r["issues"].append(f"pages differ: bib={pg_bib} crossref page={cr.get('page')} article-number={cr.get('article-number')}")
            r["crossref_title"] = cr.get("title")
            r["crossref_container"] = cr.get("container-title")
            r["crossref_page"] = cr.get("page")
            r["crossref_article_number"] = cr.get("article-number")
            r["crossref_volume"] = cr.get("volume")
            r["crossref_issue"] = cr.get("issue")
            r["crossref_type"] = cr.get("type")
            r["crossref_doi"] = cr.get("DOI")
            r["crossref_first_author"] = cr.get("author", [{}])[0].get("family")
            r["crossref_last_author"] = cr.get("author", [{}])[-1].get("family")
            r["crossref_updated_by"] = cr.get("updated-by")
            r["crossref_update_to"] = cr.get("update-to")
            r["crossref_relation"] = {k: [x.get("id") for x in v] for k, v in (cr.get("relation") or {}).items()}
            if cr.get("updated-by"):
                r["issues"].append(f"Crossref updated-by notices: {cr['updated-by']}")
            if cr.get("update-to"):
                r["issues"].append(f"This record is itself an update notice: {cr['update-to']}")
        # ---- compare vs EuropePMC
        if ep:
            r["epmc_title"] = ep.get("title")
            r["epmc_doi"] = ep.get("doi")
            r["epmc_authorString"] = ep.get("authorString")
            r["epmc_journal"] = (ep.get("journalInfo") or {}).get("journal", {}).get("title")
            r["epmc_volume"] = (ep.get("journalInfo") or {}).get("volume")
            r["epmc_issue"] = (ep.get("journalInfo") or {}).get("issue")
            r["epmc_pageInfo"] = ep.get("pageInfo")
            r["epmc_pubYear"] = ep.get("pubYear")
            r["epmc_pubType"] = ep.get("pubTypeList", {}).get("pubType")
            r["epmc_pmid"] = ep.get("pmid")
            r["epmc_id"] = ep.get("id")
            r["epmc_source"] = ep.get("source")
            r["epmc_authorList_n"] = len((ep.get("authorList") or {}).get("author", []))
            r["epmc_commentCorrections"] = ep.get("commentCorrectionList")
            cc = (ep.get("commentCorrectionList") or {}).get("commentCorrection", [])
            for c in cc:
                if c.get("type") in ("Erratum in", "Retraction in", "Retracted", "Expression of concern in", "Corrected and republished in", "Update in", "Correction in"):
                    r["issues"].append(f"EuropePMC notice: {c.get('type')} {c.get('id')} {c.get('source')}")
            pt = r["epmc_pubType"] or []
            if any("Retract" in x for x in pt):
                r["issues"].append(f"EuropePMC pubType retraction: {pt}")
            if norm(ep.get("title", "")) != norm(e["title"]):
                r["issues"].append(f"EuropePMC title differs: '{ep.get('title')}'")
            if e.get("doi") and ep.get("doi") and ep["doi"].lower() != e["doi"].lower():
                r["issues"].append(f"EuropePMC DOI differs: {ep['doi']}")
            if e.get("pmid") and ep.get("pmid") and ep["pmid"] != e["pmid"]:
                r["issues"].append(f"PMID differs: epmc={ep['pmid']}")
            if str(ep.get("pubYear")) != str(e["year"]):
                r["notes"].append(f"EuropePMC pubYear={ep.get('pubYear')} vs bib {e['year']}")
            n_epmc = r["epmc_authorList_n"]
            if n_epmc and n_epmc != r.get("bib_nauthors"):
                r["issues"].append(f"EuropePMC author count {n_epmc} != bib {r.get('bib_nauthors')}")
        print(key, "CR" if cr else "--", "EP" if ep else "--", len(r["issues"]), "issue(s)", flush=True)
        return r

    with ThreadPoolExecutor(6) as ex:
        results = list(ex.map(work, entries))
    (EV / "comparison.json").write_text(json.dumps(results, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()

