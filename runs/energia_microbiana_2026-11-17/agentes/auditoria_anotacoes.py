"""Bounded public UniProt/KEGG annotation audit; no genomes or phenotype predictions."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import time
from pathlib import Path
from typing import Any

import requests

OUT = Path(__file__).resolve().parent
JSON_PATH = OUT / "anotacoes.json"
SELECTION = [
    ("CGA009", "O83005", "Fotossistema", "PufL"),
    ("CGA009", "A0AAE9XYU2", "Fotossistema", "PuhA"),
    ("CGA009", "Q6N9K0", "Bacterioclorofila", "BchL"),
    ("CGA009", "Q6N9J0", "Fixação de carbono", "CbbL / RuBisCO forma I"),
    ("CGA009", "Q6NB67", "Oxidação de H2 candidata", "Uptake hydrogenase grande"),
    ("CGA009", "Q6NB68", "Oxidação de H2 candidata", "Uptake hydrogenase pequena"),
    ("CGA009", "A0AAE9XUC8", "Transferência extracelular candidata", "PioA-like"),
    ("CGA009", "Q6NBT3", "Transferência extracelular candidata", "PioB-like"),
    ("CGA009", "Q6NBT4", "Transferência extracelular candidata", "PioC-like"),
    ("TIE-1", "B3Q7E1", "Fixação de carbono", "CbbL / RuBisCO forma I"),
    ("TIE-1", "B3Q7C4", "Bacterioclorofila", "BchN"),
    ("H16", "P31891", "Oxidação de H2", "HoxG / hidrogenase de membrana"),
    ("H16", "P31892", "Oxidação de H2", "HoxK / hidrogenase de membrana"),
    ("H16", "P22320", "Redução de NAD por H2", "HoxH / hidrogenase solúvel"),
    ("H16", "P22317", "Redução de NAD por H2", "HoxF / hidrogenase solúvel"),
    ("H16", "Q0K1E0", "Fixação de carbono", "CbbL1"),
    ("H16", "Q0KEI3", "CODH aeróbica apenas homóloga", "CoxL1-like"),
    ("ATCC11170", "Q2RQ25", "Fotossistema", "PufL-like"),
    ("ATCC11170", "Q2RWS4", "Fotossistema", "PuhA-like"),
    ("ATCC11170", "Q2RWS1", "Bacterioclorofila", "BchL"),
    ("ATCC11170", "Q2RRP5", "Fixação de carbono", "CbbM / RuBisCO forma II"),
    ("ATCC11170", "Q2RV82", "Oxidação de H2 candidata", "Uptake hydrogenase grande"),
    ("ATCC11170", "Q2RUG7", "CODH anaeróbica", "CooS-like / Ni-CODH"),
    ("DSM2944", "A0A1I5IP77", "Sox", "SoxA"),
    ("DSM2944", "A0A1I5INT1", "Sox", "SoxY"),
    ("DSM2944", "A0AAE6TS95", "Sox", "SoxB"),
    ("DSM2944", "A0A1I5IPK4", "Sox", "SoxC"),
    ("DSM2944", "A0A1I5INZ2", "Sox", "SoxD"),
    ("DSM2944", "A0A1I5J2X8", "Fixação de carbono", "CbbL"),
    ("DSM2944", "A0AAE6TRU3", "Oxidação de H2 candidata", "Uptake hydrogenase grande"),
    ("PD1222", "A1B9M4", "Sox", "SoxB"),
    ("PD1222", "A1B2Q2", "Fixação de carbono", "CbbL"),
    ("PD1222", "A1B6N5", "Oxidação de H2 candidata", "Uptake hydrogenase grande"),
    ("PD1222", "A1B9T9", "Desnitrificação, apenas uma etapa", "NosZ"),
    ("J-10-fl", "P11695", "Fotossistema", "PufL"),
    ("J-10-fl", "Q9F6X3", "Bacterioclorofila", "BchL"),
    ("J-10-fl", "A9WIU3", "3HP bi-cycle candidato", "Mcr segundo KEGG K14468"),
    ("J-10-fl", "A9WEU4", "3HP bi-cycle candidato", "Pcs segundo KEGG K14469"),
    ("SCM1", "A9A4U2", "Oxidação de amônia candidata", "AmoA-like"),
    ("SCM1", "A9A4U5", "Oxidação de amônia candidata", "AmoB-like"),
    ("SCM1", "A9A4U4", "Oxidação de amônia candidata", "AmoC-like"),
    ("SCM1", "A9A1Y1", "HP/HB, etapa experimental", "4HB-CoA ligase ADP-forming"),
]
KEGG_SELECTION = [
    ("Rpal_1716", "Fotossistema", "PufL / K08928"),
    ("Rpal_1154", "Oxidação de H2 candidata", "Hidrogenase pequena / K06282"),
    ("Rpal_1155", "Oxidação de H2 candidata", "Hidrogenase grande / K06281"),
    ("Rpal_0817", "Transferência extracelular", "PioA / K28459"),
    ("Rpal_0816", "Transferência extracelular", "PioB / K28460"),
    ("Rpal_0815", "Transferência extracelular", "PioC / K28461"),
]
TARGETS = [
    ("CGA009", 258594, "UP000001426"),
    ("TIE-1", 395960, "UP000001725"),
    ("H16", 381666, "UP000008210"),
    ("ATCC11170", 269796, "UP000001929"),
    ("DSM2944", 82367, "UP000273626"),
    ("PD1222", 318586, "UP000000361"),
    ("J-10-fl", 324602, "UP000002008"),
    ("SCM1", 436308, "UP000000792"),
]
MODULE_QUERY = (
    "gene:pufL OR gene:puhA OR gene:bchL OR gene:cbbL* OR gene:cbbM "
    'OR protein_name:"reaction center" OR protein_name:"protochlorophyllide reductase" '
    'OR protein_name:hydrogenase OR protein_name:"carbon monoxide" '
    'OR protein_name:Sox OR protein_name:"sulfur oxidation" '
    'OR protein_name:"thiosulfate-binding" OR protein_name:"thiosulfate oxidation" '
    'OR protein_name:"sulfur-oxidizing" OR protein_name:"ammonia monooxygenase" '
    'OR protein_name:"malonyl-CoA reductase" OR protein_name:"propionyl-CoA synthase" '
    'OR protein_name:"4-hydroxybutyr" OR gene:amoA OR gene:amoB OR gene:amoC '
    "OR gene:cooS OR gene:cooH OR gene:pioA OR gene:pioB OR gene:pioC "
    "OR gene:nosZ"
)


def evidence_nodes(value: Any, path: str = "") -> list[dict[str, Any]]:
    """Extract evidence objects with their exact annotation-field path."""
    result = []
    if isinstance(value, dict):
        for key, item in value.items():
            field_path = f"{path}.{key}" if path else key
            if key == "evidences":
                result.extend({"field": path, **e} for e in item)
            else:
                result.extend(evidence_nodes(item, field_path))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            result.extend(evidence_nodes(item, f"{path}[{i}]"))
    return result


def summarize_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Keep identifiers, provenance and claim-specific evidence, not sequence."""
    comments = [
        c
        for c in entry.get("comments", [])
        if c.get("commentType")
        in {
            "FUNCTION",
            "CATALYTIC ACTIVITY",
            "COFACTOR",
            "PATHWAY",
            "CAUTION",
            "SIMILARITY",
            "SUBCELLULAR LOCATION",
        }
    ]
    description = entry.get("proteinDescription", {})
    full = description.get("recommendedName") or next(
        iter(description.get("submissionNames", [])), {}
    )
    selected = {
        "proteinDescription": description,
        "comments": comments,
        "genes": entry.get("genes", []),
    }
    refs = []
    for ref in entry.get("references", []):
        citation = ref.get("citation", {})
        refs.append(
            {
                "citation": citation,
                "referencePositions": ref.get("referencePositions", []),
                "referenceNumber": ref.get("referenceNumber"),
                "referenceComments": ref.get("referenceComments", []),
            }
        )
    return {
        "accession": entry["primaryAccession"],
        "entry_id": entry.get("uniProtkbId"),
        "source_database": "UniProtKB",
        "source_url": f"https://rest.uniprot.org/uniprotkb/{entry['primaryAccession']}.json",
        "entry_type": entry.get("entryType"),
        "reviewed": entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
        "protein_name": full.get("fullName", {}).get("value"),
        "organism": entry.get("organism"),
        "genes": entry.get("genes", []),
        "protein_existence": entry.get("proteinExistence"),
        "entry_audit": entry.get("entryAudit"),
        "length_aa": entry.get("sequence", {}).get("length"),
        "annotation_evidence": evidence_nodes(selected),
        "protein_description": description,
        "comments": comments,
        "references": refs,
        "cross_references": [
            x
            for x in entry.get("uniProtKBCrossReferences", [])
            if x.get("database") in {"Proteomes", "KEGG", "EMBL", "RefSeq", "GO", "InterPro"}
        ],
    }


class Client:
    """Small synchronous API client with logged URLs and a hard request budget."""

    def __init__(self, log: list[dict[str, Any]]) -> None:
        self.log = log
        self.session = requests.Session()
        self.session.headers.update(
            {"User-Agent": "bio-stack-academic-functional-annotation-audit/1.0"}
        )
        self.calls = 0

    def get(self, url: str, params: dict[str, Any] | None = None) -> Any:
        """Retrieve public data; retain errors without converting them to absence."""
        self.calls += 1
        if self.calls > 55:
            raise RuntimeError("Local 55-call budget exhausted.")
        rec: dict[str, Any] = {"requested_url": url, "params": params}
        try:
            r = self.session.get(url, params=params, timeout=50)
            rec.update(
                {
                    "url": r.url,
                    "http_status": r.status_code,
                    "accessed_utc": dt.datetime.now(dt.UTC).isoformat(),
                    "uniprot_release": r.headers.get("X-UniProt-Release"),
                    "uniprot_release_date": r.headers.get("X-UniProt-Release-Date"),
                    "total_results": r.headers.get("X-Total-Results"),
                    "next_link": r.links.get("next", {}).get("url"),
                    "bytes": len(r.content),
                }
            )
            r.raise_for_status()
            data = r.json() if "json" in r.headers.get("Content-Type", "") else r.text
            rec["returned_results"] = (
                len(data.get("results", [])) if isinstance(data, dict) else None
            )
            self.log.append(rec)
            return data
        except (requests.RequestException, ValueError) as exc:
            rec["error"] = str(exc)
            self.log.append(rec)
            return {"audit_error": str(exc)}
        finally:
            time.sleep(0.38)


def collect() -> None:
    """Resolve exact proteomes and gather targeted annotation candidates."""
    data: dict[str, Any] = {
        "audit_date": "2026-11-17",
        "scope": "Targeted functional annotation audit, not exhaustive pathway completeness.",
        "primary_source": "UniProt official REST API; KEGG only for coverage gaps",
        "request_log": [],
        "taxa": [],
        "candidate_records": [],
        "errors_are_not_absences": True,
    }
    client = Client(data["request_log"])
    for label, taxon, proteome in TARGETS:
        p = client.get(f"https://rest.uniprot.org/proteomes/{proteome}.json")
        data["taxa"].append(
            {"target": label, "expected_taxon": taxon, "proteome_id": proteome, "record": p}
        )
        # DSM2944 is species-level in taxonomy; exact proteome selection is essential.
        scope = f"proteome:{proteome}" if label == "DSM2944" else f"organism_id:{taxon}"
        query = f"({scope}) AND ({MODULE_QUERY})"
        url: str | None = "https://rest.uniprot.org/uniprotkb/search"
        params: dict[str, Any] | None = {
            "query": query,
            "format": "json",
            "size": 200,
            "sort": "accession asc",
        }
        fetched = []
        while url:
            answer = client.get(url, params)
            if not isinstance(answer, dict) or "audit_error" in answer:
                break
            fetched.extend(answer.get("results", []))
            url = client.log[-1].get("next_link")
            params = None
        for e in fetched:
            entry = summarize_entry(e)
            entry["target"] = label
            data["candidate_records"].append(entry)
        print(f"{label}: {len(fetched)} candidates", flush=True)
    JSON_PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {JSON_PATH}", flush=True)


def inspect() -> None:
    """Print concise candidate names for manual, scientifically bounded selection."""
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    for target, _, _ in TARGETS:
        print(f"\n## {target}")
        for e in data.get("candidate_records", data.get("representatives", [])):
            if e["target"] == target:
                names = []
                for g in e["genes"]:
                    if "geneName" in g:
                        names.append(g["geneName"]["value"])
                    for key in ("orderedLocusNames", "orfNames"):
                        names.extend(n["value"] for n in g.get(key, []))
                print(
                    e["accession"],
                    "R" if e["reviewed"] else "U",
                    ",".join(names),
                    e["protein_name"],
                    sep="\t",
                )


def parse_kegg(text: str) -> list[dict[str, Any]]:
    """Parse flat-file metadata without retaining DNA or protein sequences."""
    result = []
    for block in text.split("///"):
        fields: dict[str, list[str]] = {}
        current = ""
        for line in block.splitlines():
            key, value = line[:12].strip(), line[12:].strip()
            if key:
                current = key
            if current and current not in {"AASEQ", "NTSEQ"}:
                fields.setdefault(current, []).append(value)
        if "ENTRY" in fields:
            result.append(fields)
    return result


def supplement() -> None:
    """Fill functional-name/strain-coverage gaps with bounded targeted queries."""
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    client = Client(data["request_log"])
    queries = [
        (
            "neighbors_sox",
            "((proteome:UP000273626) AND (gene:BDE18_340* OR gene:BDE18_341*)) "
            "OR ((organism_id:318586) AND gene:Pden_415*)",
        ),
        (
            "chlor_carbon",
            "(taxonomy_id:1108) AND (protein_name:malonyl* OR "
            "protein_name:propionyl* OR gene:mcr OR gene:pcs)",
        ),
        ("SCM1", "accession:A9A1Y1"),
    ]
    for tag, query in queries:
        ans = client.get(
            "https://rest.uniprot.org/uniprotkb/search",
            {"query": query, "format": "json", "size": 100, "sort": "accession asc"},
        )
        for raw in ans.get("results", []):
            entry = summarize_entry(raw)
            taxon = raw.get("organism", {}).get("taxonId")
            entry["target"] = (
                "DSM2944"
                if taxon == 82367
                else "PD1222"
                if taxon == 318586
                else "J-10-fl"
                if tag == "chlor_carbon"
                else tag
            )
            if entry["accession"] not in {x["accession"] for x in data["candidate_records"]}:
                data["candidate_records"].append(entry)
        print(f"{tag}: {len(ans.get('results', []))} hits", flush=True)
    entries = "rpt:Rpal_0815+rpt:Rpal_0816+rpt:Rpal_0817+rpt:Rpal_1715+rpt:Rpal_1154+rpt:Rpal_1155"
    ans = client.get(f"https://rest.kegg.jp/get/{entries}")
    if isinstance(ans, str):
        data["kegg_records"] = parse_kegg(ans)
        print("KEGG:", len(data["kegg_records"]), "records", flush=True)
    else:
        data["kegg_error"] = ans
    JSON_PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def show_evidence(accessions: list[str]) -> None:
    """Inspect claim evidence, not whole biological files."""
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    for e in data.get("candidate_records", data.get("representatives", [])):
        if not accessions or e["accession"] in accessions:
            print(
                json.dumps(
                    {
                        k: e.get(k)
                        for k in (
                            "accession",
                            "protein_name",
                            "organism",
                            "protein_existence",
                            "protein_description",
                            "comments",
                            "references",
                            "cross_references",
                        )
                    },
                    ensure_ascii=False,
                    indent=2,
                )
            )
    if accessions == ["KEGG"]:
        print(json.dumps(data.get("kegg_records"), ensure_ascii=False, indent=2))


def resolve_gaps() -> None:
    """Resolve exact locus names and retain conflicting or wrong-strain annotations."""
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    client = Client(data["request_log"])
    query = (
        "accession:A1B9M4 OR accession:A9WIU3 OR accession:A9WEU4 OR "
        "accession:A0A1I5IPK4 OR accession:A0A1I5INZ2 OR BDE18_3408"
    )
    ans = client.get(
        "https://rest.uniprot.org/uniprotkb/search",
        {"query": query, "format": "json", "size": 100, "sort": "accession asc"},
    )
    for raw in ans.get("results", []):
        entry = summarize_entry(raw)
        taxon = raw.get("organism", {}).get("taxonId")
        entry["target"] = {82367: "DSM2944", 318586: "PD1222", 324602: "J-10-fl"}[taxon]
        if entry["accession"] not in {x["accession"] for x in data["candidate_records"]}:
            data["candidate_records"].append(entry)
        print(entry["accession"], entry["protein_name"], flush=True)
    entries = "rpt:Rpal_1716+rpt:Rpal_1717+cau:Caur_2614+cau:Caur_0613+rru:Rru_A1427"
    ans = client.get(f"https://rest.kegg.jp/get/{entries}")
    if isinstance(ans, str):
        data.setdefault("kegg_records", []).extend(parse_kegg(ans))
        print("KEGG extra:", len(parse_kegg(ans)), "records", flush=True)
    else:
        data["kegg_gap_error"] = ans
    JSON_PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def pio_compare() -> None:
    """Check both R. palustris strains; do not turn TIE-1 precedent into exclusivity."""
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    client = Client(data["request_log"])
    query = "(proteome:UP000001426) AND (TX73_003855 OR TX73_003850 OR TX73_003845)"
    ans = client.get(
        "https://rest.uniprot.org/uniprotkb/search",
        {"query": query, "format": "json", "size": 100, "sort": "accession asc"},
    )
    for raw in ans.get("results", []):
        entry = summarize_entry(raw)
        entry["target"] = "CGA009"
        data["candidate_records"].append(entry)
        print(entry["accession"], entry["protein_name"], flush=True)
    entries = "ko:K28459+ko:K28460+ko:K28461+rpa:TX73_003855+rpa:TX73_003850+rpa:TX73_003845"
    ans = client.get(f"https://rest.kegg.jp/get/{entries}")
    if isinstance(ans, str):
        data.setdefault("kegg_records", []).extend(parse_kegg(ans))
        print("Pio KEGG:", len(parse_kegg(ans)), "records", flush=True)
    else:
        data["kegg_pio_error"] = ans
    JSON_PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def gene_values(entry: dict[str, Any]) -> list[str]:
    """Return named genes and locus tags while preserving annotation spelling."""
    result = []
    for gene in entry.get("genes", []):
        if gene.get("geneName"):
            result.append(gene["geneName"]["value"])
        for key in ("orderedLocusNames", "orfNames"):
            result.extend(x["value"] for x in gene.get(key, []))
    return result


def claim_codes(entry: dict[str, Any]) -> tuple[list[str], list[str]]:
    """Separate naming evidence from FUNCTION/CATALYTIC ACTIVITY evidence."""
    name_codes = {e["evidenceCode"] for e in evidence_nodes(entry.get("protein_description", {}))}
    functional = [
        c
        for c in entry.get("comments", [])
        if c.get("commentType") in {"FUNCTION", "CATALYTIC ACTIVITY"}
    ]
    function_codes = {e["evidenceCode"] for e in evidence_nodes(functional)}
    return sorted(name_codes), sorted(function_codes)


def enrich_representative(entry: dict[str, Any]) -> None:
    """Add transparent audit labels without changing the source annotation."""
    name_codes, function_codes = claim_codes(entry)
    entry["name_evidence_codes"] = name_codes
    entry["function_catalysis_evidence_codes"] = function_codes
    entry["evidence_level"] = (
        "A3" if "ECO:0000269" in function_codes else "A2" if entry["reviewed"] else "A1"
    )
    explicit = any(x["database"] == "Proteomes" for x in entry.get("cross_references", []))
    entry["strain_link"] = (
        "Exact taxon plus explicit proteome cross-reference"
        if explicit
        else "Exact strain taxon; no current proteome cross-reference in this record"
    )
    entry["inference_limit"] = (
        "An annotated representative is not a complete pathway, expression, "
        "growth phenotype or organism-level energetic efficiency measurement."
    )


def finalize() -> None:
    """Retain 48 representatives and two cautions; produce the audit report."""
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    candidates = {x["accession"]: x for x in data.get("candidate_records", [])}
    candidate_total = len(candidates) or data.get("counts", {}).get(
        "candidate_records_before_selection"
    )
    if not candidates:
        candidates = {x["accession"]: x for x in data.get("representatives", [])}
        candidates.update({x["accession"]: x for x in data.get("excluded_cautions", [])})
    kegg = {x["ENTRY"][0].split()[0]: x for x in data.get("kegg_records", [])}
    reps = []
    for target, accession, module, label in SELECTION:
        entry = candidates[accession]
        if entry["target"] != target:
            raise ValueError(f"Target mismatch: {accession}")
        expected_taxon = next(t for n, t, _ in TARGETS if n == target)
        if entry["organism"]["taxonId"] != expected_taxon:
            raise ValueError(f"Taxonomy mismatch: {accession}")
        entry.update({"module": module, "audit_label": label})
        enrich_representative(entry)
        entry.pop("kegg_crosschecks", None)
        for locus in gene_values(entry):
            if locus in kegg:
                entry.setdefault("kegg_crosschecks", []).append(kegg[locus])
        reps.append(entry)
    for locus, module, label in KEGG_SELECTION:
        raw = kegg[locus]
        if "TIE-1" not in raw["ORGANISM"][0]:
            raise ValueError(f"KEGG strain mismatch: {locus}")
        reps.append(
            {
                "accession": f"rpt:{locus}",
                "target": "TIE-1",
                "source_database": "KEGG Genes",
                "source_url": f"https://rest.kegg.jp/get/rpt:{locus}",
                "organism": {
                    "scientificName": raw["ORGANISM"][0],
                    "taxonId": 395960,
                    "taxon_id_source": "Exact-strain UniProt proteome UP000001725",
                },
                "module": module,
                "audit_label": label,
                "protein_name": raw["NAME"][0],
                "reviewed": None,
                "evidence_level": "K1",
                "name_evidence_codes": [],
                "function_catalysis_evidence_codes": [],
                "genes": [{"orderedLocusNames": [{"value": locus}]}],
                "protein_existence": None,
                "kegg_record": raw,
                "inference_limit": "KO assignment; KEGG does not supply UniProt ECO or reviewed status.",
            }
        )
    assert len(reps) == 48 and len({x["accession"] for x in reps}) == 48
    cautions = []
    for accession, reason in [
        (
            "Q6QQP7",
            "WRONG STRAIN: referenceComments STRAIN=OK-70-fl, species taxon 1108; "
            "must not be used as native J-10-fl Mcr evidence.",
        ),
        (
            "Q0K5E4",
            "AMBIGUOUS ANNOTATION: putative ammonia monooxygenase is an unreviewed "
            "single-protein label; not evidence for an AmoABC nitrification module in H16.",
        ),
    ]:
        entry = candidates[accession]
        entry["reason_excluded_from_module_confirmation"] = reason
        cautions.append(entry)
    data["counts"] = {
        "candidate_records_before_selection": candidate_total,
        "module_representatives": len(reps),
        "module_representatives_uniprot": len(SELECTION),
        "module_representatives_kegg": len(KEGG_SELECTION),
        "excluded_cautions": len(cautions),
        "script_http_requests": len(data["request_log"]),
    }
    # The initial-hit count is derived from the eight initial logged searches.
    initial_searches = [
        r for r in data["request_log"][:16] if r["requested_url"].endswith("/uniprotkb/search")
    ]
    data["counts"]["initial_targeted_uniprot_hits"] = sum(
        int(r.get("total_results") or 0) for r in initial_searches
    )
    data["representatives"] = sorted(
        reps, key=lambda x: ([n for n, _, _ in TARGETS].index(x["target"]), x["accession"])
    )
    data["excluded_cautions"] = cautions
    data["evidence_levels"] = {
        "A1": "Unreviewed imported/model-based annotation; no functional validation implied.",
        "A2": "Reviewed annotation, without ECO:0000269 on FUNCTION/CATALYTIC ACTIVITY in the fetched fields.",
        "A3": "ECO:0000269 on FUNCTION/CATALYTIC ACTIVITY, scoped to this protein/reaction.",
        "K1": "KEGG exact-strain gene and KO assignment; not a UniProt evidence code.",
    }
    data["interpretation"] = [
        "Both CGA009 and TIE-1 have PioA/B/C candidate annotations; TIE-1 precedent is not gene exclusivity.",
        "H16 CoxL-like homolog is not equivalent to the Ni-CODH candidate of R. rubrum.",
        "Chloroflexus J-10-fl Mcr/Pcs assignments need KEGG cross-check because UniProt names are generic/different.",
        "SCM1 A9A1Y1 has experimental reaction evidence for ATP-to-ADP 4HB-CoA activation, not whole-organism efficiency.",
        "No phenotype, pathway-completeness score or energy-efficiency percentage was computed.",
    ]
    data["validation"] = {
        "all_48_representatives_have_exact_target_taxon": all(
            e["organism"]["taxonId"] == next(t for n, t, _ in TARGETS if n == e["target"])
            for e in reps
        ),
        "all_script_requests_http_200": all(
            r.get("http_status") == 200 for r in data["request_log"]
        ),
        "all_logged_searches_fully_retrieved": all(
            r.get("next_link") is None and int(r["total_results"]) == r["returned_results"]
            for r in data["request_log"]
            if r.get("total_results") is not None
        ),
        "representatives_reviewed_uniprot": sum(e.get("reviewed") is True for e in reps),
        "representatives_unreviewed_uniprot": sum(e.get("reviewed") is False for e in reps),
        "representatives_function_catalysis_ECO_0000269": sum(
            "ECO:0000269" in e.get("function_catalysis_evidence_codes", []) for e in reps
        ),
        "unique_representatives_plus_cautions": len({e["accession"] for e in reps + cautions}),
        "not_claimed": [
            "complete pathway",
            "biological absence",
            "expression",
            "organism-level efficiency",
            "evolutionary loss",
        ],
    }
    data.pop("candidate_records", None)
    JSON_PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    write_markdown(data)
    print(json.dumps(data["counts"], indent=2), flush=True)


def write_markdown(data: dict[str, Any]) -> None:
    """Write the bounded audit with exact accessions and claim-specific ECO fields."""
    lines = [
        "# Auditoria funcional genômica da shortlist energética",
        "",
        "**Consulta: 17/11/2026. Escopo:** representantes de módulos, não inventário "
        "exaustivo de genes, não revisão ampla de literatura e não comparação de eficiência.",
        "",
        "## Conclusões principais",
        "",
        "1. **Pio não é exclusivo de TIE-1.** KEGG atribui PioA/B/C tanto à CGA009 "
        "(TX73_003855/003850/003845) quanto à TIE-1 (Rpal_0817/0816/0815). "
        "Os acessos UniProt da CGA009 são A0AAE9XUC8, Q6NBT3 e Q6NBT4. "
        "Isso é evidência de anotação candidata em CGA009, não demonstração de "
        "fotoferrotrofia ou eletrotrofia nessa cepa. Os próprios KOs citam trabalhos "
        "funcionais de TIE-1: PMID 17189359 e 31690680; não os transferir entre cepas.",
        "2. **H16 tem sistemas de H2 distintos.** HoxG/K, P31891/P31892, são "
        "subunidades da hidrogenase de membrana; HoxH/F, P22320/P22317, da "
        "hidrogenase redutora de NAD. Não contar os quatro como quatro fontes de energia.",
        "3. **Duas CODHs não equivalentes:** Q2RUG7 em R. rubrum é anotada como "
        "família Ni-CODH, EC 1.2.7.4, com KO K00198 anaeróbico; Q0KEI3 em H16 é "
        "apenas 'aerobic-type carbon monoxide dehydrogenase homolog', predita, "
        "com domínio da família aldeído/xantina-oxidoredutase. Não demonstra CO como "
        "substrato fisiológico, centro Mo/Cu funcional nem carboxidotrofia de H16.",
        "4. **Erro de cepa concreto:** Q6QQP7, frequentemente recuperada por 'malonyl-CoA "
        "reductase', tem STRAIN=OK-70-fl nas referências, não J-10-fl. Os candidatos "
        "nativos J-10-fl são A9WIU3/Caur_2614 e A9WEU4/Caur_0613, com "
        "K14468/K14469 no KEGG, apesar de nomes UniProt genéricos/divergentes.",
        "5. **A evidência energética molecular mais direta desta seleção é localizada:** "
        "A9A1Y1 de SCM1 tem ECO:0000269, PMID 24843170, para a ligase "
        "4HB-CoA ADP-forming do ciclo HP/HB. Sua reação consome ATP e forma ADP+Pi. "
        "Isto informa o custo de uma etapa, não uma eficiência global em %, rendimento "
        "de biomassa ou superioridade geral do organismo.",
        "",
        "## Identidade das cepas e proteomas",
        "",
        "| Alvo | Identidade registrada | Taxon ID | Proteoma | Montagem vinculada |",
        "|---|---|---:|---|---|",
    ]
    for taxon in data["taxa"]:
        rec = taxon["record"]
        name = rec.get("taxonomy", {}).get("scientificName", "?")
        strain = rec.get("strain")
        if strain and strain not in name:
            name += f"; strain={strain}"
        assembly = rec.get("genomeAssembly", {}).get("assemblyId", "?")
        pid = taxon["proteome_id"]
        lines.append(
            f"| {taxon['target']} | {name} | {taxon['expected_taxon']} | "
            f"[{pid}](https://rest.uniprot.org/proteomes/{pid}.json) | {assembly} |"
        )
    lines += [
        "",
        "- DSM 2944 e GB17 são aliases no registro UP000273626; o taxon 82367 é "
        "de espécie, portanto o filtro foi **proteome:UP000273626**, não a espécie "
        "inteira. PD1222 é outra espécie/cepa, taxon 318586 ('Pd 1222'), e foi auditada "
        "separadamente. Não transferir os ensaios Sox entre elas.",
        "- TIE-1 é taxon **395960**. O proteoma está cadastrado com 5.246 proteínas, "
        "mas a consulta funcional em UniProtKB recuperou somente cinco entradas, "
        "todas revisadas. Pio e PufL foram recuperados via KEGG Genes, não declarados "
        "ausentes. Cadastro de proteoma e disponibilidade atual no UniProtKB não são sinônimos.",
        "- Entradas revisadas históricas de CGA009 e TIE-1 podem não ter ligação ao "
        "proteoma corrente. O JSON distingue ligação explícita de apenas identidade "
        "taxonômica exata. Tags RPA versus TX73 não foram alinhadas entre versões; "
        "entradas históricas e reanotações não foram contadas como duplicações gênicas.",
        "",
        "## Leitura da evidência",
        "",
        "**48 representantes**, mais **dois registros de cautela/exclusão**. "
        "Os níveis abaixo são de **anotação funcional**, não a escala de evidência "
        "de perda/pseudogenização: nenhuma perda biológica foi inferida.",
        "",
        "- **A1:** entrada não revisada; nome importado ou inferência automática.",
        "- **A2:** entrada revisada; sem ECO:0000269 nos campos FUNCTION/CATALYTIC ACTIVITY "
        "recuperados. Isso não nega experimentos existentes; limita esta extração.",
        "- **A3:** ECO:0000269 nesses campos, com fonte experimental identificável, "
        "limitada à proteína/reação correspondente.",
        "- **K1:** gene/KO do KEGG na cepa correta; o KEGG não fornece o estado "
        "Swiss-Prot/TrEMBL nem códigos ECO equivalentes.",
        "",
        "**R/U:** revisado/não revisado. **PE:** campo 'protein existence' da entrada; "
        "PE1 é evidência em nível proteico, não comprovação de todo o comentário funcional. "
        "Nas colunas ECO, N=nome/descrição, F=função/catálise; '—' significa código "
        "não fornecido nesse campo, não ausência de função ou inexistência de experimento.",
        "",
        "Códigos encontrados: ECO:0000313=informação importada em anotação automática; "
        "ECO:0000256=inferência por modelo/sequência em anotação automática; "
        "ECO:0000255=modelo/sequência em anotação manual; "
        "ECO:0000250=similaridade de sequência; ECO:0000305=inferência de curador; "
        "ECO:0000303=afirmação de autor; ECO:0000269=evidência experimental em anotação manual. "
        "O nome de Nmar_1500 traz **ECO:0008006, Google/ProtNLM**, não ensaio bioquímico.",
        "",
        "## Representantes confirmados nas bases",
    ]
    for target, _, _ in TARGETS:
        lines += [
            "",
            f"### {target}",
            "",
            "| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |",
            "|---|---|---|---|---|",
        ]
        for e in data["representatives"]:
            if e["target"] != target:
                continue
            if e["source_database"] == "UniProtKB":
                status = "R" if e["reviewed"] else "U"
                pe = (e.get("protein_existence") or "?").split(":")[0]
                n_codes = ", ".join(e["name_evidence_codes"]) or "—"
                f_codes = ", ".join(e["function_catalysis_evidence_codes"]) or "—"
                evidence = f"N: {n_codes}; F: {f_codes}"
            else:
                status, pe = "n/a", "n/a"
                evidence = "; ".join(e["kegg_record"].get("ORTHOLOGY", []))
            genes = ", ".join(gene_values(e)) or "não nomeado"
            lines.append(
                f"| [{e['accession']}]({e['source_url']}) | {e['audit_label']}; "
                f"{e['module']} | {genes} | {status}; PE{pe}; {e['evidence_level']} | "
                f"{evidence} |"
            )
        if target == "J-10-fl":
            lines += [
                "",
                "**Discordância entre bases:** UniProt chama A9WEU4 de "
                "'acetate--CoA ligase', EC 6.2.1.1, enquanto KEGG K14469 a atribui "
                "à Pcs trifuncional. Os códigos ECO da coluna F pertencem ao "
                "texto original UniProt, não validam automaticamente a atribuição "
                "alternativa KEGG. A9WIU3 recebe somente nome SDR no UniProt e "
                "Mcr/3HP-desidrogenase no KEGG K14468. Preservam-se as duas "
                "anotações e sua incerteza, sem inventar equivalência bioquímica.",
            ]
    lines += [
        "",
        "## Auditoria por módulo: afirmação permitida e extrapolação proibida",
        "",
        "### Luz: Puf/Puh e Bch",
        "Há representantes de centros de reação anoxigênicos e de biossíntese de "
        "bacterioclorofila nas bactérias fototróficas selecionadas. PufL/PuhA são "
        "componentes, e BchL/BchN também: nenhum isoladamente representa fotossistema "
        "completo, quantidade de pigmento, rendimento quântico ou crescimento à luz. "
        "Não exigir a mesma arquitetura Puh de bactérias púrpuras em Chloroflexus. "
        "A associação KEGG a um módulo/pathway é mapeamento de membro, não certificado "
        "de completude. Não foi inferida fotossíntese oxigênica.",
        "",
        "### H2: distinguir consumo, produção e regulação",
        "Os acessos de uptake selecionados são candidatos a consumo de H2. Não "
        "confundir hidrogenases com nitrogenase produtora de H2, nem fatores Hyp "
        "de maturação com enzimas de obtenção de energia. Em CGA009, a anotação "
        "estrutural Q6NB67/Q6NB68 não resolve integridade/regulação de hupUV, "
        "promotores, expressão ou crescimento usando H2; nenhuma lesão foi mapeada "
        "nesta auditoria. Em TIE-1, os genes KEGG Rpal_1154/1155 constituem evidência "
        "anotacional, não uma medição de atividade. H16 contém os dois sistemas "
        "distintos selecionados; não se mediu acoplamento, tolerância a O2 ou fluxo.",
        "",
        "**Cautela mesmo com Swiss-Prot:** P31891/P31892 possuem comentário genérico "
        "sobre reciclagem de H2 de nitrogenase, sem ECO nesse campo. Não usar esse "
        "comentário para afirmar diazotrofia da cepa H16. Revisado e PE1 não tornam "
        "toda frase uma medição específica da cepa.",
        "",
        "### Fixação de carbono",
        "CbbL forma I e CbbM forma II têm representantes, mas RuBisCO isolada não "
        "certifica o ciclo CBB completo, nem fixa seu fluxo ou prova autotrofia "
        "obrigatória. A auditoria de Chloroflexus aponta Mcr/Pcs do **3HP bi-cycle**, "
        "não exige CBB por analogia com bactérias púrpuras. O ciclo HP/HB de SCM1 "
        "também não deve ser substituído por um modelo CBB genérico. São módulos "
        "de assimilação de carbono, não fontes adicionais de energia por si mesmos.",
        "",
        "### Sox e aceptores",
        "Em DSM2944/GB17, SoxA/Y/B/C/D foram recuperados no proteoma correto; "
        "SoxB de PD1222 foi recuperada separadamente. São evidências candidatas "
        "para oxidação de compostos reduzidos de enxofre. Não equivalem a usar "
        "sulfato como fonte energética, e assimilação de sulfato não é Sox. "
        "Não se auditou completude de todos os pares SoxAX/SoxYZ ou sua maturação. "
        "NosZ A1B9T9 de PD1222 representa somente redução de N2O; não prova "
        "sozinho a sequência completa nitrato→N2, desnitrificação aeróbica ou taxas.",
        "",
        "### CODH",
        "R. rubrum Q2RUG7 tem inferência automática para família Ni-CODH, reação "
        "com ferredoxina, GO:0043885 com IEA, e KEGG K00198. H16 Q0KEI3 tem somente "
        "nome de homólogo tipo aeróbico importado do genoma; seus domínios são "
        "compartilhados com oxidoredutases de aldeído/xantina. A classificação "
        "bioquímica de CODH aeróbica Mo/Cu versus anaeróbica Ni/Fe-S não autoriza "
        "substituir uma pela outra, afirmar reversibilidade fisiológica, completar "
        "Wood–Ljungdahl, ou atribuir crescimento com CO a toda cepa com um homólogo. "
        "Não se validaram aqui motivos diagnósticos/tipo I versus tipo II, sítio "
        "catalítico, maturação ou especificidade de substrato de CoxL de H16.",
        "",
        "### Pio, ferro e eletrotrofia",
        "PioA é um citocromo multi-heme, PioB uma proteína de membrana externa "
        "associada e PioC uma HiPIP; **PioC não é uma oxidase de ferro isolada**. "
        "TIE-1 tem suporte de anotação KO e referências funcionais específicas "
        "nos KOs. CGA009 tem homólogos candidatos no proteoma atual. Não converter "
        "homologia Pio em fenótipo universal nem declarar Pio ausente da CGA009. "
        "Oxidação de Fe(II), transferência extracelular e entrada de elétrons de "
        "eletrodo não são três fontes independentes automaticamente demonstradas. "
        "Potencial do eletrodo, mediadores, luz, expressão e outros componentes "
        "precisariam de dados próprios; nada aqui mede eficiência coulômbica.",
        "",
        "### Amo em SCM1",
        "Foram retidas as três entradas genômicas Nmar_1500/1503/1502, com nomes "
        "Amo/AmoB/AmoC candidatos. Nmar_1500 tem nome inferido por ProtNLM e domínio "
        "InterPro AmoA_arc; o nome de Nmar_1503 também vem de ProtNLM. "
        "Isso não é evidência de atividade enzimática direta; "
        "não tratar similaridade à família amônia/metano monooxigenase como prova "
        "de metanotrofia. O produto energético é uma hipótese específica do sistema "
        "de oxidação de amônia, não de versatilidade extrema com muitos substratos. "
        "A9A1Y1 é uma etapa de fixação de carbono com anotação experimental separada.",
        "",
        "## Dois registros deliberadamente não usados para confirmar módulos",
        "",
        "| Acesso | Problema verificável | Consequência |",
        "|---|---|---|",
        "| [Q6QQP7](https://rest.uniprot.org/uniprotkb/Q6QQP7.json) | "
        "Mcr de OK-70-fl, taxon de espécie 1108; referências STRAIN=OK-70-fl | "
        "Não é evidência nativa da cepa J-10-fl; usar A9WIU3 com KO K14468 |",
        "| [Q0K5E4](https://rest.uniprot.org/uniprotkb/Q0K5E4.json) | "
        "H16: rótulo não revisado 'putative ammonia monooxygenase', gene abrB, "
        "sem demonstração do complexo AmoABC | Não promove H16 a nitrificador; "
        "não se conclui que a função rotulada seja biologicamente impossível |",
        "",
        "## O que estes dados não podem dizer",
        "",
        "A presença de uma entrada ou KO não demonstra ORF íntegra em todas as "
        "subculturas, via completa, disponibilidade de cofatores, expressão, fluxo, "
        "crescimento, rendimento, P/O, ATP por elétron, rendimento quântico, eficiência "
        "coulômbica ou porcentagem de conversão de energia. Não foi feita contagem "
        "de 'fontes de energia' a partir do número de genes. Não houve testes "
        "de ausência com controle positivo, genomas/leituras não foram reanalisados, "
        "e nenhum resultado zero foi tratado como ausência biológica. Não se "
        "inferiram perdas, ancestralidade, ramos, ordem de mutações ou trajetória "
        "evolutiva real. Esta é uma auditoria funcional de registros existentes, "
        "não uma reconstrução filogenética nem um desenho experimental.",
        "",
        "## Métodos e proveniência",
        "",
        "- Descoberta de ferramentas: ToolUniverse **find_tools**, consulta de "
        "UniProtRef_search_proteomes para R. palustris e documentação de "
        "UniProt_search; extração final pelo REST oficial público.",
        "- Endpoints: `https://rest.uniprot.org/proteomes/{UPID}.json`, "
        "`https://rest.uniprot.org/uniprotkb/search`; KEGG `/get/{IDs}` em lotes "
        "de até dez. URLs, parâmetros, horários, release, totais e status HTTP "
        "das requisições do script estão em `anotacoes.json:request_log`.",
        "- Consultas iniciais filtradas pela taxonomia **exata da cepa**; exceção "
        "DSM2944, que exigiu filtro por proteoma. Nomes/genéricos, tags e evidências "
        "foram verificados localmente. A busca exploratória de Mcr no taxon de "
        "espécie recuperou Q6QQP7 de outra cepa, explicitamente excluída.",
        "- Uma busca `gene:` com prefixos BDE18/Pden devolveu zero; buscas por "
        "acesso/tag textual e descrições recuperaram SoxB/C/D. Essa falha de "
        "recuperação não foi convertida em ausência dos genes.",
        "- O KEGG separou corretamente PioB/C (Rpal_0816/0815) de transportadores "
        "de sulfato adjacentes Rpal_0818/0819; não se usou contiguidade numérica "
        "como identidade funcional. Rpal_1715 é antena, não PufL; o representante "
        "validado é Rpal_1716.",
        "- ECO foi extraído por campo e fonte, não inventado a partir de 'reviewed'. "
        "Foram preservados os comentários de função/catálise/cofator, versões de "
        "entrada, referências, RefSeq/EMBL/KEGG/InterPro/GO e ligação ao proteoma. "
        "GO-IEA não foi tratado como ensaio direto.",
        "- Não foram lidos ou alterados `scripts/bioenergetics_efficiency_analysis.py` "
        "nem `runs/bioenergetics_efficiency_atlas.json`; nenhum número fixo desses "
        "arquivos foi reutilizado. Sem genomas grandes, NIM, NotebookLM ou previsões "
        "de eficiência.",
    ]
    count = data["counts"]
    releases = sorted(
        {r["uniprot_release"] for r in data["request_log"] if r.get("uniprot_release")}
    )
    lines += [
        "",
        f"**Reconciliação:** {count['initial_targeted_uniprot_hits']} hits nas oito "
        f"buscas iniciais; {count['candidate_records_before_selection']} acessos "
        "candidatos únicos após buscas suplementares; "
        f"{count['module_representatives_uniprot']} representantes UniProt e "
        f"{count['module_representatives_kegg']} KEGG retidos; "
        f"{count['excluded_cautions']} cautelas separadas. "
        f"{count['script_http_requests']} requisições HTTP no script. "
        f"Release UniProt registrada: {', '.join(releases) or 'header indisponível'}. "
        "Todos os conjuntos de busca do script couberam na primeira página, sem "
        "Link next. Isto demonstra completude das consultas registradas, "
        "**não** completude de cada via no genoma.",
        "",
        "### Fontes-chave, com o escopo correto",
        "",
        "- [KO K28459](https://rest.kegg.jp/get/ko:K28459), "
        "[K28460](https://rest.kegg.jp/get/ko:K28460), "
        "[K28461](https://rest.kegg.jp/get/ko:K28461): identidades Pio e "
        "referências de TIE-1; listam também genes CGA009.",
        "- [A9A1Y1](https://rest.uniprot.org/uniprotkb/A9A1Y1.json), "
        "[PMID 24843170](https://pubmed.ncbi.nlm.nih.gov/24843170/): "
        "anotação experimental da ligase ADP-forming de SCM1.",
        "- [Caur_2614](https://rest.kegg.jp/get/cau:Caur_2614) e "
        "[Caur_0613](https://rest.kegg.jp/get/cau:Caur_0613): "
        "Mcr/Pcs candidatos na linhagem genômica J-10-fl.",
        "- [Q2RUG7](https://rest.uniprot.org/uniprotkb/Q2RUG7.json) e "
        "[Rru_A1427](https://rest.kegg.jp/get/rru:Rru_A1427): Ni-CODH anaeróbica "
        "candidata de ATCC11170; [Q0KEI3](https://rest.uniprot.org/uniprotkb/Q0KEI3.json): "
        "homólogo tipo aeróbico de H16, com incerteza de substrato mantida.",
        "",
    ]
    (OUT / "anotacoes.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        choices=(
            "collect",
            "inspect",
            "supplement",
            "evidence",
            "resolve-gaps",
            "pio-compare",
            "finalize",
        ),
    )
    parser.add_argument("accessions", nargs="*")
    args = parser.parse_args()
    if args.mode == "evidence":
        show_evidence(args.accessions)
    else:
        {
            "collect": collect,
            "inspect": inspect,
            "supplement": supplement,
            "resolve-gaps": resolve_gaps,
            "pio-compare": pio_compare,
            "finalize": finalize,
        }[args.mode]()
