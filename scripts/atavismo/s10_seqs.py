"""Baixa proteinas (elink gene->protein + efetch) e regioes genomicas (efetch nuccore)."""

from __future__ import annotations

import json
import time
from pathlib import Path

import requests

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
OUT.mkdir(parents=True, exist_ok=True)
SEQ = OUT / "seqs"
SEQ.mkdir(exist_ok=True)
CM = {"tool": "bio-stack", "email": "limatheus458@gmail.com"}
S = requests.Session()

GENES = {  # rotulo -> gene_id NCBI
    "Amis_ENAM": 102563791,
    "Amis_AMBN": 102564019,
    "Amis_AMTN": 102564246,
    "Amis_AMELX": 102577340,
    "Amis_MMP20": 102563202,
    "Amis_DSPP": 102566454,
    "Amis_ACP4": 102558661,
    "Acar_ENAM": 100499197,
    "Acar_AMBN": 100554339,
    "Acar_AMTN": 100554538,
    "Acar_AMELX": 100552162,
    "Acar_MMP20": 100553435,
    "Cpic_ENAM": 112060187,
    "Cmac_ENAM": 147596333,
    "Scam_AMELX": 104149831,
    "Afor_AMELX": 103906029,
    "Afor_MMP20": 103895117,
    "Cmon_AMBN": 116455309,
    "Shab_HES7": 115603040,
    "Fche_HES7": 114014481,
    "Ggal_LPL": 396219,
    "Ggal_SULT1E1": 422662,
}


def loc(gid: int):
    r = S.get(f"{API}/gene/id/{gid}", params={"returned_content": "COMPLETE"}, timeout=90)
    r.raise_for_status()
    g = r.json()["reports"][0]["gene"]
    a = g["annotations"][0]
    gl = a["genomic_locations"][0]
    gr = gl["genomic_range"]
    return {
        "gene_id": gid,
        "symbol": g["symbol"],
        "taxname": g["taxname"],
        "asm": a["assembly_accession"],
        "asm_name": a["assembly_name"],
        "acc": gl["genomic_accession_version"],
        "begin": int(gr["begin"]),
        "end": int(gr["end"]),
        "orient": gr["orientation"],
        "chrom": gl.get("sequence_name"),
        "n_prot": g.get("protein_count"),
        "desc": g.get("description"),
    }


def prots(gid: int):
    r = S.get(
        f"{EU}/elink.fcgi",
        params={**CM, "dbfrom": "gene", "db": "protein", "id": gid, "retmode": "json"},
        timeout=90,
    )
    r.raise_for_status()
    d = r.json()
    ids = []
    for ls in d["linksets"]:
        for db in ls.get("linksetdbs", []):
            if db["linkname"] in ("gene_protein", "gene_protein_refseq"):
                ids += db["links"]
    return sorted(set(ids))


def efetch_prot(ids):
    r = S.get(
        f"{EU}/efetch.fcgi",
        params={**CM, "db": "protein", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"},
        timeout=180,
    )
    r.raise_for_status()
    return r.text


meta = {}
for lab, gid in GENES.items():
    try:
        L = loc(gid)
        meta[lab] = L
        ids = prots(gid)
        fa = efetch_prot(ids[:6]) if ids else ""
        (SEQ / f"prot_{lab}.fa").write_text(fa, encoding="utf-8")
        lens = [len("".join(b.splitlines()[1:])) for b in fa.split(">") if b.strip()]
        meta[lab]["prot_ids"] = ids[:6]
        meta[lab]["prot_lens"] = lens
        span = L["end"] - L["begin"] + 1
        print(
            f"{lab:14} {L['taxname'][:26]:26} {L['asm_name']:<22} {L['acc']}:{L['begin']}-{L['end']} ({span / 1000:.1f} kb) {L['orient'][:5]} prots={lens}"
        )
    except Exception as e:
        meta[lab] = {"erro": str(e)}
        print(lab, "ERRO", e)
    time.sleep(0.4)
(OUT / "s10_meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
