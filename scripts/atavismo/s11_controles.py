from __future__ import annotations

import json
import time
from pathlib import Path

import requests

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
SEQ = OUT / "seqs"
CM = {"tool": "bio-stack", "email": "limatheus458@gmail.com"}
S = requests.Session()
# controles positivos internos: ortologos de Alligator de genes que existem no frango
# e que ficam DENTRO da mesma janela de consulta
ALVOS = {
    "Amis_SULT1E1": ("SULT1E1", 8496),
    "Amis_FAT4": ("FAT4", 8496),
    "Amis_JCHAIN": ("JCHAIN", 8496),
}
meta = json.loads((OUT / "s10_meta.json").read_text(encoding="utf-8"))
for lab, (sym, tax) in ALVOS.items():
    r = S.get(
        f"{API}/gene/symbol/{sym}/taxon/{tax}", params={"returned_content": "COMPLETE"}, timeout=90
    )
    reps = [
        x["gene"] for x in r.json().get("reports", []) if x["gene"]["symbol"].upper() == sym.upper()
    ]
    if not reps:
        print(lab, "NENHUM")
        continue
    g = reps[0]
    gid = g["gene_id"]
    a = g["annotations"][0]
    gl = a["genomic_locations"][0]
    gr = gl["genomic_range"]
    el = S.get(
        f"{EU}/elink.fcgi",
        params={**CM, "dbfrom": "gene", "db": "protein", "id": gid, "retmode": "json"},
        timeout=90,
    ).json()
    ids = sorted(
        {
            i
            for ls in el["linksets"]
            for db in ls.get("linksetdbs", [])
            if "protein" in db["linkname"]
            for i in db["links"]
        }
    )
    fa = S.get(
        f"{EU}/efetch.fcgi",
        params={
            **CM,
            "db": "protein",
            "id": ",".join(ids[:3]),
            "rettype": "fasta",
            "retmode": "text",
        },
        timeout=180,
    ).text
    (SEQ / f"prot_{lab}.fa").write_text(fa, encoding="utf-8")
    lens = [len("".join(b.splitlines()[1:])) for b in fa.split(">") if b.strip()]
    meta[lab] = {
        "gene_id": gid,
        "symbol": sym,
        "taxname": g["taxname"],
        "acc": gl["genomic_accession_version"],
        "begin": int(gr["begin"]),
        "end": int(gr["end"]),
        "orient": gr["orientation"],
        "prot_lens": lens,
    }
    print(
        f"{lab:14} id={gid} {gl['genomic_accession_version']}:{gr['begin']}-{gr['end']} prots={lens}"
    )
    time.sleep(0.4)
(OUT / "s10_meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
