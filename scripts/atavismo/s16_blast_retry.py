from __future__ import annotations

import re
import time
from pathlib import Path

import requests

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
SEQ = Path(r"C:\Users\Manec\bio-stack\runs\atavismo\seqs")
S = requests.Session()


def fa(lab):
    for b in (SEQ / "ref_proteinas.fa").read_text(encoding="utf-8").split(">"):
        if b.startswith(lab):
            return ">" + b
    raise KeyError(lab)


def run(lab, db, entrez, tag):
    p = {
        "CMD": "Put",
        "PROGRAM": "tblastn",
        "DATABASE": db,
        "QUERY": fa(lab),
        "EXPECT": "1e-3",
        "HITLIST_SIZE": "50",
        "FILTER": "F",
        "EMAIL": "limatheus458@gmail.com",
        "TOOL": "bio-stack",
    }
    if entrez:
        p["ENTREZ_QUERY"] = entrez
    r = S.post(URL, data=p, timeout=180)
    m = re.search(r"RID = (\S+)", r.text)
    if not m:
        print(f"{tag:34} SUBMISSAO FALHOU: {re.findall(r'Error:.*', r.text)[:2] or r.text[:200]}")
        return None
    rid = m.group(1)
    for _ in range(45):
        time.sleep(20)
        s = S.get(
            URL, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}, timeout=120
        ).text
        st = re.search(r"Status=(\w+)", s)
        if st and st.group(1) == "READY":
            break
        if st and st.group(1) == "UNKNOWN":
            print(f"{tag:34} RID expirou")
            return None
    t = S.get(
        URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "Tabular",
            "ALIGNMENTS": "100",
            "DESCRIPTIONS": "100",
        },
        timeout=180,
    ).text
    rows = [
        linha for linha in t.splitlines() if linha.count("\t") >= 10 and not linha.startswith("#")
    ]
    print(f"{tag:34} db={db:32} RID={rid} -> {len(rows)} HSPs")
    for linha in rows[:6]:
        f = linha.split("\t")
        print(
            f"      {f[1][:30]:30} pid={f[2]:>6} len={f[3]:>5} s={f[8]}-{f[9]} E={f[10]:>10} bits={f[11]}"
        )
    Path(rf"C:\Users\Manec\bio-stack\runs\atavismo\blast2_{tag}.txt").write_text(
        t, encoding="utf-8"
    )
    return len(rows)


for db in ["refseq_representative_genomes", "refseq_genomes", "wgs"]:
    n = run("Amis_FAT4", db, "txid9031[ORGN]", f"CTRL_FAT4_{db}")
    if n:
        print(f"--> banco {db} funciona; rodando os testes")
        for lab in ["Amis_ENAM", "Cmac_ENAM", "Amis_AMBN", "Scam_AMELX", "Shab_HES7"]:
            run(lab, db, "txid9031[ORGN]", f"TESTE_{lab}_{db}")
        break
