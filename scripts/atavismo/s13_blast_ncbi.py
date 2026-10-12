"""tBLASTn genoma-inteiro de Gallus gallus via NCBI BLAST URL API, com controle positivo."""

from __future__ import annotations

import re
import time
from pathlib import Path

import requests

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
SEQ = Path(r"C:\Users\Manec\bio-stack\runs\atavismo\seqs")
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S = requests.Session()


def fa(lab):
    txt = (SEQ / "ref_proteinas.fa").read_text(encoding="utf-8")
    for b in txt.split(">"):
        if b.startswith(lab):
            return ">" + b
    raise KeyError(lab)


def submit(seq, db, entrez):
    p = {
        "CMD": "Put",
        "PROGRAM": "tblastn",
        "DATABASE": db,
        "QUERY": seq,
        "ENTREZ_QUERY": entrez,
        "EXPECT": "1e-3",
        "HITLIST_SIZE": "50",
        "FILTER": "F",
        "WORD_SIZE": "3",
        "EMAIL": "limatheus458@gmail.com",
        "TOOL": "bio-stack",
    }
    r = S.post(URL, data=p, timeout=180)
    r.raise_for_status()
    m = re.search(r"RID = (\S+)", r.text)
    rt = re.search(r"RTOE = (\d+)", r.text)
    if not m:
        err = re.findall(r"(?:Message ID#\d+\s*)?Error:.*", r.text)[:3]
        raise RuntimeError(f"sem RID. {err or r.text[:300]}")
    return m.group(1), int(rt.group(1)) if rt else 30


def wait(rid, maxmin=25):
    t0 = time.time()
    while time.time() - t0 < maxmin * 60:
        time.sleep(20)
        r = S.get(
            URL, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}, timeout=120
        )
        st = re.search(r"Status=(\w+)", r.text)
        st = st.group(1) if st else "?"
        if st == "READY":
            hits = re.search(r"ThereAreHits=(\w+)", r.text)
            return True, (hits.group(1) if hits else "?")
        if st == "UNKNOWN":
            return False, "UNKNOWN/expirado"
        print(f"  [{rid}] {st} ({int(time.time() - t0)}s)", flush=True)
    return False, "timeout"


def tab(rid):
    r = S.get(
        URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "Tabular",
            "ALIGNMENTS": "50",
            "DESCRIPTIONS": "50",
            "NCBI_GI": "F",
        },
        timeout=180,
    )
    return r.text


JOBS = [
    ("CTRL_FAT4_Amis_vs_Gallus", "Amis_FAT4"),
    ("TESTE_ENAM_Amis_vs_Gallus", "Amis_ENAM"),
    ("TESTE_ENAM_Cmac_vs_Gallus", "Cmac_ENAM"),
    ("TESTE_AMBN_Amis_vs_Gallus", "Amis_AMBN"),
    ("TESTE_AMELX_Scam_vs_Gallus", "Scam_AMELX"),
]
DB = "core_nt"
ENTREZ = "txid9031[ORGN]"
rids = {}
for nome, lab in JOBS:
    try:
        rid, rtoe = submit(fa(lab), DB, ENTREZ)
        rids[nome] = rid
        print(f"submetido {nome:30} query={lab:12} RID={rid} (RTOE~{rtoe}s)", flush=True)
    except Exception as e:
        print(f"FALHOU {nome}: {e}", flush=True)
    time.sleep(6)
for nome, rid in rids.items():
    ok, info = wait(rid)
    print(f"\n===== {nome}  RID={rid}  ready={ok} hits={info} =====", flush=True)
    if ok:
        t = tab(rid)
        (OUT / f"blast_{nome}.tsv").write_text(t, encoding="utf-8")
        lines = [linha for linha in t.splitlines() if linha and not linha.startswith("#")]
        print(f"  {len(lines)} linhas de hit; primeiras 8:")
        for linha in lines[:8]:
            print("   ", linha)
