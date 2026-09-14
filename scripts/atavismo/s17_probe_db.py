"""Sonda quais bancos do NCBI BLAST contem sequencia genomica de Gallus gallus."""
from __future__ import annotations
import re, time
from pathlib import Path
import requests
URL="https://blast.ncbi.nlm.nih.gov/Blast.cgi"; SEQ=Path(r"C:\Users\Manec\bio-stack\runs\atavismo\seqs")
S=requests.Session()
q=[">"+b for b in (SEQ/"ref_proteinas.fa").read_text(encoding="utf-8").split(">") if b.startswith("Ggal_SULT1E1")][0]
CAND=["GCF_016699485.2","GPIPE/9031/current/all_top_level","GPIPE/9031/current/ref_top_level",
      "refseq_genomes","wgs","refseq_select_rna","core_nt"]
rids={}
for db in CAND:
    p={"CMD":"Put","PROGRAM":"tblastn","DATABASE":db,"QUERY":q,"EXPECT":"1e-5","HITLIST_SIZE":"10",
       "FILTER":"F","EMAIL":"limatheus458@gmail.com","TOOL":"bio-stack"}
    if db in ("refseq_genomes","wgs","core_nt"): p["ENTREZ_QUERY"]="txid9031[ORGN]"
    try:
        r=S.post(URL,data=p,timeout=180)
        m=re.search(r"RID = (\S+)",r.text)
        if m: rids[db]=m.group(1); print(f"{db:40} submetido RID={m.group(1)}",flush=True)
        else: print(f"{db:40} REJEITADO",flush=True)
    except Exception as e: print(f"{db:40} ERRO {e}",flush=True)
    time.sleep(5)
for db,rid in rids.items():
    for _ in range(30):
        time.sleep(15)
        s=S.get(URL,params={"CMD":"Get","FORMAT_OBJECT":"SearchInfo","RID":rid},timeout=120).text
        st=re.search(r"Status=(\w+)",s)
        if st and st.group(1) in ("READY","UNKNOWN"): break
    t=S.get(URL,params={"CMD":"Get","RID":rid,"FORMAT_TYPE":"Tabular","ALIGNMENTS":"10","DESCRIPTIONS":"10"},timeout=180).text
    rows=[l for l in t.splitlines() if l.count("\t")>=10 and not l.startswith("#")]
    print(f"\n### {db:40} {len(rows)} HSPs")
    for l in rows[:3]:
        f=l.split("\t"); print(f"    {f[1][:40]:40} pid={f[2]} E={f[10]}")
