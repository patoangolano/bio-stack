"""Diagnostico: o cliente BLAST do NCBI funciona? Controle que OBRIGATORIAMENTE tem acerto."""
from __future__ import annotations
import re, time
from pathlib import Path
import requests
URL="https://blast.ncbi.nlm.nih.gov/Blast.cgi"; SEQ=Path(r"C:\Users\Manec\bio-stack\runs\atavismo\seqs")
S=requests.Session()
def getfa(lab):
    for b in (SEQ/"ref_proteinas.fa").read_text(encoding="utf-8").split(">"):
        if b.startswith(lab): return ">"+b
    raise KeyError(lab)
FRAG=">frag42\nLLLYHGYYEYGEPFNIIWCYGFLSPYSSYPWPRQRAQKYETQ\n"
JOBS=[("CTRL_SULT1E1_nr",getfa("Ggal_SULT1E1"),"nr","3","10"),
      ("CTRL_AMBN_Amis_nr",getfa("Amis_AMBN"),"nr","3","10"),
      ("FRAG42_nr",FRAG,"nr","2","200")]
for nome,q,db,ws,ev in JOBS:
    p={"CMD":"Put","PROGRAM":"blastp","DATABASE":db,"QUERY":q,"EXPECT":ev,"HITLIST_SIZE":"20",
       "FILTER":"F","WORD_SIZE":ws,"EMAIL":"limatheus458@gmail.com","TOOL":"bio-stack"}
    r=S.post(URL,data=p,timeout=200)
    m=re.search(r"RID = (\S+)",r.text)
    if not m:
        print(f"{nome:22} SEM RID | erros: {re.findall(r'Error:.*',r.text)[:2]}"); continue
    rid=m.group(1); print(f"{nome:22} RID={rid}",flush=True)
    status="?"
    for _ in range(55):
        time.sleep(20)
        s=S.get(URL,params={"CMD":"Get","FORMAT_OBJECT":"SearchInfo","RID":rid},timeout=150).text
        mm=re.search(r"Status=(\w+)",s)
        status=mm.group(1) if mm else "?"
        if status in ("READY","UNKNOWN"): break
    hits=re.search(r"ThereAreHits=(\w+)",s)
    t=S.get(URL,params={"CMD":"Get","RID":rid,"FORMAT_TYPE":"Tabular","ALIGNMENTS":"20","DESCRIPTIONS":"20"},timeout=200).text
    rows=[l for l in t.splitlines() if l.count("\t")>=10 and not l.startswith("#")]
    Path(rf"C:\Users\Manec\bio-stack\runs\atavismo\diag_{nome}.txt").write_text(t,encoding="utf-8")
    print(f"{'':22} status={status} ThereAreHits={hits.group(1) if hits else '?'} "
          f"linhas_tabulares={len(rows)} bytes={len(t)}")
    print(f"{'':22} amostra bruta: {t[:150].replace(chr(10),' | ')}")
    for l in rows[:4]:
        f=l.split("\t"); print(f"{'':24}{f[1][:26]:26} pid={f[2]:>6} E={f[10]:>10}")
