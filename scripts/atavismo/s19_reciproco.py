"""Teste reciproco: o fragmento de 42 aa do chrZ de Gallus tem AMBN como melhor acerto?"""
from __future__ import annotations
import re, time
import requests
URL="https://blast.ncbi.nlm.nih.gov/Blast.cgi"; S=requests.Session()
FRAG=("LLLYHGYYEYGEPFNIIWCYGFLSPYSSYPWPRQRAQKYETQ")
ENAM=(">psi_ENAM_frag_Gallus_chrZ\n")
def run(seq,name,db="nr",entrez=None,ws="2"):
    p={"CMD":"Put","PROGRAM":"blastp","DATABASE":db,"QUERY":f">{name}\n{seq}\n","EXPECT":"200",
       "HITLIST_SIZE":"30","FILTER":"F","WORD_SIZE":ws,"MATRIX_NAME":"BLOSUM62",
       "EMAIL":"limatheus458@gmail.com","TOOL":"bio-stack"}
    if entrez: p["ENTREZ_QUERY"]=entrez
    r=S.post(URL,data=p,timeout=180)
    m=re.search(r"RID = (\S+)",r.text)
    if not m: print(name,"submissao falhou"); return
    rid=m.group(1); print(f"{name}: RID={rid}",flush=True)
    for _ in range(50):
        time.sleep(20)
        s=S.get(URL,params={"CMD":"Get","FORMAT_OBJECT":"SearchInfo","RID":rid},timeout=120).text
        st=re.search(r"Status=(\w+)",s)
        if st and st.group(1) in ("READY","UNKNOWN"): break
    t=S.get(URL,params={"CMD":"Get","RID":rid,"FORMAT_TYPE":"Tabular","ALIGNMENTS":"30","DESCRIPTIONS":"30"},timeout=180).text
    rows=[l.split("\t") for l in t.splitlines() if l.count("\t")>=10 and not l.startswith("#")]
    print(f"\n### {name}: {len(rows)} HSPs (blastp vs {db})")
    # resolver titulos
    accs=[r_[1] for r_ in rows[:15]]
    tit={}
    if accs:
        e=S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db":"protein","id":",".join(accs),"retmode":"json",
                        "tool":"bio-stack","email":"limatheus458@gmail.com"},timeout=120).json().get("result",{})
        for k,v in e.items():
            if isinstance(v,dict): tit[v.get("accessionversion","")]=f"{v.get('title','')[:62]}"
    for r_ in rows[:15]:
        print(f"   {r_[1][:18]:18} pid={r_[2]:>6} len={r_[3]:>4} E={r_[10]:>10} bits={r_[11]:>6}  {tit.get(r_[1],'')}")
run(FRAG,"Ggal_chrZ_55168855_55168983_f2")
