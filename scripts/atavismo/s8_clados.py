from __future__ import annotations
import json, time
from pathlib import Path
import requests
API="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT=Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S=requests.Session()
res=json.loads((OUT/"s7_ortologos_ncbi.json").read_text(encoding="utf-8"))

taxids=sorted({int(s["taxid"]) for g in res.values() for s in g.get("especies",[]) if s.get("taxid")})
print("taxids unicos:",len(taxids))

CLADOS=[("Mammalia",40674),("Aves",8782),("Crocodylia",1294634),("Testudines",8459),
        ("Lepidosauria",8504),("Amphibia",8292),("Actinopterygii",7898),
        ("Chondrichthyes",7777),("Coelacanthimorpha",118072),("Dipnomorpha",7878)]
lineage={}
CH=200
for i in range(0,len(taxids),CH):
    chunk=[str(t) for t in taxids[i:i+CH]]
    for a in range(5):
        try:
            r=S.post(f"{API}/taxonomy",json={"taxons":chunk,"returned_content":"COMPLETE"},timeout=180)
            r.raise_for_status(); d=r.json(); break
        except Exception as exc:
            print("retry",exc); time.sleep(3*(a+1)); d=None
    if not d: continue
    for rep in d.get("taxonomy_nodes",[]):
        t=rep.get("taxonomy",{})
        lineage[t.get("tax_id")]=set(t.get("lineage",[])) | {t.get("tax_id")}
    time.sleep(0.3)
print("linhagens obtidas:",len(lineage))

def clado(tid:int)->str:
    lin=lineage.get(tid,set())
    for nome,ct in CLADOS:
        if ct in lin: return nome
    return "outro/desconhecido"

tabela={}
print(f"\n{'gene':10} " + " ".join(f"{c[:6]:>7}" for c,_ in CLADOS))
for sym,g in res.items():
    cont={c:0 for c,_ in CLADOS}; cont["outro/desconhecido"]=0
    aves=[]
    for s in g.get("especies",[]):
        c=clado(int(s["taxid"])); cont[c]=cont.get(c,0)+1
        if c=="Aves": aves.append(f"{s['taxname']} [{s['symbol']}|{s['gene_id']}|{s.get('type')}]")
    tabela[sym]={"contagem":cont,"aves":sorted(aves)}
    print(f"{sym:10} " + " ".join(f"{cont[c]:>7}" for c,_ in CLADOS) + f"   outro={cont['outro/desconhecido']}")

print("\n=== AVES com ortologo anotado (por gene) ===")
for sym,v in tabela.items():
    if v["aves"]:
        print(f"\n-- {sym}: {len(v['aves'])} aves")
        for a in v["aves"][:40]: print("   ",a)

(OUT/"s8_clados.json").write_text(json.dumps(tabela,indent=2,ensure_ascii=False),encoding="utf-8")
