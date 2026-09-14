"""Baixa relatorios de gene (coordenadas + acesso proteico) do NCBI Datasets."""
from __future__ import annotations
import json, time
from pathlib import Path
import requests
API="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT=Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S=requests.Session(); S.headers.update({"Accept":"application/json"})

ALVOS=[
 ("Alligator mississippiensis",8496,["ENAM","AMBN","AMTN","AMELX","MMP20","DSPP","ACP4"]),
 ("Anolis carolinensis",28377,["ENAM","AMBN","AMTN","AMELX","MMP20"]),
 ("Chrysemys picta bellii",8478,["ENAM","AMBN"]),
 ("Chlamydotis macqueenii",187382,["ENAM"]),
 ("Struthio camelus",8801,["AMELX"]),
 ("Aptenodytes forsteri",9233,["AMELX","MMP20"]),
 ("Corvus moneduloides",1196302,["AMBN"]),
 ("Strigops habroptila",2489341,["HES7"]),
 ("Falco cherrug",345164,["HES7"]),
 ("Gallus gallus",9031,["SULT1E1","LPL","JCHAIN","FAM83H","SDAD1","FAT4"]),
]
rep={}
for sp,tax,syms in ALVOS:
    for sym in syms:
        key=f"{sp}|{sym}"
        try:
            r=S.get(f"{API}/gene/symbol/{sym}/taxon/{tax}",params={"returned_content":"COMPLETE"},timeout=90)
            r.raise_for_status(); d=r.json()
        except Exception as e:
            rep[key]={"erro":str(e)}; print(key,"ERRO",e); continue
        recs=d.get("reports",[])
        keep=[]
        for rc in recs:
            g=rc.get("gene",{})
            if g.get("symbol","").upper()!=sym.upper(): continue
            ranges=g.get("genomic_ranges",[])
            loc=None
            if ranges:
                rr=ranges[0]; rgs=rr.get("range",[{}])[0]
                loc={"acc":rr.get("accession_version"),
                     "begin":int(rgs.get("begin",0)),"end":int(rgs.get("end",0)),
                     "orient":rgs.get("orientation")}
            prots=[]
            for t in g.get("transcripts",[]) or []:
                p=t.get("protein")
                if p and p.get("accession_version"):
                    prots.append({"acc":p["accession_version"],"len":p.get("length"),
                                  "tx":t.get("accession_version"),"name":t.get("name")})
            keep.append({"gene_id":g.get("gene_id"),"taxname":g.get("taxname"),
                         "type":g.get("type"),"desc":g.get("description"),
                         "loc":loc,"prots":prots,
                         "annot":(g.get("annotations") or [{}])[0].get("assembly_name")})
        rep[key]=keep
        for k in keep:
            pl=max([p["len"] or 0 for p in k["prots"]],default=0)
            print(f"{key:42} id={k['gene_id']:<12} {k['type']:<14} {k['loc']['acc'] if k['loc'] else '?'}:"
                  f"{k['loc']['begin'] if k['loc'] else 0}-{k['loc']['end'] if k['loc'] else 0} "
                  f"{k['loc']['orient'] if k['loc'] else ''} prot_max={pl}aa n_prot={len(k['prots'])}")
        if not keep: print(f"{key:42} NENHUM")
        time.sleep(0.3)
(OUT/"s9_genes.json").write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding="utf-8")
