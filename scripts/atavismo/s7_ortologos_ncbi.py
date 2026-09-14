"""Distribuicao taxonomica de ortologos (NCBI Datasets v2) para testar se ODAM/KLK4
sao inovacoes de mamifero."""
from __future__ import annotations
import json, time
from pathlib import Path
import requests
API="https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT=Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S=requests.Session(); S.headers.update({"Accept":"application/json"})

GENES={"ENAM":10117,"AMBN":258,"AMELX":265,"AMTN":401138,"ODAM":54959,
       "KLK4":9622,"MMP20":9313,"DSPP":1834,"ACP4":93650,"FAM83H":286077,
       "HES7":84667,"SLC24A4":123041,"WDR72":256764}
# grupos taxonomicos de interesse (ordem = do mais distante ao mamifero)
CLADOS={"Mammalia":40674,"Aves":8782,"Crocodylia":1294634,"Testudines":8459,
        "Lepidosauria":8504,"Amphibia":8292,"Actinopterygii":7898}

res={}
for sym,gid in GENES.items():
    try:
        r=S.get(f"{API}/gene/id/{gid}/orthologs",params={"returned_content":"COMPLETE","page_size":1000},timeout=120)
        r.raise_for_status(); d=r.json()
    except Exception as exc:
        res[sym]={"erro":str(exc)}; print(sym,"ERRO",exc); continue
    genes=d.get("reports",[])
    spp=[]
    for g in genes:
        gg=g.get("gene",{})
        spp.append({"taxname":gg.get("taxname"),"taxid":gg.get("tax_id"),
                    "symbol":gg.get("symbol"),"gene_id":gg.get("gene_id"),
                    "type":gg.get("type")})
    res[sym]={"n_ortologos":len(spp),"especies":spp}
    names=sorted({s["taxname"] for s in spp})
    print(f"\n### {sym} (HS gene {gid}): {len(spp)} ortologos em NCBI")
    print("   ", "; ".join(names[:60]))
    time.sleep(0.4)

OUT.joinpath("s7_ortologos_ncbi.json").write_text(json.dumps(res,indent=2),encoding="utf-8")
