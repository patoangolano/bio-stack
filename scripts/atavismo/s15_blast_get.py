"""Recupera os resultados tabulares dos RIDs do NCBI BLAST e resume."""
from __future__ import annotations
import re, sys, time
from pathlib import Path
import requests
URL="https://blast.ncbi.nlm.nih.gov/Blast.cgi"; OUT=Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S=requests.Session()
RIDS={"CTRL_FAT4_Amis_vs_Gallus":"AGZATVHG014","TESTE_ENAM_Amis_vs_Gallus":"AGZB13JW014",
      "TESTE_ENAM_Cmac_vs_Gallus":"AGZBA0A9014","TESTE_AMBN_Amis_vs_Gallus":"AGZBJB6J014",
      "TESTE_AMELX_Scam_vs_Gallus":"AGZBUUUD014"}
for nome,rid in RIDS.items():
    for a in range(40):
        s=S.get(URL,params={"CMD":"Get","FORMAT_OBJECT":"SearchInfo","RID":rid},timeout=120).text
        st=(re.search(r"Status=(\w+)",s) or [None,"?"])[1] if re.search(r"Status=(\w+)",s) else "?"
        if st=="READY": break
        time.sleep(20)
    r=S.get(URL,params={"CMD":"Get","RID":rid,"FORMAT_TYPE":"Tabular","ALIGNMENTS":"100",
                        "DESCRIPTIONS":"100"},timeout=180)
    txt=r.text
    (OUT/f"blast_{nome}.txt").write_text(txt,encoding="utf-8")
    # extrair bloco tabular real
    rows=[l for l in txt.splitlines() if l.count("\t")>=10 and not l.startswith("#")]
    print(f"\n===== {nome}  RID={rid}  status={st}  {len(rows)} HSPs (tblastn, core_nt, txid9031, E<=1e-3) =====")
    if not rows: 
        print("   SEM HITS")
        continue
    for l in rows[:10]:
        f=l.split("\t")
        print(f"   {f[1][:34]:34} pid={f[2]:>6} len={f[3]:>5} s={f[8]}-{f[9]} E={f[10]:>10} bits={f[11]}")
