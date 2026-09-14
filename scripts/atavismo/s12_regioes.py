"""Baixa as regioes genomicas (consultas) e monta o banco de proteinas de referencia."""
from __future__ import annotations
import json, time
from pathlib import Path
import requests
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT=Path(r"C:\Users\Manec\bio-stack\runs\atavismo"); SEQ=OUT/"seqs"; SEQ.mkdir(exist_ok=True)
CM={"tool":"bio-stack","email":"limatheus458@gmail.com"}; S=requests.Session()

REGIOES={
 # rotulo: (accession, inicio, fim, descricao)
 "GG_void":      ("NC_052535.1",51_250_000,52_200_000,"Gallus chr4 SULT1E1->FAT4 (locus sintenico do cluster SCPP)"),
 "GG_juncao":    ("NC_052535.1",49_450_000,49_560_000,"Gallus chr4 SDAD1-UTP3-JCHAIN-RUFY3"),
 "AMIS_ctrl":    ("NC_081826.1",157_270_000,157_400_000,"Alligator ENAM-AMBN-AMTN (controle +, mesmo taxon)"),
 "ACAR_ctrl":    ("NC_085842.1",199_180_000,199_320_000,"Anolis ENAM-AMBN-AMTN (controle +, ~280 Ma)"),
 "CPIC_ctrl":    ("NC_088796.1",129_090_000,129_160_000,"Chrysemys ENAM (controle +, ~255 Ma, edentulo)"),
 "CMAC_ctrl":    ("NC_141593.1",50_555_000,50_625_000,"Chlamydotis ENAM (controle +, AVE com reliquia anotada)"),
 "SCAM_ctrl":    ("NC_090942.1",141_825_000,141_865_000,"Struthio AMELX (controle +, AVE paleognata)"),
}
for lab,(acc,s,e,desc) in REGIOES.items():
    f=SEQ/f"reg_{lab}.fa"
    if f.exists() and f.stat().st_size>1000:
        print(f"{lab:12} ja existe ({f.stat().st_size/1e6:.2f} MB)"); continue
    r=S.get(f"{EU}/efetch.fcgi",params={**CM,"db":"nuccore","id":acc,"seq_start":s,"seq_stop":e,
            "rettype":"fasta","retmode":"text"},timeout=600)
    r.raise_for_status()
    body="".join(r.text.splitlines()[1:])
    hdr=f">{lab} {acc}:{s}-{e} {desc}"
    f.write_text(hdr+"\n"+"\n".join(body[i:i+80] for i in range(0,len(body),80))+"\n",encoding="utf-8")
    nN=body.upper().count("N")
    print(f"{lab:12} {acc}:{s}-{e}  {len(body):>9} bp  N={nN} ({100*nN/max(len(body),1):.2f}%)")
    time.sleep(0.5)

# banco de proteinas de referencia
ref=[]
for p in sorted(SEQ.glob("prot_*.fa")):
    lab=p.stem.replace("prot_","")
    txt=p.read_text(encoding="utf-8")
    blocks=[b for b in txt.split(">") if b.strip()]
    if not blocks: print("VAZIO:",lab); continue
    best=max(blocks,key=lambda b: len("".join(b.splitlines()[1:])))
    acc=best.splitlines()[0].split()[0]
    sq="".join(best.splitlines()[1:])
    ref.append(f">{lab}|{acc}|{len(sq)}aa\n"+"\n".join(sq[i:i+80] for i in range(0,len(sq),80)))
(SEQ/"ref_proteinas.fa").write_text("\n".join(ref)+"\n",encoding="utf-8")
print(f"\nbanco de referencia: {len(ref)} proteinas -> {SEQ/'ref_proteinas.fa'}")
for r_ in ref: print("  ",r_.splitlines()[0])
