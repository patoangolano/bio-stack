"""HES7 em Gallus: sintenia humana -> regiao do frango -> busca com HES7 de outras aves."""
from __future__ import annotations
import json, time
from pathlib import Path
import requests
BASE="https://rest.ensembl.org"; EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT=Path(r"C:\Users\Manec\bio-stack\runs\atavismo"); SEQ=OUT/"seqs"
CM={"tool":"bio-stack","email":"limatheus458@gmail.com"}; S=requests.Session()
def jget(path,**p):
    S.headers.update({"Content-Type":"application/json"})
    for a in range(6):
        try: r=S.get(BASE+path,params=p,timeout=120)
        except requests.RequestException: time.sleep(3*(a+1)); continue
        if r.status_code in (429,500,502,503,504): time.sleep(3*(a+1)); continue
        if r.status_code==400: return None
        r.raise_for_status(); return r.json()
    return None

print("=== Homo sapiens chr17:8.00-8.25 Mb (vizinhanca de HES7) ===")
hs=sorted([g for g in jget("/overlap/region/homo_sapiens/17:8000000-8250000",feature="gene")
           if g.get("biotype")=="protein_coding"],key=lambda g:g["start"])
for g in hs: print(f"  {g['start']:>10}-{g['end']:<10} {g['strand']:+d} {g.get('external_name')}")

viz=[g["external_name"] for g in hs if g.get("external_name")]
print("\n=== ortologos em Gallus ===")
gg=[]
for sym in viz:
    lk=jget(f"/lookup/symbol/gallus_gallus/{sym}")
    if lk:
        gg.append((sym,lk["seq_region_name"],lk["start"],lk["end"],lk["strand"],lk["id"]))
        print(f"  {sym:12} -> chr{lk['seq_region_name']}:{lk['start']}-{lk['end']} ({lk['strand']:+d})")
    else:
        print(f"  {sym:12} -> sem simbolo em Gallus")
# regiao do frango que cobre os vizinhos sintenicos
cands=[x for x in gg if x[1] in {c for _,c,*_ in gg}]
from collections import Counter
chrom=Counter(x[1] for x in gg).most_common(1)[0][0]
sel=[x for x in gg if x[1]==chrom]
lo=min(x[2] for x in sel)-60_000; hi=max(x[3] for x in sel)+60_000
print(f"\n-> regiao sintenica em Gallus: chr{chrom}:{lo}-{hi}  ({(hi-lo)/1000:.0f} kb)")
allg=sorted(jget(f"/overlap/region/gallus_gallus/{chrom}:{lo}-{hi}",feature="gene"),key=lambda g:g["start"])
print("   conteudo anotado:")
for g in allg:
    print(f"     {g['start']:>10}-{g['end']:<10} {g['strand']:+d} {g.get('external_name') or g['id']:<22} {g.get('biotype')}")
# baixar a regiao (usar acesso RefSeq do cromossomo via Ensembl seq)
r=S.get(f"{BASE}/sequence/region/gallus_gallus/{chrom}:{lo}-{hi}",
        headers={"Content-Type":"text/x-fasta"},timeout=300); r.raise_for_status()
body="".join(r.text.splitlines()[1:])
(SEQ/"reg_GG_hes7.fa").write_text(f">GG_hes7 chr{chrom}:{lo}-{hi} (Ensembl GRCg7b)\n"+
  "\n".join(body[i:i+80] for i in range(0,len(body),80))+"\n",encoding="utf-8")
print(f"   sequencia salva: {len(body)} bp, N={body.upper().count('N')}")
json.dump({"chrom":chrom,"lo":lo,"hi":hi,"vizinhos_gallus":gg,
           "conteudo":[{k:g.get(k) for k in ('start','end','strand','external_name','id','biotype')} for g in allg]},
          open(OUT/"s14_hes7.json","w",encoding="utf-8"),indent=2,ensure_ascii=False)
