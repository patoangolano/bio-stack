"""Os 877 genes do universo de D fora do eixo do arc: quem estima o deslocamento?"""
import json
from pathlib import Path
import numpy as np, pandas as pd
ARR=Path("data/raw/pie/arrays"); METAS=Path("data/raw/pie"); VCC=Path("data/raw/vcc/controls-final")
SC=Path(r"C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad")
LIM=5.0
genes=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist(); ig={g:i for i,g in enumerate(genes)}
pg=pd.read_csv(Path("runs/vcc/ctrl_bulk")/"por_gene.tsv",sep="\t")
univ=np.asarray(pg["no_universo_D"]).astype(bool)
cob=dict(np.load(SC/"cob.npz")); arc=cob["arc_vcc_25.ARC_H1"]
desl=np.load(SC/"desl_D.npy")
# cobertura jiang somada nos 11 contextos usados
rd="PIE_jiang"; m=json.loads((METAS/rd/"meta.json").read_text(encoding="utf-8"))
gr=[str(g) for g in m["genes"] if str(g) in ig]
eixo={str(g):i for i,g in enumerate(m["genes"])}
cr=np.array([eixo[g] for g in gr]); cp=np.array([ig[g] for g in gr])
cx=np.load(ARR/rd/"preprocessed"/"ctx_ids.npy")
fc=np.load(ARR/rd/"preprocessed"/"fold_changes.npy",mmap_mode="r"); te=np.load(ARR/rd/"preprocessed"/"tested.npy",mmap_mode="r")
usados=["a549_ifnb","a549_tnfa","bxpc3_ifnb","bxpc3_tnfa","hap1_ifnb","hap1_tnfa","ht29_ifnb","ht29_tnfa","k562_ifnb","k562_tnfa","mcf7_tnfa"]
tot=np.zeros(cr.size,dtype=np.int64); nper=0
for ctx in usados:
    lin=np.flatnonzero(cx==int(m["context_to_id"][ctx])); nper+=lin.size
    f=np.asarray(fc[lin])[:,cr].astype(np.float64); t=np.asarray(te[lin])[:,cr]
    with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
    tot+=(t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM)).sum(0)
cobj=np.full(len(genes),-1.0); cobj[cp]=tot/nper
alvo=univ&(arc<0)
print(f"877 genes do universo de D fora do eixo do arc: {int(alvo.sum())}")
print(f"  no eixo do jiang: {int((alvo&(cobj>=0)).sum())} | so no orion: {int((alvo&(cobj<0)).sum())}")
print(f"  perturbacoes jiang somadas nos 11 contextos: {nper}")
for lo,hi in [(0,.1),(.1,.5),(.5,.9),(.9,1.01)]:
    s=alvo&(cobj>=lo)&(cobj<hi)
    if s.sum()==0: continue
    print(f"  cobertura jiang {lo:.1f}-{hi:.2f}: {int(s.sum()):4d} genes | |desl| RMS {np.sqrt((desl[s]**2).mean()):.5f} max {np.abs(desl[s]).max():.4f}")
s=alvo&(cobj<0)
if s.any(): print(f"  so orion (peso 4):        {int(s.sum()):4d} genes | |desl| RMS {np.sqrt((desl[s]**2).mean()):.5f} max {np.abs(desl[s]).max():.4f}")
print("\njiang entrou em alguma das 24 direcoes validadas?",
      any("jiang" in p.stem for p in Path("runs/vcc/transfer").glob("beta_*.npz")))
