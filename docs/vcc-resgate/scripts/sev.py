import json
from pathlib import Path
import numpy as np, pandas as pd
SC=Path(r"C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad")
T=Path("runs/vcc/transfer"); VCC=Path("data/raw/vcc/controls-final")
genes=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist(); ig={g:i for i,g in enumerate(genes)}
cob=dict(np.load(SC/"cob.npz"))
pool=np.load(SC/"pool.npz" if (SC/"pool.npz").exists() else SC/"pool.npy")
pg=pd.read_csv(Path("runs/vcc/ctrl_bulk")/"por_gene.tsv",sep="\t")
univ=np.asarray(pg["no_universo_D"]).astype(bool)
desl=np.load(SC/"desl_D.npy")
# cobertura "efetiva": a do doador de maior peso que cobre o gene
arc=cob["arc_vcc_25.ARC_H1"]; ori=cob["orion.hct116"]
cov=np.where(arc>=0, arc, np.nan)   # arc domina com 93% do peso
print("=== dentro do universo de DE de D (%d genes) ===" % univ.sum())
em_arc=univ&(arc>=0)
print(f"  no eixo do arc_vcc_25: {int(em_arc.sum())} de {int(univ.sum())}")
print(f"  fora do eixo do arc  : {int((univ&(arc<0)).sum())}  <- pool vem de jiang/orion")
print(f"\n  {'cobertura no arc':>18s} {'genes':>6s} {'|desl| RMS':>11s} {'|desl| p99':>11s} {'max':>7s} {'massa %':>8s}")
tot=float((desl[univ]**2).sum())
for lo,hi,rot in [(-2,0,'fora do eixo'),(0,.5,'0.00-0.50'),(.5,.9,'0.50-0.90'),(.9,1.01,'0.90-1.00')]:
    m=univ&(arc>=lo)&(arc<hi) if lo>=0 else univ&(arc<0)
    if m.sum()==0: continue
    print(f"  {rot:>18s} {int(m.sum()):6d} {np.sqrt((desl[m]**2).mean()):11.5f} {np.percentile(np.abs(desl[m]),99):11.5f} "
          f"{np.abs(desl[m]).max():7.4f} {100*float((desl[m]**2).sum())/tot:7.1f}%")
# robustez: orion->hek293t restrito a genes bem expressos
print("\n=== robustez: a correlacao negativa orion->orion e efeito de cobertura ou de expressao? ===")
pk=np.load(T/"beta_orion_hct116_para_orion_hek293t.npz",allow_pickle=True)
gb=[str(g) for g in pk["genes"]]; idx=np.array([ig[g] for g in gb])
d=np.asarray(pk["media_por_gene_doador"],float); a=np.asarray(pk["media_por_gene_alvo"],float)
ok=(d!=0)&(a!=0)
ch=cob["orion.hct116"][idx]; ck=cob["orion.hek293t"][idx]
cpm=np.asarray(pg["cpm_cel_D"],float)[idx]
for rot,m in [("tudo (como o script mede)",ok),
              ("cob>=0.9 nos dois",ok&(ch>=0.9)&(ck>=0.9)),
              ("cob>=0.9 nos dois + universo D",ok&(ch>=0.9)&(ck>=0.9)&np.asarray(pg['no_universo_D']).astype(bool)[idx]),
              ("cob>=0.99 nos dois + CPM>20",ok&(ch>=0.99)&(ck>=0.99)&(cpm>20))]:
    if m.sum()<10: continue
    r=np.corrcoef(d[m],a[m])[0,1]; c=(d[m]@a[m])/(d[m]@d[m])
    sse=float(((c*d[m]-a[m])**2).sum()); sz=float((a[m]**2).sum())
    print(f"  {rot:34s} n={int(m.sum()):6d} r={r:+.4f} coef={c:+.4f} ganho={1-sse/sz:+.4f}")
print("\n=== o mesmo para arc->orion (os dois doadores que constroem o pool) ===")
pk=np.load(T/"beta_arc_vcc_25_ARC_H1_para_orion_hct116.npz",allow_pickle=True)
gb=[str(g) for g in pk["genes"]]; idx=np.array([ig[g] for g in gb])
d=np.asarray(pk["media_por_gene_doador"],float); a=np.asarray(pk["media_por_gene_alvo"],float)
ok=(d!=0)&(a!=0); ca=cob["arc_vcc_25.ARC_H1"][idx]; ch=cob["orion.hct116"][idx]
u=np.asarray(pg["no_universo_D"]).astype(bool)[idx]
for rot,m in [("tudo",ok),("cob>=0.9 nos dois",ok&(ca>=0.9)&(ch>=0.9)),("cob>=0.9 + universo D",ok&(ca>=0.9)&(ch>=0.9)&u)]:
    if m.sum()<10: continue
    r=np.corrcoef(d[m],a[m])[0,1]; c=(d[m]@a[m])/(d[m]@d[m])
    sse=float(((c*d[m]-a[m])**2).sum()); sz=float((a[m]**2).sum())
    print(f"  {rot:34s} n={int(m.sum()):6d} r={r:+.4f} coef={c:+.4f} ganho={1-sse/sz:+.4f}")
