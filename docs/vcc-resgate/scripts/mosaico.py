import json
from pathlib import Path
import numpy as np, pandas as pd
METAS=Path("data/raw/pie"); VCC=Path("data/raw/vcc/controls-final"); SC=Path(r"C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad")
genes=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist()
ig={g:i for i,g in enumerate(genes)}
def axis(rd):
    m=json.loads((METAS/rd/"meta.json").read_text(encoding="utf-8"))
    return np.array([ig[str(g)] for g in m["genes"] if str(g) in ig])
arc=np.zeros(len(genes),bool); arc[axis("PIE_arc_vcc_25")]=True
jia=np.zeros(len(genes),bool); jia[axis("PIE_jiang")]=True
ori=np.zeros(len(genes),bool); ori[axis("PIE_x_atlas_orion")]=True
pool=np.load(SC/"pool.npy")
pg=pd.read_csv(Path("runs/vcc/ctrl_bulk")/"por_gene.tsv",sep="\t")
univ={c:np.asarray(pg[f"no_universo_{c}"]).astype(bool) for c in "DEF"}
tiers={"arc (peso 2016)":arc,"sem arc, com jiang (peso 143)":(~arc)&jia,"so orion (peso 4)":(~arc)&(~jia)&ori,"nenhum doador":(~arc)&(~jia)&(~ori)}
print(f"{'tier':34s} {'genes':>6s} {'com pool':>8s} {'RMS':>8s} {'p99|.|':>8s} {'max':>7s} {'no univ D':>10s}")
for nome,m in tiers.items():
    s=m&(pool!=0)
    print(f"{nome:34s} {int(m.sum()):6d} {int(s.sum()):8d} {np.sqrt((pool[s]**2).mean()) if s.any() else 0:8.4f} "
          f"{np.percentile(np.abs(pool[s]),99) if s.any() else 0:8.4f} {np.abs(pool[s]).max() if s.any() else 0:7.4f} {int((m&univ['D']).sum()):10d}")
# top movers do delta so-media
z=np.load(Path("runs/vcc/transfer")/"delta_previsto.npz",allow_pickle=True)
alv=[str(a) for a in z["alvos"]]
sem=set(json.loads((Path("runs/vcc/eixo")/"resumo.json").read_text(encoding="utf-8"))["evidencia_doadora"]["alvos_sem_evidencia"])
i0=[i for i,a in enumerate(alv) if a in sem][0]
for c in ("D",):
    d=z[f"delta_{c}"][i0].astype(np.float64)
    o=np.argsort(-np.abs(d))[:20]
    print(f"\ntop 20 |delta| do caminho so-resposta-media em {c}:")
    for i in o:
        tier=next(n for n,m in tiers.items() if m[i])
        print(f"  {genes[i]:12s} delta={d[i]:+.4f} pool={pool[i]:+.4f} univ={bool(univ[c][i])} tier={tier}")
    print(f"\n  negativos entre os 20: {int((d[o]<0).sum())} | fora do universo: {int((~univ[c][o]).sum())}")
