"""Comparacao limpa de tres vias no MESMO conjunto de genes."""
import json
from pathlib import Path
import numpy as np, pandas as pd
SC=Path(r"C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad")
T=Path("runs/vcc/transfer"); VCC=Path("data/raw/vcc/controls-final")
genes=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist(); ig={g:i for i,g in enumerate(genes)}
cob=dict(np.load(SC/"cob.npz"))
pg=pd.read_csv(Path("runs/vcc/ctrl_bulk")/"por_gene.tsv",sep="\t")
univ=np.asarray(pg["no_universo_D"]).astype(bool)
# medias por contexto no eixo do painel, tiradas dos npz
def media(f,qual):
    pk=np.load(T/f,allow_pickle=True); gb=[str(g) for g in pk["genes"]]
    out=np.zeros(len(genes)); out[[ig[g] for g in gb]]=np.asarray(pk[f"media_por_gene_{qual}"],float)
    return out
m_h=media("beta_orion_hct116_para_orion_hek293t.npz","doador")
m_k=media("beta_orion_hct116_para_orion_hek293t.npz","alvo")
m_a=media("beta_arc_vcc_25_ARC_H1_para_orion_hct116.npz","doador")
m_h2=media("beta_arc_vcc_25_ARC_H1_para_orion_hct116.npz","alvo")
print("hct116 identico nos dois npz (no suporte comum):",
      np.allclose(m_h[(m_h!=0)&(m_h2!=0)],m_h2[(m_h!=0)&(m_h2!=0)]))
m=(cob["orion.hct116"]>=0.9)&(cob["orion.hek293t"]>=0.9)&(cob["arc_vcc_25.ARC_H1"]>=0.9)&univ
m&=(m_h!=0)&(m_k!=0)&(m_a!=0)
print(f"\ngenes bem cobertos nos tres contextos E no universo de DE de D: {int(m.sum())}")
def par(x,y,nx,ny):
    r=np.corrcoef(x[m],y[m])[0,1]; c=(x[m]@y[m])/(x[m]@x[m])
    g=1-float(((c*x[m]-y[m])**2).sum())/float((y[m]**2).sum())
    print(f"  {nx:22s} -> {ny:22s} r={r:+.4f} coef={c:+.4f} ganho={g:+.4f}")
par(m_h,m_k,"orion.hct116","orion.hek293t")
par(m_k,m_h,"orion.hek293t","orion.hct116")
par(m_a,m_h,"arc.ARC_H1","orion.hct116")
par(m_h,m_a,"orion.hct116","arc.ARC_H1")
par(m_a,m_k,"arc.ARC_H1","orion.hek293t")
par(m_k,m_a,"orion.hek293t","arc.ARC_H1")
for n,v in (("orion.hct116",m_h),("orion.hek293t",m_k),("arc.ARC_H1",m_a)):
    print(f"  RMS da media em {n:16s} = {np.sqrt((v[m]**2).mean()):.5f}")
