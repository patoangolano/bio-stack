"""O coeficiente 0,408 foi medido nos genes em que e aplicado?"""
import json
from pathlib import Path
import numpy as np, pandas as pd
ARR=Path("data/raw/pie/arrays"); METAS=Path("data/raw/pie"); VCC=Path("data/raw/vcc/controls-final")
T=Path("runs/vcc/transfer"); SC=Path(r"C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad")
LIM=5.0
RD={"orion":"PIE_x_atlas_orion","replogle":"PIE_replogle_nadig_essential","arc_vcc_25":"PIE_arc_vcc_25","jiang":"PIE_jiang"}
genes=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist(); ig={g:i for i,g in enumerate(genes)}
def cobertura(repo,ctx):
    """fracao de perturbacoes em que cada gene do painel foi testado, no painel inteiro (-1 = fora do eixo)"""
    rd=RD[repo]; m=json.loads((METAS/rd/"meta.json").read_text(encoding="utf-8"))
    gr=[str(g) for g in m["genes"] if str(g) in ig]
    eixo={str(g):i for i,g in enumerate(m["genes"])}
    cr=np.array([eixo[g] for g in gr]); cp=np.array([ig[g] for g in gr])
    cx=np.load(ARR/rd/"preprocessed"/"ctx_ids.npy")
    lin=np.flatnonzero(cx==int(m["context_to_id"][ctx]))
    fc=np.load(ARR/rd/"preprocessed"/"fold_changes.npy",mmap_mode="r")
    te=np.load(ARR/rd/"preprocessed"/"tested.npy",mmap_mode="r")
    c=np.zeros(cr.size,dtype=np.int64)
    for i0 in range(0,lin.size,512):
        bl=lin[i0:i0+512]
        f=np.asarray(fc[bl])[:,cr].astype(np.float64); t=np.asarray(te[bl])[:,cr]
        with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
        c+=(t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM)).sum(0)
    out=np.full(len(genes),-1.0); out[cp]=c/lin.size
    return out
ctxs=[("orion","hct116"),("orion","hek293t"),("arc_vcc_25","ARC_H1"),
      ("replogle","hepg2"),("replogle","jurkat"),("replogle","k562"),("replogle","rpe1")]
cob={f"{r}.{c}":cobertura(r,c) for r,c in ctxs}
np.savez(SC/"cob.npz", **cob)
print("=== por direcao: o coeficiente foi medido em que genes? ===")
print(f"{'direcao':46s} {'comp':>6s} {'<0.9':>6s} {'r(>=0.9)':>9s} {'r(<0.9)':>9s} {'coef(>=.9)':>10s}")
for f in sorted(T.glob("beta_*.npz")):
    nome=f.stem.removeprefix("beta_")
    doa,alv=nome.split("_para_")
    def key(s):
        for k in cob:
            if s==k.replace(".","_"): return k
        return None
    kd,ka=key(doa),key(alv)
    if kd is None or ka is None: continue
    pk=np.load(f,allow_pickle=True)
    gb=[str(g) for g in pk["genes"]]; idx=np.array([ig[g] for g in gb])
    d=np.asarray(pk["media_por_gene_doador"],float); a=np.asarray(pk["media_por_gene_alvo"],float)
    ok=(d!=0)&(a!=0)
    cd=cob[kd][idx]
    alta=ok&(cd>=0.9); baixa=ok&(cd<0.9)&(cd>=0)
    r_a=np.corrcoef(d[alta],a[alta])[0,1] if alta.sum()>10 else np.nan
    r_b=np.corrcoef(d[baixa],a[baixa])[0,1] if baixa.sum()>10 else np.nan
    c_a=(d[alta]@a[alta])/(d[alta]@d[alta]) if alta.sum()>10 else np.nan
    print(f"{nome:46s} {int(ok.sum()):6d} {int(baixa.sum()):6d} {r_a:+9.4f} {r_b:+9.4f} {c_a:+10.4f}")
