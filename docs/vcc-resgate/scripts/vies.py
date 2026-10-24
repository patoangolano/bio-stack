"""A media por gene de log2(fold_changes) mascarada por `tested` e enviesada
positivamente pela propria selecao? Teste: media contra cobertura."""
import json
from pathlib import Path
import numpy as np, pandas as pd
ARR=Path("data/raw/pie/arrays"); METAS=Path("data/raw/pie"); VCC=Path("data/raw/vcc/controls-final")
LIM=5.0; KM=20.0
genes=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist(); ig={g:i for i,g in enumerate(genes)}
def medir(rd,ctx):
    m=json.loads((METAS/rd/"meta.json").read_text(encoding="utf-8"))
    gr=[str(g) for g in m["genes"] if str(g) in ig]
    eixo={str(g):i for i,g in enumerate(m["genes"])}
    cr=np.array([eixo[g] for g in gr])
    cx=np.load(ARR/rd/"preprocessed"/"ctx_ids.npy")
    lin=np.flatnonzero(cx==int(m["context_to_id"][ctx]))
    fc=np.load(ARR/rd/"preprocessed"/"fold_changes.npy",mmap_mode="r")
    te=np.load(ARR/rd/"preprocessed"/"tested.npy",mmap_mode="r")
    s=np.zeros(cr.size); c=np.zeros(cr.size,dtype=np.int64)
    for i0 in range(0,lin.size,512):
        bl=lin[i0:i0+512]
        f=np.asarray(fc[bl])[:,cr].astype(np.float64); t=np.asarray(te[bl])[:,cr]
        with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
        ok=t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM)
        s+=np.where(ok,v,0.0).sum(0); c+=ok.sum(0)
    bruta=np.where(c>0,s/np.maximum(c,1),0.0)
    return gr,bruta,c,lin.size
for rd,ctx in (("PIE_x_atlas_orion","hct116"),("PIE_arc_vcc_25","ARC_H1"),("PIE_replogle_nadig_essential","k562")):
    gr,bruta,c,n=medir(rd,ctx)
    frac=c/n
    ok=c>0
    print(f"\n=== {rd.removeprefix('PIE_')}.{ctx} | {n} perts | fracao testada mediana {np.median(frac):.3f} ===")
    print(f"  corr(fracao_testada, media_bruta) = {np.corrcoef(frac[ok],bruta[ok])[0,1]:+.4f}")
    print(f"  {'faixa de cobertura':>22s} {'genes':>6s} {'media':>8s} {'mediana':>8s} {'%>0':>6s} {'RMS':>7s}")
    cortes=[(0,.01),(.01,.05),(.05,.2),(.2,.5),(.5,.9),(.9,1.01)]
    for lo,hi in cortes:
        m=ok&(frac>=lo)&(frac<hi)
        if m.sum()<5: continue
        print(f"  {f'{lo:.2f}-{hi:.2f}':>22s} {int(m.sum()):6d} {bruta[m].mean():+8.4f} {np.median(bruta[m]):+8.4f} "
              f"{100*(bruta[m]>0).mean():5.1f}% {np.sqrt((bruta[m]**2).mean()):7.4f}")
    pes=bruta*(c/(c+KM)); mp=ok
    print(f"  pos-encolhimento n/(n+20): media global {pes[mp].mean():+.4f}, %>0 {100*(pes[mp]>0).mean():.1f}%, RMS {np.sqrt((pes[mp]**2).mean()):.4f}")
