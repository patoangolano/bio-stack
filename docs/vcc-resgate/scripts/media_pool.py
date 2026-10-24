import json
from pathlib import Path
import numpy as np, pandas as pd
ARR=Path("data/raw/pie/arrays"); METAS=Path("data/raw/pie"); VCC=Path("data/raw/vcc/controls-final")
T=Path("runs/vcc/transfer"); SC=Path(r"C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad")
LIM=5.0; KM=20.0
REPOS={"orion":"PIE_x_atlas_orion","replogle":"PIE_replogle_nadig_essential","arc_vcc_25":"PIE_arc_vcc_25","jiang":"PIE_jiang"}
prof={"arc_vcc_25":2016.0,"replogle":46.0,"jiang":13.0,"orion":2.0}
genes_painel=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist()
alvos=set(pd.read_csv(VCC/"pert_counts.csv")["target_gene"].astype(str))
ig={g:i for i,g in enumerate(genes_painel)}
split=json.loads((METAS/"PIE_splits"/"replogle_xdataset"/"train.json").read_text(encoding="utf-8"))
tr={}
for k,v in split.items(): tr.setdefault(str(k).split(".",1)[0],set()).update(str(p) for p in v)
num=np.zeros(len(genes_painel)); den=np.zeros(len(genes_painel))
for repo,rd in REPOS.items():
    meta=json.loads((METAS/rd/"meta.json").read_text(encoding="utf-8")); ds=str(meta["dataset"])
    gr=[str(g) for g in meta["genes"] if str(g) in ig]
    eixo={str(g):i for i,g in enumerate(meta["genes"])}
    cr=np.array([eixo[g] for g in gr]); cp=np.array([ig[g] for g in gr])
    id2p={int(i):str(p) for p,i in meta["pert_to_id"].items()}
    cx=np.load(ARR/rd/"preprocessed"/"ctx_ids.npy"); px=np.load(ARR/rd/"preprocessed"/"pert_ids.npy")
    fc=np.load(ARR/rd/"preprocessed"/"fold_changes.npy",mmap_mode="r")
    te=np.load(ARR/rd/"preprocessed"/"tested.npy",mmap_mode="r")
    for ctx,cidv in sorted(meta["context_to_id"].items()):
        lin=np.flatnonzero(cx==int(cidv))
        disp={id2p[int(i)] for i in np.unique(px[lin]).tolist()}
        if not (alvos & disp & tr.get(ds,set())): continue
        s=np.zeros(cr.size); c=np.zeros(cr.size,dtype=np.int64)
        for i0 in range(0,lin.size,512):
            bl=lin[i0:i0+512]
            f=np.asarray(fc[bl])[:,cr].astype(np.float64); t=np.asarray(te[bl])[:,cr]
            with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
            ok=t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM)
            s+=np.where(ok,v,0.0).sum(0); c+=ok.sum(0)
        m=np.where(c>0,s/np.maximum(c,1),0.0)*(c/(c+KM))
        tem=(m!=0.0).astype(np.float64)
        num[cp]+=prof[repo]*m*tem; den[cp]+=prof[repo]*tem
        print(f"  {repo}.{ctx:12s} perts={lin.size:5d} RMS m={np.sqrt((m**2).mean()):.4f} genes com media={int(tem.sum())}")
pool=np.where(den>0,num/np.maximum(den,1e-12),0.0)
rm=json.loads((T/"resposta_media.json").read_text(encoding="utf-8"))
coef=float(rm["coeficiente_mediano"])
print(f"\ncoeficiente_mediano={coef:.6f} usar={rm['usar_resposta_media_transferida']}")
print(f"pool: genes!=0 {int((pool!=0).sum())} RMS {np.sqrt((pool**2).mean()):.5f} max|.| {np.abs(pool).max():.4f}")
mine=coef*pool
np.save(SC/"pool.npy",pool)
for ctx in ("D","E","F"):
    desl=np.load(SC/f"desl_{ctx}.npy")
    sup=desl!=0.0
    d=np.abs(mine[sup]-desl[sup])
    print(f"{ctx}: suporte {int(sup.sum())} | max|diff| {d.max():.3e} rms {np.sqrt((d**2).mean()):.3e}")
    fora=(~sup)&(np.abs(mine)>1e-9)
    print(f"   genes com desl=0 no arquivo mas mine!=0: {int(fora.sum())} (b_ctrl==0, delta forcado a 0)")
