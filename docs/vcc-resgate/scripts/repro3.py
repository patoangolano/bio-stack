import json
from pathlib import Path
import numpy as np, pandas as pd
A=Path("data/raw/pie/arrays/PIE_x_atlas_orion/preprocessed"); M=Path("data/raw/pie/PIE_x_atlas_orion")
VCC=Path("data/raw/vcc/controls-final"); T=Path("runs/vcc/transfer"); CTRL=Path("runs/vcc/ctrl_bulk")
SC=Path(r"C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad")
LIM=5.0; KM=20.0
genes_painel=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist()
alvos=pd.read_csv(VCC/"pert_counts.csv")["target_gene"].astype(str).tolist()
ig={g:i for i,g in enumerate(genes_painel)}
meta=json.loads((M/"meta.json").read_text(encoding="utf-8"))
genes_repo=[str(g) for g in meta["genes"] if str(g) in ig]
eixo={str(g):i for i,g in enumerate(meta["genes"])}
cols_repo=np.array([eixo[g] for g in genes_repo]); colunas=np.array([ig[g] for g in genes_repo])
p2i={str(p):int(i) for p,i in meta["pert_to_id"].items()}; id2p={v:k for k,v in p2i.items()}
cid=int(meta["context_to_id"]["hct116"]); cid2=int(meta["context_to_id"]["hek293t"])
cx=np.load(A/"ctx_ids.npy"); px=np.load(A/"pert_ids.npy")
lin=np.flatnonzero(cx==cid)
fc=np.load(A/"fold_changes.npy",mmap_mode="r"); te=np.load(A/"tested.npy",mmap_mode="r")
soma=np.zeros(cols_repo.size); cont=np.zeros(cols_repo.size,dtype=np.int64)
for i0 in range(0,lin.size,512):
    bl=lin[i0:i0+512]
    f=np.asarray(fc[bl])[:,cols_repo].astype(np.float64); t=np.asarray(te[bl])[:,cols_repo]
    with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
    ok=t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM)
    soma+=np.where(ok,v,0.0).sum(0); cont+=ok.sum(0)
media=np.where(cont>0,soma/np.maximum(cont,1),0.0)*(cont/(cont+KM))
np.save(SC/"media_hct116.npy",media); np.save(SC/"cols_repo_orion.npy",colunas)
split=json.loads((Path("data/raw/pie")/"PIE_splits"/"replogle_xdataset"/"train.json").read_text(encoding="utf-8"))
tr={str(p) for k,v in split.items() if str(k).split(".",1)[0]=="orion" for p in v}
disp_h={id2p[int(i)] for i in np.unique(px[lin]).tolist()}
disp_k={id2p[int(i)] for i in np.unique(px[np.flatnonzero(cx==cid2)]).tolist()}
multi={"AKT2","MTA1","NIT1","PRCP","IFNAR2","IKBKE"}
cands=[c for c in sorted((set(alvos)&disp_h&tr)-disp_k) if c not in multi]
pk=np.load(T/"beta_orion_hct116_para_replogle_rpe1.npz",allow_pickle=True)
genes_beta=[str(g) for g in pk["genes"]]; bvals=np.asarray(pk["beta"],float)
bglob=float(pk["beta_global"]); base=np.asarray(pk["base_programas"],float); gama=np.asarray(pk["gama"],float)
beta=np.full(len(genes_painel),bglob)
for g,val in zip(genes_beta,bvals,strict=True):
    if g in ig: beta[ig[g]]=val
pos={g:i for i,g in enumerate(genes_repo)}
pares=[(pos[g],j) for j,g in enumerate(genes_beta) if g in pos]
lr=np.array([p for p,_ in pares]); lb=np.array([j for _,j in pares]); basea=base[lb]
z=np.load(T/"delta_previsto.npz",allow_pickle=True)
sem=set(json.loads((Path("runs/vcc/eixo")/"resumo.json").read_text(encoding="utf-8"))["evidencia_doadora"]["alvos_sem_evidencia"])
isem=[i for i,a in enumerate(alvos) if a in sem][0]
pg=pd.read_csv(CTRL/"por_gene.tsv",sep="\t")
for ctx in ("D","E","F"):
    b=np.asarray(pg[f"b_ctrl_{ctx}"],float); prop=np.expm1(b); pos_=prop>0
    d=z[f"delta_{ctx}"].astype(np.float64)
    u=np.full(prop.size,np.nan); u[pos_]=np.log2(np.expm1(b[pos_]+d[isem][pos_])/prop[pos_])
    uu=u[pos_]; vals,cnts=np.unique(np.round(uu,9),return_counts=True); log2k=vals[np.argmax(cnts)]
    desl=np.where(pos_,u-log2k,0.0)
    print(f"{ctx}: log2k={log2k:+.6f} moda {cnts.max()}/{uu.size} | desl!=0 {int((np.abs(desl)>1e-9).sum())} | |desl| rms {np.sqrt((desl**2).mean()):.5f} max {np.abs(desl).max():.4f}")
    np.save(SC/f"desl_{ctx}.npy",desl)
    if ctx!="D": continue
    piores=[]
    for ALVO in cands[:6]:
        lrow=int(np.flatnonzero((cx==cid)&(px==p2i[ALVO]))[0])
        f=np.asarray(fc[lrow])[cols_repo].astype(np.float64); t=np.asarray(te[lrow])[cols_repo]
        with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
        val=t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM); v=np.where(val,v,0.0)
        cen=np.where(val,v-media,0.0)
        diag=beta[colunas]*cen*val
        carga=(cen[lr]*val[lr])@basea; recon=basea@(carga*gama)
        num=diag.copy(); den=val.astype(float).copy(); num[lr]+=0.5*(recon*val[lr]-diag[lr])
        prev=np.zeros(len(genes_painel)); prev[colunas]=np.where(den>0,num/np.maximum(den,1e-12),0.0)
        tot=desl+prev; pert=prop*np.exp2(tot); pert*=prop.sum()/pert.sum()
        mine=np.log1p(pert)-b; theirs=d[alvos.index(ALVO)]
        err=np.abs(mine-theirs)
        piores.append((ALVO,err.max(),np.sqrt((err**2).mean()),np.linalg.norm(mine),np.linalg.norm(theirs)))
    for a,mx,rms,nm,nt in piores:
        print(f"   {a:10s} max|diff|={mx:.3e} rms={rms:.3e} |mine|={nm:.4f} |theirs|={nt:.4f}")
