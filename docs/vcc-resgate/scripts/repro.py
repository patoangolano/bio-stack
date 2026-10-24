"""Reproducao independente de ajustar() para um alvo de doador unico."""
import json
from pathlib import Path
import numpy as np, pandas as pd

A=Path("data/raw/pie/arrays/PIE_x_atlas_orion/preprocessed")
M=Path("data/raw/pie/PIE_x_atlas_orion")
VCC=Path("data/raw/vcc/controls-final"); T=Path("runs/vcc/transfer"); CTRL=Path("runs/vcc/ctrl_bulk")
LIM=5.0; KM=20.0

genes_painel=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist()
alvos=pd.read_csv(VCC/"pert_counts.csv")["target_gene"].astype(str).tolist()
ig={g:i for i,g in enumerate(genes_painel)}
meta=json.loads((M/"meta.json").read_text(encoding="utf-8"))
genes_repo=[str(g) for g in meta["genes"] if str(g) in ig]
eixo={str(g):i for i,g in enumerate(meta["genes"])}
cols_repo=np.array([eixo[g] for g in genes_repo],dtype=np.int64)
colunas=np.array([ig[g] for g in genes_repo],dtype=np.int64)
p2i={str(p):int(i) for p,i in meta["pert_to_id"].items()}
cid=int(meta["context_to_id"]["hct116"])

cx=np.load(A/"ctx_ids.npy"); px=np.load(A/"pert_ids.npy")
lin_ctx=np.flatnonzero(cx==cid)
fc=np.load(A/"fold_changes.npy",mmap_mode="r"); te=np.load(A/"tested.npy",mmap_mode="r")

# media do contexto hct116 sobre TODAS as perts, encolhida por cobertura
soma=np.zeros(cols_repo.size); cont=np.zeros(cols_repo.size,dtype=np.int64)
for i0 in range(0,lin_ctx.size,512):
    b=lin_ctx[i0:i0+512]
    f=np.asarray(fc[b])[:,cols_repo].astype(np.float64)
    t=np.asarray(te[b])[:,cols_repo]
    with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
    ok=t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM)
    soma+=np.where(ok,v,0.0).sum(0); cont+=ok.sum(0)
bruta=np.where(cont>0,soma/np.maximum(cont,1),0.0)
media=bruta*(cont/(cont+KM))
print("media hct116: RMS bruta %.4f -> encolhida %.4f ; cobertura mediana %.0f"%(
    np.sqrt((bruta**2).mean()), np.sqrt((media**2).mean()), np.median(cont)))

# alvo de doador unico (orion.hct116 so)
split=json.loads((Path("data/raw/pie")/"PIE_splits"/"replogle_xdataset"/"train.json").read_text(encoding="utf-8"))
tr=set()
for k,v in split.items():
    if str(k).split(".",1)[0]=="orion": tr.update(str(p) for p in v)
disp_h={str(p) for p in np.unique(px[lin_ctx]).tolist() for _ in [0]}
id2p={int(i):str(p) for p,i in meta["pert_to_id"].items()}
disp_h={id2p[int(i)] for i in np.unique(px[lin_ctx]).tolist()}
cid2=int(meta["context_to_id"]["hek293t"])
disp_k={id2p[int(i)] for i in np.unique(px[np.flatnonzero(cx==cid2)]).tolist()}
cands=sorted((set(alvos)&disp_h&tr)-disp_k)
# tirar os que tambem estao em arc/jiang
multi={"AKT2","MTA1","NIT1","PRCP","IFNAR2","IKBKE"}
cands=[c for c in cands if c not in multi]
ALVO=cands[0]
print("alvo escolhido:",ALVO,"de",len(cands),"candidatos de doador unico")

lrow=int(np.flatnonzero((cx==cid)&(px==p2i[ALVO]))[0])
f=np.asarray(fc[lrow])[cols_repo].astype(np.float64)
t=np.asarray(te[lrow])[cols_repo]
with np.errstate(divide="ignore",invalid="ignore"): v=np.log2(f)
valido=t&(f>0)&np.isfinite(v)&(np.abs(v)<=LIM)
v=np.where(valido,v,0.0)
centrado=np.where(valido,v-media,0.0)

# beta / base / gama da direcao de referencia
pk=np.load(T/"beta_orion_hct116_para_replogle_rpe1.npz",allow_pickle=True)
genes_beta=[str(g) for g in pk["genes"]]; bvals=np.asarray(pk["beta"],dtype=np.float64)
bglob=float(pk["beta_global"]); base=np.asarray(pk["base_programas"],dtype=np.float64)
gama=np.asarray(pk["gama"],dtype=np.float64)
beta=np.full(len(genes_painel),bglob)
for g,val in zip(genes_beta,bvals,strict=True):
    if g in ig: beta[ig[g]]=val
pos={g:i for i,g in enumerate(genes_repo)}
pares=[(pos[g],j) for j,g in enumerate(genes_beta) if g in pos]
lr=np.array([p for p,_ in pares]); lb=np.array([j for _,j in pares])
print("genes de base alinhados:",len(pares),"de",len(genes_repo),"| base shape",base.shape)
basea=base[lb]
diag=beta[colunas]*centrado*valido
carga=(centrado[lr]*valido[lr])@basea
recon=basea@(carga*gama)
num=diag.copy(); den=valido.astype(np.float64).copy()
num[lr]+= 0.5*(recon*valido[lr]-diag[lr])
prev=np.zeros(len(genes_painel))
prev[colunas]=np.where(den>0,num/np.maximum(den,1e-12),0.0)

# deslocamento recuperado do proprio arquivo (via alvos sem evidencia)
z=np.load(T/"delta_previsto.npz",allow_pickle=True)
res=json.loads((Path("runs/vcc/eixo")/"resumo.json").read_text(encoding="utf-8"))
sem=set(res["evidencia_doadora"]["alvos_sem_evidencia"])
pg=pd.read_csv(CTRL/"por_gene.tsv",sep="\t")
for ctx in ("D",):
    b=np.asarray(pg[f"b_ctrl_{ctx}"],dtype=np.float64); prop=np.expm1(b)
    d=z[f"delta_{ctx}"]
    isem=[i for i,a in enumerate(alvos) if a in sem][0]
    u=np.log2(np.expm1(b+d[isem].astype(np.float64))/prop)
    # genes sem estimativa de media -> desl=0 -> u = log2 k constante
    k=np.median(u[np.abs(u-np.median(u))<1e-9]) if False else None
    vals,cnts=np.unique(np.round(u,10),return_counts=True)
    log2k=vals[np.argmax(cnts)]
    print(f"{ctx}: log2k={log2k:+.6f} (moda em {cnts.max()} de {u.size} genes)")
    desl=u-log2k
    # confere o alvo
    tot=desl+prev
    pert=prop*np.exp2(tot); pert*=prop.sum()/pert.sum()
    mine=(np.log1p(pert)-b)
    theirs=d[alvos.index(ALVO)].astype(np.float64)
    err=np.abs(mine-theirs)
    print(f"  {ALVO}: max|diff|={err.max():.3e}  rms={np.sqrt((err**2).mean()):.3e}  "
          f"corr={np.corrcoef(mine,theirs)[0,1]:.8f}  norma_minha={np.linalg.norm(mine):.4f} norma_deles={np.linalg.norm(theirs):.4f}")
    np.save("/c/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/07a4432b-599e-4679-83d3-c71280471116/scratchpad/desl_D.npy", desl)
    print(f"  desl: genes!=0 {int((desl!=0).sum())} | |desl| tipico {np.linalg.norm(desl)/np.sqrt(desl.size):.4f}")
