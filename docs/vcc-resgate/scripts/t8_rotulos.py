import sys, numpy as np, pandas as pd, anndata as ad, h5py
from pathlib import Path
OUT=Path(sys.argv[1])
genes=pd.read_csv("data/raw/vcc/controls-final/gene_names.csv")["gene_name"].astype(str).tolist()
alvos=pd.read_csv("data/raw/vcc/controls-final/pert_counts.csv")["target_gene"].astype(str).tolist()
NG,NC=len(genes),400; nobs=len(alvos)*NC; K=20
rng=np.random.default_rng(7)
def xdata():
    cols=np.zeros((nobs,K),dtype=np.int32)
    for r in range(nobs): cols[r]=np.sort(rng.permutation(NG)[:K])
    return rng.integers(1,50,size=(nobs,K)).astype(np.float32).ravel(), cols.ravel(), np.arange(0,nobs*K+1,K,dtype=np.int64)
def escrever(v,caminho):
    ca,cc=[],[]
    for a in alvos: ca.extend([a]*NC); cc.extend(["D"]*NC)
    if v=="ntc": ca[0:5]=["non-targeting"]*5          # 5 linhas de controle
    if v=="contagem": ca[0:5]=[alvos[1]]*5            # alvos[0] com 395, alvos[1] com 405
    cats=sorted(set(ca))
    obs=pd.DataFrame({"target_gene":pd.Categorical(ca,categories=cats),
                      "context":pd.Categorical(cc,categories=["D"])},
                     index=pd.Index([f"c{i}" for i in range(nobs)],dtype=object))
    ad.AnnData(obs=obs,var=pd.DataFrame(index=pd.Index(genes,dtype=object))).write_h5ad(caminho,compression=None)
    dat,ind,ip=xdata()
    with h5py.File(caminho,"r+") as f:
        if "X" in f: del f["X"]
        g=f.create_group("X"); g.attrs["encoding-type"]="csr_matrix"; g.attrs["encoding-version"]="0.1.0"
        g.attrs["shape"]=np.array([nobs,NG],dtype=np.int64)
        g.create_dataset("data",data=dat); g.create_dataset("indices",data=ind); g.create_dataset("indptr",data=ip)
for v in sys.argv[2:]:
    p=OUT/f"forja_{v}.h5ad"; escrever(v,p); print("escrito",p.name)
