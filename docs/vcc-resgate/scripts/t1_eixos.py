import h5py, numpy as np, pandas as pd, json
VCC="data/raw/vcc/controls-final"
genes=pd.read_csv(f"{VCC}/gene_names.csv")["gene_name"].astype(str).tolist()
print("painel:",len(genes),"unicos:",len(set(genes)))
alvos=pd.read_csv(f"{VCC}/pert_counts.csv")["target_gene"].astype(str).tolist()
print("alvos:",len(alvos),"unicos:",len(set(alvos)),"nt presente:", "non-targeting" in alvos)
def txt(g,name):
    o=g[name]
    if isinstance(o,h5py.Group):
        if o.attrs.get("encoding-type","")=="categorical":
            cats=o["categories"][:]; cod=o["codes"][:]
            cats=np.array([c.decode() if isinstance(c,bytes) else str(c) for c in cats],dtype=object)
            return cats[cod]
        return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in o["values"][:]],dtype=object)
    v=o[:]
    return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in v],dtype=object)
for c in ("D","E","F"):
    with h5py.File(f"{VCC}/context_{c}.h5ad","r") as f:
        var=f["var"]; idx=str(var.attrs.get("_index","_index"))
        vn=txt(var,idx).tolist()
        print(f"ctx {c}: n_vars={len(vn)} identico_na_ordem={vn==genes} mesmo_conjunto={set(vn)==set(genes)}")
        if vn!=genes:
            for i,(a,b) in enumerate(zip(vn,genes)):
                if a!=b: print("   primeira dif pos",i,repr(a),"vs",repr(b)); break
        print("   obs cols:", [str(x) for x in np.asarray(f["obs"].attrs.get("column-order",[])).ravel()])
        print("   X enc:", dict(f["X"].attrs) if isinstance(f["X"],h5py.Group) else "dense")
