import h5py
import numpy as np
p = "C:/Users/Manec/bio-stack/data/raw/vcc/controls-final/context_D.h5ad"
with h5py.File(p, "r") as f:
    def walk(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(f"  D {name}: shape={obj.shape} dtype={obj.dtype}")
        else:
            print(f"  G {name}: attrs={dict(obj.attrs)}")
    f.visititems(walk)
    print("root attrs:", dict(f.attrs))
    print("X attrs:", dict(f["X"].attrs))
    obs = f["obs"]
    print("obs keys:", list(obs.keys()), "attrs:", dict(obs.attrs))
    for k in obs.keys():
        o = obs[k]
        if isinstance(o, h5py.Dataset):
            print(f"   obs/{k}: {o.shape} {o.dtype} sample={o[:3]}")
    print("var keys:", list(f["var"].keys()), dict(f["var"].attrs))
    print("X data sample:", f["X/data"][:8], f["X/indptr"][:5])
