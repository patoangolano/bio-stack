import json
from pathlib import Path
import numpy as np, pandas as pd
T=Path("runs/vcc/transfer"); EIXO=Path("runs/vcc/eixo"); CTRL=Path("runs/vcc/ctrl_bulk")
z=np.load(T/"delta_previsto.npz", allow_pickle=True)
genes=[str(g) for g in z["genes"]]; alvos=[str(a) for a in z["alvos"]]
print("genes",len(genes),"alvos",len(alvos))
res=json.loads((EIXO/"resumo.json").read_text(encoding="utf-8"))
sem=set(res["evidencia_doadora"]["alvos_sem_evidencia"]); resv=set(res["evidencia_doadora"]["reservados_para_validacao"])
print("sem_evidencia",len(sem),"reservados",len(resv),"reservados subset de sem:",resv<=sem)
gate=pd.read_csv(EIXO/"gate_alvos.tsv", sep="\t")
print("gate linhas",len(gate),"dup (ctx,alvo):",gate.duplicated(["contexto","alvo"]).sum())
for c in ("D","E","F"):
    d=z[f"delta_{c}"]
    nz=np.abs(d).sum(axis=1)
    print(f"\n--- {c} --- shape {d.shape} dtype {d.dtype}")
    print("  alvos com delta exatamente zero:", int((nz==0).sum()))
    print("  nan/inf:", int(np.isnan(d).sum()), int(np.isinf(d).sum()))
    idx_sem=[i for i,a in enumerate(alvos) if a in sem]
    sub=d[idx_sem]
    print(f"  {len(idx_sem)} sem-evidencia: todas iguais? ", bool(np.all(sub==sub[0])))
    print("  |delta| sem-evid norma:", float(np.linalg.norm(sub[0])))
    gc=gate.loc[gate.contexto==c].set_index("alvo")
    gated=[a for a in alvos if a not in sem and str(gc.loc[a,"gate"])!="ativo"]
    print("  gateados com evidencia:", gated)
    for a in gated:
        i=alvos.index(a)
        print("    ", a, "igual ao sem-evid?", bool(np.all(d[i]==sub[0])))
