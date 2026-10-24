import json, numpy as np, pandas as pd
alvos = pd.read_csv("data/raw/vcc/controls-final/pert_counts.csv")["target_gene"].astype(str).tolist()
res = json.load(open("runs/vcc/eixo/resumo.json", encoding="utf-8"))
sem = set(res["evidencia_doadora"]["alvos_sem_evidencia"]); resv = set(res["evidencia_doadora"]["reservados_para_validacao"])
gate = pd.read_csv("runs/vcc/eixo/gate_alvos.tsv", sep="\t")
with np.load("runs/vcc/transfer/delta_previsto.npz", allow_pickle=True) as z:
    for c in ("D","E","F"):
        d = np.asarray(z[f"delta_{c}"], dtype=np.float64)
        norma = np.abs(d).sum(axis=1); rms = np.sqrt((d**2).mean(axis=1))
        zeros = int((norma == 0).sum())
        g = gate[gate.contexto==c].set_index("alvo")["gate"]
        inativos = [a for a in alvos if str(g.loc[a])!="ativo"]
        idx_sem = [i for i,a in enumerate(alvos) if a in sem or a in resv]
        idx_com = [i for i,a in enumerate(alvos) if a not in sem and a not in resv]
        # linhas identicas entre si (os sem evidencia devem compartilhar o mesmo delta)
        uniq = np.unique(d[idx_sem], axis=0).shape[0]
        print(f"{c}: linhas exatamente zero={zeros} | rms|delta| min={rms.min():.3e} "
              f"p50={np.median(rms):.3e} | sem_evid/reserv n={len(idx_sem)} linhas distintas={uniq} "
              f"rms={rms[idx_sem][0]:.3e} | com_evid rms p50={np.median(rms[idx_com]):.3e} "
              f"| gate inativo={inativos}")
        # o alvo gateado-com-evidencia deve ter o MESMO delta dos sem evidencia
        for a in inativos:
            i = alvos.index(a)
            igual = np.array_equal(d[i], d[idx_sem[0]])
            print(f"    gateado {a}: delta identico ao de 'so resposta media'? {igual} rms={rms[i]:.3e}")
