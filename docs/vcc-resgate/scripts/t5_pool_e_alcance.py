import sys, numpy as np, pandas as pd
sys.path.insert(0, "scripts/vcc")
import emissor_400 as E

genes, alvos = E.ler_eixos(); n = len(genes)
with np.load(E.CTRL_BULK, allow_pickle=True) as z:
    b = np.asarray(z["b_ctrl_D"], dtype=np.float64).ravel()
with np.load("runs/vcc/transfer/delta_previsto.npz", allow_pickle=True) as z:
    dD = np.asarray(z["delta_D"], dtype=np.float64)
ctrl = E.carregar_controle("D", n)

# A) o pool sobrevive a muitos blocos? (eliminate_zeros muta in-place)
assin0 = (int(ctrl.nnz), float(ctrl.data.sum()), int(ctrl.indices.sum()), int(ctrl.indptr[-1]))
esc0,_,_ = E.escala_de_gene(b, dD[0], delta_e_log_razao=False)
ref = E.emitir_bloco(ctrl, esc0, np.random.default_rng([0,0,0]), 400)
alvo_ref = (int(ref.nnz), float(ref.data.sum()))
for i in range(1, 25):
    e,_,_ = E.escala_de_gene(b, dD[i], delta_e_log_razao=False)
    E.emitir_bloco(ctrl, e, np.random.default_rng([0,0,i]), 400)
assin1 = (int(ctrl.nnz), float(ctrl.data.sum()), int(ctrl.indices.sum()), int(ctrl.indptr[-1]))
dep = E.emitir_bloco(ctrl, esc0, np.random.default_rng([0,0,0]), 400)
print("pool intacto apos 25 blocos:", assin0 == assin1, assin0, assin1)
print("bloco reproduzivel (mesma semente, depois de 25 blocos):",
      alvo_ref == (int(dep.nnz), float(dep.data.sum())), alvo_ref, (int(dep.nnz), float(dep.data.sum())))

# B) massa inalcancavel REAL: p_alvo em genes ausentes das 400 celulas sorteadas
print("\nmassa que o relatorio mede (p_ctrl==0) vs massa realmente inalcancavel:")
for i in (0, 150, int(np.abs(dD).sum(axis=1).argmax())):
    delta = dD[i]
    esc, inalc_relatado, _ = E.escala_de_gene(b, delta, delta_e_log_razao=False)
    p_alvo = np.expm1(b + delta); p_alvo = np.maximum(p_alvo,0.0); p_alvo /= p_alvo.sum()
    rng = np.random.default_rng([0,0,i])
    fonte = rng.integers(0, E.n_linhas(ctrl), size=400)
    presentes = np.zeros(n, dtype=bool)
    presentes[np.unique(ctrl[fonte].indices)] = True
    real = float(p_alvo[~presentes].sum())
    print(f"  {alvos[i]:>10s}: relatado={inalc_relatado:.6f}  real={real:.6f}  "
          f"genes ausentes do bloco={int((~presentes).sum())} de {n}")
