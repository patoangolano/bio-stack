"""Fecha o ciclo: delta -> escala -> 400 celulas -> pseudobulk -> b_previsto."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "scripts/vcc")
import emissor_400 as E

genes, alvos = E.ler_eixos()
n = len(genes)
with np.load(E.CTRL_BULK, allow_pickle=True) as z:
    b_ctrl = np.asarray(z["b_ctrl_D"], dtype=np.float64).ravel()
with np.load("runs/vcc/transfer/delta_previsto.npz", allow_pickle=True) as z:
    gnpz = [str(x) for x in np.asarray(z["genes"]).ravel()]
    anpz = [str(x) for x in np.asarray(z["alvos"]).ravel()]
    print("npz genes == painel na ordem:", gnpz == genes)
    print("npz alvos == pert_counts na ordem:", anpz == alvos)
    dD = np.asarray(z["delta_D"], dtype=np.float64)
print("delta_D forma", dD.shape, "zeros por linha min", int((np.abs(dD).sum(axis=1)==0).sum()))

ctrl = E.carregar_controle("D", n)
print("pool D:", E.n_linhas(ctrl), "celulas,", ctrl.nnz, "nnz")

# escolhe alvos com delta grande e pequeno
norma = np.abs(dD).sum(axis=1)
escolhidos = [int(np.argmax(norma)), int(np.argmin(norma)), 0, 150]
for i in escolhidos:
    delta = dD[i]
    escala, inalc, emax = E.escala_de_gene(b_ctrl, delta, delta_e_log_razao=False)
    # 1) escala deve ser exatamente 2^log2r implicito: p_alvo/p_ctrl
    razao_esperada = np.where(np.expm1(b_ctrl) > 0, np.expm1(b_ctrl + delta) / np.maximum(np.expm1(b_ctrl), 1e-300), 0.0)
    # renormalizacao: expm1(b_ctrl+delta) ja soma 5e4 (transfer_fit renormalizou)
    s_alvo = np.expm1(b_ctrl + delta).sum()
    # 2) emite muitas celulas e mede o pseudobulk
    rng = np.random.default_rng([0, 0, i])
    bloco = E.emitir_bloco(ctrl, escala, rng, 400)
    pb = np.asarray(bloco.sum(axis=0), dtype=np.float64).ravel()
    b_emitido = np.log1p(pb / pb.sum() * E.ESCALA_PSEUDOBULK)
    b_alvo = b_ctrl + delta
    m = np.expm1(b_ctrl) > 0
    err = b_emitido[m] - b_alvo[m]
    print(f"alvo {alvos[i]:>10s} |delta|_1={norma[i]:.3f} soma_expm1(b+d)={s_alvo:.3f} "
          f"escala[min,max]=[{escala[m].min():.4f},{escala.max():.4f}] "
          f"max|escala-razao|={np.abs(escala-razao_esperada).max():.3e} "
          f"vies_medio_b={err.mean():+.5f} rms_b={np.sqrt((err**2).mean()):.5f}")
