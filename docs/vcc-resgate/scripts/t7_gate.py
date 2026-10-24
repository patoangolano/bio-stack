import sys, numpy as np, pandas as pd
sys.path.insert(0, "scripts/vcc")
import emissor_400 as E
genes, alvos = E.ler_eixos(); n=len(genes)
with np.load(E.CTRL_BULK, allow_pickle=True) as z:
    b = np.asarray(z["b_ctrl_D"],dtype=np.float64).ravel()
with np.load("runs/vcc/transfer/delta_previsto.npz", allow_pickle=True) as z:
    dD = np.asarray(z["delta_D"],dtype=np.float64)
gate = E.ler_gate(alvos, ("D",))
gateados = [a for a in alvos if not gate[("D",a)]]
print("pares gateados em D:", gateados)
for a in gateados:
    i = alvos.index(a)
    # com --sem-gate (default do emissor NAO e este): usa o delta do modelo
    esc_ok,_,_ = E.escala_de_gene(b, dD[i], delta_e_log_razao=False)
    # com o gate ligado (DEFAULT): delta zerado
    esc_gate,_,_ = E.escala_de_gene(b, np.zeros(n), delta_e_log_razao=False)
    print(f"  {a}: gate LIGADO (default) -> escala identica a 1? "
          f"{np.allclose(esc_gate[np.expm1(b)>0],1.0)} (max desvio {np.abs(esc_gate[np.expm1(b)>0]-1).max():.2e})"
          f" | gate DESLIGADO -> escala[min,max]=[{esc_ok[np.expm1(b)>0].min():.4f},{esc_ok.max():.4f}]")
    print(f"     rms do delta que o gate joga fora: {np.sqrt((dD[i]**2).mean()):.3e}")
