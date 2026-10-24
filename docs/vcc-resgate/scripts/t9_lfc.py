import sys, numpy as np
sys.path.insert(0,"scripts/vcc"); import emissor_400 as E
rng=np.random.default_rng(1); n=500
b=np.log1p(rng.random(n)*50); b[:20]=0.0
lfc=rng.normal(0,0.5,n)                       # log2 da razao
# rota lfc: o leitor multiplica por ln2 e escala_de_gene faz exp()
esc_lfc,_,_=E.escala_de_gene(b, lfc*np.log(2.0), delta_e_log_razao=True)
# rota delta equivalente, construida a mao com a mesma renormalizacao
p=np.expm1(b); alvo=p*np.exp2(lfc); alvo*=p.sum()/alvo.sum()
delta=np.log1p(alvo)-b
esc_delta,_,_=E.escala_de_gene(b, delta, delta_e_log_razao=False)
m=p>0
print("rota lfc == rota delta:", np.allclose(esc_lfc[m],esc_delta[m]),
      "max dif", np.abs(esc_lfc[m]-esc_delta[m]).max())
print("escala == 2^lfc renormalizado:", np.allclose(esc_lfc[m], (np.exp2(lfc)*p.sum()/ (p*np.exp2(lfc)).sum())[m]))
