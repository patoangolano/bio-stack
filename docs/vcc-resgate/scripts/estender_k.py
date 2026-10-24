"""Estende a varredura para k=6 e k=8 sem editar scripts/vcc/varredura_k.py.

O otimo nao ficou entre parenteses: `avg_score` ainda subia em k=4 (0,6467), e
so `expr_mse` virou (pico em k=3). Sem um k acima do pico de `avg_score` a
resposta e "4 e o melhor dos testados", que e mais fraca que "o otimo esta em
k≈N".

Monkeypatch em vez de edicao: a constante do autor fica como esta, e a tupla
estendida INCLUI os cinco originais, logo `varredura.json` e
`contagem_por_k.tsv` saem cobrindo os sete em vez de perder os cinco. Os npz dos
cinco sao reescritos com o mesmo conteudo (a decomposicao e deterministica).
"""

import sys
from pathlib import Path

RAIZ = Path("C:/Users/Manec/bio-stack")
sys.path.insert(0, str(RAIZ / "scripts" / "vcc"))

import varredura_k  # noqa: E402

varredura_k.ESCALAS_DESVIO = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0)
print("escalas:", varredura_k.ESCALAS_DESVIO)
varredura_k.main()
