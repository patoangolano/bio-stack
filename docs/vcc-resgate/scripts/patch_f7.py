"""Troca literais obsoletos de transfer_fit.py por leitura dos arquivos medidos."""
from pathlib import Path

alvo = Path("scripts/vcc/transfer_fit.py")
texto = alvo.read_text(encoding="utf-8")
trocas = [
    # F7a: o comentario cita 80 enquanto a constante e 60.
    (
        "# Minimo de perturbacoes em comum para um par entrar na validacao. 80 e baixo",
        "# Minimo de perturbacoes em comum para um par entrar na validacao. 60 e baixo",
    ),
    # F7b: comentario com numeros que envelheceram duas vezes. Agora aponta a fonte.
    (
        """    # A RESPOSTA MEDIA do contexto, que e o zero da escala. Medida em
    # runs/vcc/transfer/resposta_media.json: transfere em 24 de 24 direcoes,
    # correlacao mediana 0,436, ganho mediano 0,164 sobre predizer o controle.""",
        """    # A RESPOSTA MEDIA do contexto, que e o zero da escala. Medida em
    # runs/vcc/transfer/resposta_media.json, de onde os numeros sao LIDOS e nao
    # copiados: duas revisoes acharam este comentario com valores de uma corrida
    # anterior, e a segunda trouxe a propria correcao tambem obsoleta.""",
    ),
    # F7c: o resumo carregava 0,164 fixo e a faixa de uma corrida antiga.
    (
        """        "limite_honesto": (
            "Duas parcelas, de tamanhos muito diferentes. A resposta media do contexto "
            "transferiu em 24 de 24 direcoes, com ganho mediano de 0,164 do erro "
            "quadratico sobre predizer o controle — e a parcela grande, e vale para os "
            "300 alvos.""",
        """        "limite_honesto": (
            "Duas parcelas, de tamanhos muito diferentes. A resposta media do contexto "
            f"transferiu em {int(medida['direcoes_com_ganho_positivo'])} de "
            f"{int(medida['direcoes_totais'])} direcoes, com ganho mediano de "
            f"{float(medida['ganho_mediano_sobre_o_controle']):.3f} do erro "
            "quadratico sobre predizer o controle — e a parcela grande, e vale para os "
            "300 alvos.""",
    ),
    (
        """        "peso_diagonal_contra_programas": {
            "usado": 0.5,
            "otimo_medido_mediano": 0.393,
            "faixa_medida": [0.041, 0.788],
            "por_que_fica_em_meio_a_meio": (
                "O otimo rende 0,0003 de mse/mse0 DENTRO da amostra e a faixa entre "
                "direcoes vai de 0,04 a 0,79 — nao esta determinado. Perseguir isso "
                "custaria uma reemissao inteira por ganho abaixo do ruido."
            ),
        },""",
        """        "peso_diagonal_contra_programas": {
            "usado": 0.5,
            "otimo_medido_mediano": round(float(np.median(pesos_otimos)), 4),
            "faixa_medida": [
                round(float(pesos_otimos.min()), 4),
                round(float(pesos_otimos.max()), 4),
            ],
            "por_que_fica_em_meio_a_meio": (
                "O otimo rende 0,0003 de mse/mse0 DENTRO da amostra e a faixa entre "
                f"direcoes vai de {pesos_otimos.min():.2f} a {pesos_otimos.max():.2f} "
                "— nao esta determinado. Perseguir isso custaria uma reemissao inteira "
                "por ganho abaixo do ruido. Mediana e faixa saem de validacao_loco.json "
                "a cada corrida, nao de literal no codigo."
            ),
        },""",
    ),
]
for antes, depois in trocas:
    if texto.count(antes) != 1:
        raise SystemExit(f"ancora nao unica ({texto.count(antes)}x): {antes[:70]!r}")
    texto = texto.replace(antes, depois)
alvo.write_text(texto, encoding="utf-8")
print("4 trocas aplicadas")
