"""
A resposta MEDIA do contexto as perturbacoes — e nao o controle — e o zero da
escala do VCC 2026. Este script mede quanto dela se transfere entre contextos,
porque disso depende se os 143 alvos sem evidencia pontuam zero ou abaixo de
zero.

## O erro que este script existe para corrigir

A escala do desafio e `s = (u - b) / (r - b)`, onde `b` e o valor da metrica
alcancado por uma predicao que devolve, para toda perturbacao, a **media real
das 300 respostas daquele contexto**. Predizer essa media pontua exatamente
zero. Predizer o CONTROLE (delta = 0) pontua **abaixo** de zero, porque o
controle esta mais longe de cada resposta real do que a media delas esta — a
media carrega o componente comum a todas as 300, que num rastreio de CRISPRi
inclui programa de estresse e de proliferacao.

O `transfer_fit` prediz o desvio CENTRADO de cada perturbacao. Centrado em
relacao a media do contexto doador, que foi subtraida e descartada. Na emissao,
`delta = 0` virou "prediga o controle". Para os 143 alvos sem evidencia isso e
a predicao inteira, e ela pontua negativo em vez de zero.

A correcao e somar uma estimativa da resposta media de D/E/F a todos os 300
alvos. Para os 157 com evidencia: media mais desvio. Para os 143: media apenas,
que e o zero honesto.

## O que da para medir, e o que nao da

Nao ha uma unica celula perturbada em D, E ou F, entao a resposta media deles e
inestimavel diretamente. O que se mede e **quanto a resposta media de um
contexto prediz a de outro**, nas 24 direcoes que o `validar` ja percorreu —
`media_por_gene_doador` e `media_por_gene_alvo` ficaram gravadas em cada
`beta_*.npz`.

Se a resposta media transferir bem (correlacao alta, coeficiente perto de 1),
usa-la e seguro e o ganho e grande: tira 143 alvos de baixo de zero. Se
transferir mal, somar uma estimativa ruim pode piorar os 300 em vez de
melhorar 143 — e aí a decisao muda.

Rodar da raiz do projeto, depois de `transfer_fit.py validar`:
    ./.venv/Scripts/python.exe -I scripts/vcc/resposta_media.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

TRANSFER = Path("runs/vcc/transfer")
SAIDA = Path("runs/vcc/transfer")


def direcoes_da_corrida_atual() -> list[Path]:
    """Os beta_*.npz que a corrida ATUAL de `validar` produziu.

    Lista tirada de `validacao_loco.json`, nao de um glob do diretorio. O glob
    pegava sobra de corrida anterior com outro esquema de nome e outra variavel
    transferida: a primeira medicao desta funcao viu 31 direcoes onde existem
    24, e as 7 extras vinham de um pipeline que nao existe mais.
    """
    relatorio = TRANSFER / "validacao_loco.json"
    if not relatorio.exists():
        raise FileNotFoundError("rode `transfer_fit.py validar` antes: sem validacao_loco.json")
    pares = json.loads(relatorio.read_text(encoding="utf-8"))["pares"]
    caminhos: list[Path] = []
    for rotulo in pares:
        doador, alvo = str(rotulo).split("->")
        nome = (
            f"beta_{doador.replace('.', '_')}_para_{alvo.replace('.', '_')}.npz"
        )
        caminho = TRANSFER / nome
        if not caminho.exists():
            raise FileNotFoundError(f"{nome} declarado em validacao_loco.json e ausente")
        caminhos.append(caminho)
    return caminhos


def main() -> None:
    arquivos = direcoes_da_corrida_atual()
    print(f"== {len(arquivos)} direcoes da corrida atual ==")

    linhas: list[dict[str, float | str]] = []
    for caminho in arquivos:
        pacote = np.load(caminho, allow_pickle=True)
        doador = np.asarray(pacote["media_por_gene_doador"], dtype=np.float64)
        alvo = np.asarray(pacote["media_por_gene_alvo"], dtype=np.float64)
        # Genes onde as duas medias existem e nao sao exatamente zero (zero aqui
        # significa gene sem nenhuma medida valida no contexto, nao media nula).
        ok = (doador != 0.0) & (alvo != 0.0) & np.isfinite(doador) & np.isfinite(alvo)
        if ok.sum() < 500:
            print(f"  {caminho.stem}: so {int(ok.sum())} genes comparaveis — fora")
            continue
        d, a = doador[ok], alvo[ok]
        correlacao = float(np.corrcoef(d, a)[0, 1])
        # Coeficiente sem intercepto: quanto da media do alvo a do doador explica.
        coeficiente = float((d @ a) / max(float(d @ d), 1e-12))
        # Fracao do erro quadratico que usar a media do doador remove, contra
        # usar zero (o controle). Negativo significa que usar piora.
        sse_doador = float(np.square(coeficiente * d - a).sum())
        sse_zero = float(np.square(a).sum())
        ganho = 1.0 - sse_doador / max(sse_zero, 1e-12)
        linhas.append(
            {
                "direcao": caminho.stem.removeprefix("beta_"),
                "genes": float(ok.sum()),
                "correlacao": correlacao,
                "coeficiente": coeficiente,
                "ganho_sobre_o_controle": ganho,
                "norma_media_doador": float(np.linalg.norm(d) / np.sqrt(d.size)),
                "norma_media_alvo": float(np.linalg.norm(a) / np.sqrt(a.size)),
            }
        )
        print(
            f"  {caminho.stem.removeprefix('beta_'):52s} r={correlacao:+.4f} "
            f"coef={coeficiente:+.4f} ganho={ganho:+.4f} "
            f"|media| doador={linhas[-1]['norma_media_doador']:.4f} "
            f"alvo={linhas[-1]['norma_media_alvo']:.4f}"
        )

    if not linhas:
        raise ValueError("nenhuma direcao com genes comparaveis suficientes")

    correlacoes = np.array([float(x["correlacao"]) for x in linhas])
    coeficientes = np.array([float(x["coeficiente"]) for x in linhas])
    ganhos = np.array([float(x["ganho_sobre_o_controle"]) for x in linhas])
    normas = np.array([float(x["norma_media_alvo"]) for x in linhas])

    print("\n== resumo entre direcoes ==")
    print(f"  correlacao da resposta media:  mediana {np.median(correlacoes):+.4f} "
          f"min {correlacoes.min():+.4f} max {correlacoes.max():+.4f}")
    print(f"  coeficiente:                   mediana {np.median(coeficientes):+.4f} "
          f"min {coeficientes.min():+.4f} max {coeficientes.max():+.4f}")
    print(f"  ganho sobre predizer o controle: mediana {np.median(ganhos):+.4f} "
          f"min {ganhos.min():+.4f} max {ganhos.max():+.4f}")
    print(f"  direcoes com ganho positivo: {int((ganhos > 0).sum())} de {len(ganhos)}")
    print(f"  magnitude tipica da resposta media (log2): {np.median(normas):.4f}")

    usar = bool(np.median(ganhos) > 0 and (ganhos > 0).mean() >= 0.75)
    resumo = {
        "direcoes": linhas,
        "correlacao_mediana": float(np.median(correlacoes)),
        "coeficiente_mediano": float(np.median(coeficientes)),
        "ganho_mediano_sobre_o_controle": float(np.median(ganhos)),
        "direcoes_com_ganho_positivo": int((ganhos > 0).sum()),
        "direcoes_totais": len(ganhos),
        "magnitude_tipica_log2": float(np.median(normas)),
        "usar_resposta_media_transferida": usar,
        "regra": (
            "Usar a resposta media transferida exige ganho mediano positivo E ganho "
            "positivo em pelo menos 3 de cada 4 direcoes. Uma estimativa ruim da "
            "media piora os 300 alvos, enquanto a correcao certa so melhora 143 — a "
            "assimetria pede o teste mais exigente, nao o mais generoso."
        ),
        "consequencia": (
            "Se aprovado, somar coeficiente_mediano * media_do_doador_agrupada ao "
            "delta de TODOS os 300 alvos. Os 143 sem evidencia passam de abaixo de "
            "zero para zero; os 157 com evidencia ganham o mesmo deslocamento."
            if usar
            else
            "Se reprovado, nao somar nada: predizer o controle nos 143 e ruim, mas "
            "somar uma media mal estimada e pior, porque atinge os 300."
        ),
    }
    (SAIDA / "resposta_media.json").write_text(
        json.dumps(resumo, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\nveredito: {'USAR' if usar else 'NAO USAR'} a resposta media transferida")
    print(f"escrito em {SAIDA / 'resposta_media.json'}")


if __name__ == "__main__":
    main()
