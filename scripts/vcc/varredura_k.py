"""
Separa o delta previsto em suas duas parcelas e gera combinacoes escaladas, para
calibrar o encolhimento contra as SEIS metricas do vcc2026 — e nao contra o MSE,
que foi quem o escolheu.

## O buraco que isto ataca, agora medido

Pontuacao da v4 no `avaliador_proxy`, recorte de 40 alvos reais em D, contra o
teto do oraculo:

    pds_cosine          0,990 de 1,000    99% do alcancavel
    direction_reach     0,267 de 0,989    27%
    expr_mse            0,228 de 0,947    24%
    lfc_nmae            0,102 de 1,024    10%
    sig_jaccard         0,040 de 1,032     4%
    direction_fidelity  0,022 de 0,828     3%   <- o buraco

O modelo acerta o padrao geral e erra a DIRECAO por perturbacao. Direcao por
perturbacao e, por construcao, a parcela do **desvio especifico** — que o beta
de minimos quadrados encolhe a ~0,22 porque minimos quadrados e otimo para erro
quadratico e sub-disperso para tudo que depende de acertar o sinal.

## O resultado, para nao reler a hipotese como se fosse pergunta aberta

A varredura rodou e respondeu: `k_desvio = 4`, e `avg_score` foi de 0,275 para
0,647 contra 0,970 do oraculo. A hipotese abaixo estava certa no diagnostico e
ERRADA no mecanismo — nao ha troca ate k=3, tudo sobe junto, inclusive o
`expr_mse` que minimos quadrados teoricamente otimiza (0,228 -> 0,442). O beta
nao trocava direcao por magnitude: estava sub-disperso para os seis membros ao
mesmo tempo. Detalhe em docs/VCC-2026-EXECUCAO.md secao 21.

Duas armadilhas de leitura da tabela de ESCALAS_DESVIO:

  - `expr_mse` = 0,0000 em k=6 e k=8 e CORTE de escala, nao empate com o
    controle. O bruto vai a 1,337 e 3,092, PIOR que o baseline. O `avg_score`
    desses dois esta otimista e a queda acima de 5 e penhasco, nao rampa.
  - k=1 reproduz `delta_previsto.npz` a 7e-8 relativo. Se nao reproduzisse, a
    varredura mediria o erro da decomposicao e nao `k`.

## Duas parcelas, dois escalares

O delta de cada alvo tem duas partes:

    resposta media do contexto   vale para os 300, identica entre eles
    desvio especifico            so nos 157 com evidencia

A decomposicao e EXATA e nao precisa recalcular o modelo: os 143 alvos sem
evidencia carregam **apenas** a media, e sao identicos entre si (verificado).
Entao o delta de qualquer um deles E a parcela media.

## Por que a decomposicao acontece em log2r, e nao em delta

O delta e funcao nao linear da razao — `log1p` de uma proporcao renormalizada.
Duas parcelas NAO se somam em espaco de delta. Elas se somam em espaco de
log-razao, que e onde o modelo as construiu. Entao:

    log2r_total  = log2(expm1(b_ctrl + delta) / expm1(b_ctrl))
    log2r_media  = o mesmo, a partir do delta de um alvo sem evidencia
    log2r_desvio = log2r_total - log2r_media
    log2r_novo   = k_media * log2r_media + k_desvio * log2r_desvio

e so entao converte de volta, com renormalizacao composicional.

## O que este proxy pode e nao pode calibrar

**Pode `k_desvio`.** A verdade do proxy e Orion/HCT116 transplantado, que contem
o desvio especifico — e que portanto mede quando ele esta fraco demais.

**Nao pode `k_media`.** A mesma verdade NAO contem a resposta media de D. O
proxy e cego a essa parcela, provado: a correcao do vies de selecao mudou a
magnitude dela de 0,063 para 0,019 e moveu `avg_score` em +0,0028, dentro do
ruido. Varrer `k_media` aqui mediria nada.

Por isso o padrao e `k_media = 1` e a varredura e em `k_desvio`.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/varredura_k.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# `-I` tira o diretorio do script do sys.path; o mesmo idioma do pool_media.py.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from formato_delta import Delta, salvar_delta  # noqa: E402

TRANSFER = Path("runs/vcc/transfer")
EIXO = Path("runs/vcc/eixo")
CTRL = Path("runs/vcc/ctrl_bulk")
SAIDA = Path("runs/vcc/transfer/varredura_k")
CONTEXTOS = ("D", "E", "F")

# Escalas do DESVIO ESPECIFICO. 1,0 e o beta de minimos quadrados como esta.
#
# MEDIDO (avaliador_proxy, bundle proxy-D-orion-k40, 40 alvos reais, 11 bracos):
#   k      1     1,5    2      3     [4]     5      6      8    oraculo
#   avg  0,275  0,370  0,460  0,590  0,647  0,627  0,587  0,526  0,970
# Curva unimodal, otimo em k=4, plato largo de 3 a 5. O grid original parava em
# 4 subindo — 5, 6 e 8 foram acrescentados justamente para pos o otimo entre
# parenteses, porque "melhor dos testados" e mais fraco que "o otimo esta aqui".
ESCALAS_DESVIO = (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
# A media fica em 1,0: o proxy e cego a ela (ver docstring).
ESCALA_MEDIA = 1.0
# Genes significativos por perturbacao no dado real, para leitura local.
DE_REAIS = 2739
LIMIAR_DE = 0.010


def para_log2r(delta: np.ndarray, proporcao: np.ndarray) -> np.ndarray:
    """Razao de proporcao, em log2, implicada por um delta sobre `b_ctrl`."""
    alcancavel = proporcao > 0
    base = np.where(alcancavel, proporcao, 1.0)
    b_ctrl = np.log1p(proporcao)
    with np.errstate(divide="ignore", invalid="ignore"):
        razao = np.expm1(b_ctrl[None, :] + delta) / base[None, :]
    valido = alcancavel[None, :] & np.isfinite(razao) & (razao > 0)
    return np.where(valido, np.log2(np.where(valido, razao, 1.0)), 0.0)


def para_delta(log2r: np.ndarray, proporcao: np.ndarray) -> np.ndarray:
    """Caminho inverso, com renormalizacao composicional por linha."""
    alcancavel = proporcao > 0
    b_ctrl = np.log1p(proporcao)
    perturbada = proporcao[None, :] * np.exp2(log2r)
    massa_alvo = float(proporcao[alcancavel].sum())
    massa = perturbada[:, alcancavel].sum(axis=1, keepdims=True)
    perturbada = perturbada * np.where(massa > 0, massa_alvo / np.maximum(massa, 1e-30), 1.0)
    novo = np.log1p(perturbada) - b_ctrl[None, :]
    novo[:, ~alcancavel] = 0.0
    return novo.astype(np.float32)


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    pacote = np.load(TRANSFER / "delta_previsto.npz", allow_pickle=True)
    genes = np.array([str(g) for g in pacote["genes"]], dtype=object)
    alvos = np.array([str(a) for a in pacote["alvos"]], dtype=object)
    por_gene = pd.read_csv(CTRL / "por_gene.tsv", sep="\t")
    eixo = json.loads((EIXO / "resumo.json").read_text(encoding="utf-8"))
    sem_evidencia = set(eixo["evidencia_doadora"]["alvos_sem_evidencia"])
    posicoes_sem = [i for i, a in enumerate(alvos) if a in sem_evidencia]
    if not posicoes_sem:
        raise ValueError("nenhum alvo sem evidencia: a parcela media nao e isolavel")

    linhas: list[dict[str, float | str]] = []
    decomposicao: dict[str, dict[str, float]] = {}
    base_log2r: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for contexto in CONTEXTOS:
        proporcao = np.expm1(np.asarray(por_gene[f"b_ctrl_{contexto}"], dtype=np.float64))
        delta = np.asarray(pacote[f"delta_{contexto}"], dtype=np.float64)
        log2r_total = para_log2r(delta, proporcao)

        # A parcela media: o delta de um alvo sem evidencia carrega so ela.
        candidatos = log2r_total[posicoes_sem]
        if not np.allclose(candidatos, candidatos[0], atol=1e-6):
            raise ValueError(
                f"contexto {contexto}: alvos sem evidencia nao sao identicos; "
                "a decomposicao assume que carregam apenas a resposta media"
            )
        log2r_media = candidatos[0]
        log2r_desvio = log2r_total - log2r_media[None, :]
        base_log2r[contexto] = (proporcao, log2r_media, log2r_desvio)

        norma_media = float(np.linalg.norm(log2r_media))
        normas_desvio = np.linalg.norm(log2r_desvio, axis=1)
        com_desvio = int((normas_desvio > 1e-9).sum())
        decomposicao[contexto] = {
            "norma_da_media": norma_media,
            "norma_mediana_do_desvio": float(np.median(normas_desvio[normas_desvio > 1e-9]))
            if com_desvio
            else 0.0,
            "alvos_com_desvio": float(com_desvio),
        }
        print(
            f"{contexto}: |media| = {norma_media:.3f}, "
            f"|desvio| mediano = {decomposicao[contexto]['norma_mediana_do_desvio']:.3f}, "
            f"{com_desvio} alvos com desvio"
        )

    for k in ESCALAS_DESVIO:
        saidas: dict[str, np.ndarray] = {}
        for contexto in CONTEXTOS:
            proporcao, log2r_media, log2r_desvio = base_log2r[contexto]
            combinado = ESCALA_MEDIA * log2r_media[None, :] + k * log2r_desvio
            novo = para_delta(combinado, proporcao)
            saidas[f"delta_{contexto}"] = novo
            universo = np.asarray(por_gene[f"no_universo_{contexto}"], dtype=bool)
            no_universo = np.abs(novo[:, universo])
            acima = float(np.median((no_universo > LIMIAR_DE).sum(axis=1)))
            linhas.append(
                {
                    "k_desvio": k,
                    "contexto": contexto,
                    "genes_acima_do_limiar": acima,
                    "razao_para_de_real": acima / DE_REAIS,
                }
            )
        destino = SAIDA / f"delta_kd{k:g}.npz"
        # Pelo `formato_delta`, que e o que o `gerar_submissao` le: rotulos em
        # texto de largura fixa, carrega com `allow_pickle=False`. A versao de
        # 22/10 gravava arrays de objeto, que so carregam com pickle.
        salvar_delta(
            destino,
            Delta(
                genes=[str(g) for g in genes],
                alvos=[str(a) for a in alvos],
                contextos=list(CONTEXTOS),
                por_contexto={c: saidas[f"delta_{c}"] for c in CONTEXTOS},
            ),
        )
        medio = float(np.mean([float(x["genes_acima_do_limiar"]) for x in linhas if x["k_desvio"] == k]))
        print(
            f"k_desvio={k:4g}  genes acima de {LIMIAR_DE}: {medio:7.0f} "
            f"({medio / DE_REAIS:5.2f}x os {DE_REAIS} reais)  -> {destino.name}"
        )

    pd.DataFrame(linhas).to_csv(SAIDA / "contagem_por_k.tsv", sep="\t", index=False)
    (SAIDA / "varredura.json").write_text(
        json.dumps(
            {
                "escalas_do_desvio": list(ESCALAS_DESVIO),
                "escala_da_media": ESCALA_MEDIA,
                "decomposicao": decomposicao,
                "onde_o_k_se_aplica": "em log2 da razao de proporcao, nao no delta",
                "por_que_a_media_nao_e_varrida": (
                    "A verdade do avaliador_proxy e Orion/HCT116 transplantado, que contem "
                    "o desvio especifico e NAO contem a resposta media de D. Provado: a "
                    "correcao do vies de selecao mudou a magnitude da media de 0,063 para "
                    "0,019 e moveu avg_score em +0,0028, dentro do ruido. Varrer k_media "
                    "neste instrumento mediria nada."
                ),
                "alvo_da_varredura": (
                    "direction_fidelity, em 0,022 de um teto de 0,828 — 3% do alcancavel. "
                    "Direcao por perturbacao e a parcela do desvio especifico, que o beta "
                    "de minimos quadrados encolhe a ~0,22."
                ),
                "aviso": (
                    "A contagem por limiar e aproximacao declarada. So o avaliador decide: "
                    "emitir cada k e pontuar, registrando os SEIS membros separados."
                ),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(f"\n{len(ESCALAS_DESVIO)} arquivos em {SAIDA}")


if __name__ == "__main__":
    main()
