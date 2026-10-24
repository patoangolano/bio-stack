"""
Os 143 alvos do VCC 2026 sem nenhuma evidencia doadora recebem predicao zero,
que pontua exatamente zero na escala do desafio. Este script pergunta se e
possivel tira-los do piso usando as tres fontes de conhecimento que cobrem
300 de 300 alvos: `esm2`, `depmap_gene_effect` e `string_space`.

## A pergunta, e a regra de decisao escrita antes do resultado

Treinar nos 157 alvos COM evidencia um preditor de "carga de programa de
resposta" a partir das features do alvo, validar fora da amostra, e aplicar
aos 143 so se funcionar.

**Regra, fixada antes de rodar:** os 143 entram na submissao apenas se o
cosseno fora da amostra entre delta previsto e delta observado for positivo
com margem clara sobre zero em validacao cruzada. Cosseno negativo ou proximo
de zero significa que a predicao e pior ou igual a nao predizer nada, e nesse
caso zero continua sendo a resposta.

## Por que o prior e pessimista

O artigo do PIE mede o ganho TOTAL das cinco fontes de conhecimento em regime
de contexto e perturbacao ambos novos: 0,010 de AUPRC. Praticamente nada. E o
nosso sinal de treino vem em 148 dos 157 casos do Orion, que tem mediana de
dois genes significativos por perturbacao. Um mapa perfeito de features para
cargas ainda mostraria correlacao baixa, porque o rotulo e ruidoso.

Medir de todo jeito e barato, e a resposta negativa tambem e resultado: diz
que os 143 ficam no piso por falta de sinal, nao por falta de tentativa.

## Metodo

1. PCA do delta ja previsto para os 157 (saida de `transfer_fit.py ajustar`),
   retendo os componentes dominantes. As cargas sao o alvo da regressao.
2. Features: concatenacao das tres fontes, reduzida por PCA. A reducao e
   ajustada nos 300 alvos — e nao supervisionada, nao olha resposta nenhuma,
   entao nao vaza.
3. Ridge com validacao cruzada em k dobras sobre os 157. Ridge e nao minimos
   quadrados porque sao 2.978 features originais para 157 exemplos.
4. Relatorio por componente e reconstrucao fora da amostra.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/alvos_sem_evidencia.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

FONTES = Path("data/raw/pie/arrays/PIE_sources")
TRANSFER = Path("runs/vcc/transfer")
VCC = Path("data/raw/vcc/controls-final")
SAIDA = Path("runs/vcc/sem_evidencia")

NOMES_FONTES = ("esm2", "depmap_gene_effect", "string_space")
# Componentes retidos do delta. 24 e menos que os 32 programas do transfer_fit
# porque aqui ha 157 exemplos, nao 1.000 perturbacoes.
K_COMPONENTES = 24
# Dimensao das features depois da reducao. 48 features para 157 exemplos ja e
# uma razao apertada; mais que isso e ajuste de ruido garantido.
K_FEATURES = 48
# Dobras da validacao cruzada.
DOBRAS = 5
# Forcas de regularizacao varridas. A escolhida e a de melhor cosseno fora da
# amostra, e a varredura inteira vai para o relatorio para que se veja se o
# otimo e um pico estreito (suspeito) ou um platô (confiavel).
ALFAS = (1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0)
# Margem minima de cosseno fora da amostra para os 143 entrarem.
MARGEM_MINIMA = 0.05


def carregar_features(alvos: list[str]) -> tuple[np.ndarray, dict[str, int]]:
    """Concatena as tres fontes por alvo, padronizando cada fonte separadamente.

    Cada fonte tem escala propria (esm2 e embedding de proteina, depmap e
    efeito de essencialidade, string e espaco de rede). Padronizar por fonte
    antes de concatenar evita que a de maior variancia domine a PCA.
    """
    blocos: list[np.ndarray] = []
    dimensoes: dict[str, int] = {}
    for nome in NOMES_FONTES:
        meta = json.loads((FONTES / nome / "meta.json").read_text(encoding="utf-8"))
        chaves = {str(k): i for i, k in enumerate(meta["keys"])}
        alias_path = FONTES / nome / "aliases.yaml"
        alias: dict[str, str] = {}
        if alias_path.exists():
            for linha in alias_path.read_text(encoding="utf-8").splitlines():
                if ":" in linha:
                    esquerda, direita = linha.split(":", 1)
                    alias[esquerda.strip()] = direita.strip()
        matriz = np.load(FONTES / nome / "embeddings.npy", allow_pickle=False)
        bloco = np.zeros((len(alvos), matriz.shape[1]), dtype=np.float64)
        faltando = 0
        for posicao, alvo in enumerate(alvos):
            indice = chaves.get(alvo)
            if indice is None:
                indice = chaves.get(alias.get(alvo, ""))
            if indice is None:
                faltando += 1
                continue
            bloco[posicao] = matriz[indice]
        media = bloco.mean(axis=0, keepdims=True)
        desvio = bloco.std(axis=0, keepdims=True)
        bloco = (bloco - media) / np.where(desvio > 0, desvio, 1.0)
        blocos.append(bloco)
        dimensoes[nome] = int(matriz.shape[1])
        print(f"  {nome:20s} {matriz.shape[1]:5d} dims | alvos sem entrada: {faltando}")
    return np.hstack(blocos), dimensoes


def reduzir(matriz: np.ndarray, k: int) -> np.ndarray:
    """Projeta nas k direcoes de maior variancia, por SVD economica."""
    centrada = matriz - matriz.mean(axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centrada, full_matrices=False)
    return centrada @ vt[:k].T


def ridge(x: np.ndarray, y: np.ndarray, alfa: float) -> np.ndarray:
    """Coeficientes de ridge com intercepto implicito (entradas centradas)."""
    gram = x.T @ x + alfa * np.eye(x.shape[1])
    return np.linalg.solve(gram, x.T @ y)


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    alvos = pd.read_csv(VCC / "pert_counts.csv")["target_gene"].astype(str).tolist()
    pacote = np.load(TRANSFER / "delta_previsto.npz", allow_pickle=True)

    # A fonte autoritativa de quem tem evidencia doadora, gravada pela etapa 1.
    eixo = json.loads((Path("runs/vcc/eixo") / "resumo.json").read_text(encoding="utf-8"))
    sem_evidencia = set(eixo["evidencia_doadora"]["alvos_sem_evidencia"])
    print(f"alvos sem evidencia doadora, pelo eixo: {len(sem_evidencia)} de {len(alvos)}")

    print("== features ==")
    features, dimensoes = carregar_features(alvos)
    print(f"  concatenado: {features.shape[1]} dims -> reduzido para {K_FEATURES}")
    features_reduzidas = reduzir(features, K_FEATURES)

    relatorio: dict[str, object] = {
        "k_componentes": K_COMPONENTES,
        "k_features": K_FEATURES,
        "dobras": DOBRAS,
        "dimensoes_das_fontes": dimensoes,
        "margem_minima_de_cosseno": MARGEM_MINIMA,
        "regra_de_decisao": (
            "Os 143 alvos sem evidencia entram na submissao apenas se o cosseno fora "
            "da amostra superar a margem minima. Regra fixada antes de rodar."
        ),
        "por_contexto": {},
    }
    veredito: dict[str, bool] = {}

    for contexto in ("D", "E", "F"):
        delta = np.asarray(pacote[f"delta_{contexto}"], dtype=np.float64)
        # Quem tem evidencia vem do eixo, NAO de "delta diferente de zero".
        #
        # O criterio antigo (`abs(delta).sum() > 0`) funcionava antes de a
        # resposta media do contexto passar a valer para os 300 alvos. Depois
        # disso ele devolve 300, e 144 deles compartilham uma linha IDENTICA —
        # que entraria em treino e teste de todas as 5 dobras e faria a
        # validacao cruzada medir a capacidade de prever uma constante.
        #
        # E o mesmo padrao da secao 16 da execucao: a camada de baixo infere uma
        # propriedade que a de cima apagou de proposito, em vez de ler a fonte.
        indice_alvo = {a: i for i, a in enumerate(alvos)}
        com_sinal = np.array(
            sorted(indice_alvo[a] for a in alvos if a not in sem_evidencia),
            dtype=np.int64,
        )
        if com_sinal.size < 50:
            raise ValueError(f"contexto {contexto}: so {com_sinal.size} alvos com sinal")
        print(f"\n== contexto {contexto}: {com_sinal.size} alvos com evidencia ==")

        # Base de componentes do delta observado, ajustada so nos que tem sinal.
        observado = delta[com_sinal]
        centro = observado.mean(axis=0, keepdims=True)
        centrado = observado - centro
        _, _, vt = np.linalg.svd(centrado, full_matrices=False)
        base = vt[:K_COMPONENTES].T
        cargas = centrado @ base
        variancia_retida = float(
            np.square(cargas).sum() / max(np.square(centrado).sum(), 1e-12)
        )
        print(f"  {K_COMPONENTES} componentes retem {variancia_retida:.3f} da variancia do delta")

        x = features_reduzidas[com_sinal]
        indices = np.arange(com_sinal.size)
        dobra_de = indices % DOBRAS

        melhor = (-np.inf, 0.0)
        varredura: list[dict[str, float]] = []
        for alfa in ALFAS:
            previsto = np.zeros_like(cargas)
            for dobra in range(DOBRAS):
                treino = dobra_de != dobra
                teste = ~treino
                media_x = x[treino].mean(axis=0, keepdims=True)
                media_y = cargas[treino].mean(axis=0, keepdims=True)
                coef = ridge(x[treino] - media_x, cargas[treino] - media_y, alfa)
                previsto[teste] = (x[teste] - media_x) @ coef + media_y
            # Cosseno no espaco do delta, que e onde a metrica vive.
            delta_previsto = previsto @ base.T
            delta_observado = centrado
            normas = np.linalg.norm(delta_previsto, axis=1) * np.linalg.norm(
                delta_observado, axis=1
            )
            validos = normas > 0
            cosseno = float(
                np.mean(
                    (delta_previsto[validos] * delta_observado[validos]).sum(axis=1)
                    / normas[validos]
                )
            )
            r2 = float(
                1.0
                - np.square(previsto - cargas).sum() / max(np.square(cargas).sum(), 1e-12)
            )
            varredura.append({"alfa": alfa, "cosseno_fora_da_amostra": cosseno, "r2": r2})
            print(f"    alfa={alfa:9.0f} cosseno={cosseno:+.4f} r2={r2:+.4f}")
            if cosseno > melhor[0]:
                melhor = (cosseno, alfa)

        cosseno_melhor, alfa_melhor = melhor
        passa = cosseno_melhor > MARGEM_MINIMA
        veredito[contexto] = passa
        print(
            f"  melhor: alfa={alfa_melhor:.0f} cosseno={cosseno_melhor:+.4f} -> "
            f"{'USAR nos 143' if passa else 'NAO usar; os 143 ficam em zero'}"
        )
        assert isinstance(relatorio["por_contexto"], dict)
        relatorio["por_contexto"][contexto] = {
            "alvos_com_evidencia": int(com_sinal.size),
            "variancia_retida_pelos_componentes": variancia_retida,
            "varredura_de_alfa": varredura,
            "melhor_alfa": alfa_melhor,
            "melhor_cosseno_fora_da_amostra": cosseno_melhor,
            "aprovado": passa,
        }

    relatorio["veredito"] = veredito
    relatorio["conclusao"] = (
        "Os 143 alvos sem evidencia entram na submissao nos contextos aprovados. "
        "Onde o cosseno fora da amostra nao supera a margem, a predicao fica em zero, "
        "que pontua zero — e zero e melhor que uma predicao pior que zero."
        if any(veredito.values())
        else (
            "Nenhum contexto aprovou. Os 143 alvos sem evidencia ficam com predicao "
            "zero. As tres fontes de conhecimento cobrem 300/300 alvos, mas nao "
            "carregam sinal suficiente para prever a carga de resposta quando o "
            "rotulo de treino vem de um ensaio com mediana de dois genes "
            "significativos por perturbacao. Resultado negativo, e medido."
        )
    )
    (SAIDA / "relatorio.json").write_text(
        json.dumps(relatorio, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\n{relatorio['conclusao']}")
    print(f"escrito em {SAIDA / 'relatorio.json'}")


if __name__ == "__main__":
    main()
