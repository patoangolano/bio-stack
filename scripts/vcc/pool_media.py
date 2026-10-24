"""
Qual esquema de PONDERACAO agrupa melhor a resposta media dos contextos doadores.

## O que motivou, e e um numero, nao uma suspeita

A resposta media do contexto e a parcela maior da predicao: vale para os 300
alvos, contra 157 do desvio especifico, e `resposta_media.json` mede ganho de
+0,252 sobre predizer o controle em 24 de 24 direcoes. No `transfer_fit.ajustar`
ela e agrupada entre contextos doadores ponderada por `profundidade[repo]`.

Dois defeitos medidos nesse agrupamento:

**1. O peso esta invertido contra a medida.** Peso efetivo (profundidade vezes
numero de contextos usaveis) contra transferencia medida da media:

    contexto               peso    share   ganho medido como doador
    arc_vcc_25.ARC_H1      2016    93,2%   +0,0000  (r = -0,0001, 8.800 genes)
    jiang (11 contextos)    143     6,6%   nao medido
    orion.hct116              2     0,1%   +0,3917  (6 direcoes)
    orion.hek293t             2     0,1%   +0,6589  (1 direcao)
    replogle (4 contextos)    0     0,0%   +0,195 a +0,434

O contexto que leva 93% do peso tem media ORTOGONAL a do Orion, medida em 8.800
genes com magnitudes comparaveis (0,034 contra 0,043) — nao e falta de
sobreposicao nem artefato do filtro de cobertura.

**2. Um gate serve duas parcelas com requisitos diferentes.** O laco do
`ajustar` comeca com

    usaveis = set(alvos) & disponiveis & perts_treino[dataset]
    if not usaveis: continue

e a contribuicao da MEDIA esta dentro desse `continue`. O gate e correto para o
desvio, que precisa do nosso alvo medido naquele contexto. A media nao precisa:
`media_do_contexto()` percorre TODAS as perturbacoes do doador. Consequencia
medida: os 4 contextos do Replogle, que tem a melhor transferencia de media do
conjunto, entram com peso ZERO porque nenhum dos nossos 300 alvos esta no split
de treino do Replogle.

## O RESULTADO, para o leitor nao reler a hipotese como pergunta aberta

Corrigido o GATE, nao o peso. Deixar-um-de-fora com `--sem-mesmo-repo`, que e o
analogo honesto da nossa situacao (D/E/F nao tem doador do proprio estudo):

    profundidade   mediana +0,0983   media +0,0859
    uniforme       mediana +0,0579   media +0,0873
    atual          mediana +0,0524   media +0,0618
    concordancia   mediana +0,0403   media +0,0463
    amostra        mediana +0,0231   media +0,0448

Mediana diz profundidade, media diz peso igual, e as medias diferem em 0,0014.
Sinal misto da ordem do ruido nao escolhe vencedor, entao o PESO fica como
estava. O que muda e so o gate, que ganha em 4 de 4 configuracoes.

**Sem `--sem-mesmo-repo` a conclusao se inverte, e o motivo e vazamento.** Nessa
forma `uniforme` vencia com +0,2176 contra +0,1660 de profundidade — e eu ia
implementar peso igual. Mas nas duas linhas do ensaio da PROPRIA competicao
profundidade ganhava 5x (+0,5248 e +0,5503 contra +0,1038 e +0,1713), porque com
`ARC_H1` retido ela punha 91% do peso em `ARC_H1_VAL`: o mesmo estudo. Retirar o
contexto nao basta; tem de retirar o ESTUDO, senao mede-se atalho em vez de
generalizacao. Com o atalho fora, todos os ganhos caem por 2 a 3.

**`concordancia` perdeu nas quatro configuracoes** (cortes de 10, 40 e 200, e sem
mesmo estudo). Era o lever 4 da revisao de literatura e a hipotese que motivou
este script. A medida descarta.

## O teste, e por que ele fecha sem D/E/F

O `avaliador_proxy` e CEGO a parcela da media (a verdade dele e Orion/HCT116
transplantado, que nao contem a media de D) — entao trocar o peso nao e
mensuravel la. Mas e mensuravel aqui, em dado real: deixa-se um contexto de
FORA, agrupa-se a media dos outros com cada esquema, e ve-se qual prediz melhor
a media do contexto retido.

Metrica identica a do `resposta_media.py`, para os numeros serem comparaveis:
coeficiente sem intercepto, e ganho `1 - SSE(coef*agrupada, retida)/SSE(0, retida)`.
Negativo significa que agrupar e pior que predizer o controle.

**Sem vazamento.** O peso de concordancia de um doador e recalculado por
exclusao: para o contexto retido `h`, o peso de `d` e a mediana do ganho de `d`
nas direcoes `d -> y` com `y != h`. Usar a medida que inclui `h` como alvo seria
treinar no conjunto de teste.

Rodar da raiz do projeto, depois de `transfer_fit.py validar`:
    ./.venv/Scripts/python.exe -I scripts/vcc/pool_media.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import transfer_fit as tf  # noqa: E402

SAIDA = Path("runs/vcc/transfer")
# Profundidade por repo (genes significativos/pert, mediana) — o peso atual.
PROFUNDIDADE = {"arc_vcc_25": 2016.0, "replogle": 46.0, "jiang": 13.0, "orion": 2.0}
# Um contexto entra na comparacao se mediu pelo menos isto de perturbacoes: com
# menos, a media dele e ruido e ele polui o teste sem informar nada.
MIN_PERTS_CONTEXTO = 40
# Genes em comum exigidos para uma comparacao valer.
MIN_GENES_COMPARAVEIS = 500


def ganho_sem_intercepto(previsto: np.ndarray, real: np.ndarray) -> tuple[float, float, float]:
    """(coeficiente, ganho sobre predizer zero, correlacao), como no resposta_media."""
    denominador = float(previsto @ previsto)
    if denominador <= 1e-12:
        return 0.0, 0.0, 0.0
    coeficiente = float(previsto @ real) / denominador
    sse = float(np.square(coeficiente * previsto - real).sum())
    sse_zero = float(np.square(real).sum())
    ganho = 1.0 - sse / max(sse_zero, 1e-12)
    correlacao = float(np.corrcoef(previsto, real)[0, 1]) if previsto.size > 2 else 0.0
    return coeficiente, ganho, correlacao


def concordancia_por_doador(excluir_alvo: str | None) -> dict[str, float]:
    """Ganho mediano de cada contexto COMO DOADOR, ignorando direcoes para `excluir_alvo`."""
    caminho = SAIDA / "resposta_media.json"
    if not caminho.exists():
        raise FileNotFoundError("rode resposta_media.py antes: sem resposta_media.json")
    direcoes = json.loads(caminho.read_text(encoding="utf-8"))["direcoes"]
    acumulado: dict[str, list[float]] = {}
    for registro in direcoes:
        nome = str(registro["direcao"])
        doador, _, alvo = nome.partition("_para_")
        if excluir_alvo is not None and alvo == excluir_alvo:
            continue
        acumulado.setdefault(doador, []).append(float(registro["ganho_sobre_o_controle"]))
    return {d: float(np.median(v)) for d, v in acumulado.items()}


def main() -> None:
    analisador = argparse.ArgumentParser(description=__doc__ or "")
    analisador.add_argument(
        "--min-perts",
        type=int,
        default=MIN_PERTS_CONTEXTO,
        help=(
            "minimo de perturbacoes para um contexto entrar. O default cai UMA acima "
            "do maior contexto do Jiang (39), o que e perto demais para nao ser testado"
        ),
    )
    analisador.add_argument(
        "--sem-mesmo-repo",
        action="store_true",
        help=(
            "ao retirar um contexto, retira TAMBEM todos os outros do mesmo repo. "
            "E o analogo honesto da nossa situacao: D/E/F nao tem doador do proprio "
            "estudo, e sem isso o esquema que concentra peso num contexto irmao ganha "
            "por atalho de mesmo-estudo em vez de por generalizacao"
        ),
    )
    analisador.add_argument("--rotulo", default=None, help="sufixo do arquivo de saida")
    argumentos = analisador.parse_args()
    min_perts = int(argumentos.min_perts)
    genes_painel = pd.read_csv(tf.VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
    indice_gene = {g: i for i, g in enumerate(genes_painel)}
    alvos = set(pd.read_csv(tf.VCC / "pert_counts.csv")["target_gene"].astype(str))
    split = json.loads(
        (tf.METAS / "PIE_splits" / "replogle_xdataset" / "train.json").read_text(encoding="utf-8")
    )
    perts_treino: dict[str, set[str]] = {}
    for chave, lista in split.items():
        perts_treino.setdefault(str(chave).split(".", 1)[0], set()).update(str(p) for p in lista)

    # Media de cada contexto doador no NOSSO eixo, com mascara de validade.
    medias: dict[str, np.ndarray] = {}
    valido: dict[str, np.ndarray] = {}
    meta_ctx: dict[str, dict[str, float]] = {}
    print("== media por contexto doador, no eixo do painel ==")
    for repo, pasta in tf.REPOS.items():
        meta = json.loads((tf.caminho_meta(pasta)).read_text(encoding="utf-8"))
        dataset = str(meta["dataset"])
        catalogo = tf.perts_por_contexto(repo)
        genes_repo = [g for g in tf.genes_do_repo(repo) if g in indice_gene]
        colunas = np.array([indice_gene[g] for g in genes_repo], dtype=np.int64)
        for contexto, disponiveis in sorted(catalogo.items()):
            if len(disponiveis) < min_perts:
                continue
            rotulo = f"{repo}_{contexto}".replace(".", "_")
            bruta = tf.media_do_contexto(repo, contexto, genes_repo)
            vetor = np.zeros(len(genes_painel), dtype=np.float64)
            marca = np.zeros(len(genes_painel), dtype=bool)
            vetor[colunas] = bruta
            marca[colunas] = bruta != 0.0
            medias[rotulo] = vetor
            valido[rotulo] = marca
            usaveis = len(alvos & disponiveis & perts_treino.get(dataset, set()))
            meta_ctx[rotulo] = {
                "repo_id": float(list(tf.REPOS).index(repo)),
                "perts": float(len(disponiveis)),
                "profundidade": PROFUNDIDADE[repo],
                "alvos_usaveis": float(usaveis),
                "genes_com_estimativa": float(marca.sum()),
            }
            print(
                f"  {rotulo:30s} {len(disponiveis):5d} perts | "
                f"{int(marca.sum()):5d} genes com estimativa | "
                f"{usaveis:3d} dos nossos alvos no split | "
                f"|media| {np.linalg.norm(bruta) / max(np.sqrt(bruta.size), 1):.5f}"
            )

    rotulos = sorted(medias)
    if len(rotulos) < 3:
        raise ValueError("menos de 3 contextos: deixar-um-de-fora nao diz nada")

    esquemas = ("atual", "profundidade", "concordancia", "uniforme", "amostra")
    resultados: dict[str, list[dict[str, float | str]]] = {e: [] for e in esquemas}

    print(f"\n== deixar-um-de-fora sobre {len(rotulos)} contextos ==")
    cabecalho = "  ".join(f"{e:>13s}" for e in esquemas)
    print(f"{'contexto retido':30s} {cabecalho}")
    for retido in rotulos:
        concordancia = concordancia_por_doador(retido)
        proibidos = {retido}
        if argumentos.sem_mesmo_repo:
            for outro in rotulos:
                if outro != retido and meta_ctx[outro]["repo_id"] == meta_ctx[retido]["repo_id"]:
                    proibidos.add(outro)
        pesos_por_esquema: dict[str, dict[str, float]] = {}
        for esquema in esquemas:
            pesos: dict[str, float] = {}
            for doador in rotulos:
                if doador in proibidos:
                    continue
                info = meta_ctx[doador]
                if esquema == "atual":
                    # Reproduz o agrupamento de hoje: profundidade, e o gate do
                    # desvio aplicado tambem a media.
                    pesos[doador] = info["profundidade"] if info["alvos_usaveis"] > 0 else 0.0
                elif esquema == "profundidade":
                    pesos[doador] = info["profundidade"]
                elif esquema == "concordancia":
                    # Sem medida nao ha peso. Nao e punicao: e a assimetria de uma
                    # media mal estimada atingir os 300 alvos.
                    pesos[doador] = max(0.0, concordancia.get(doador, 0.0))
                elif esquema == "uniforme":
                    pesos[doador] = 1.0
                else:
                    pesos[doador] = info["perts"]
            pesos_por_esquema[esquema] = pesos

        real = medias[retido]
        marca_real = valido[retido]
        linha: list[str] = []
        for esquema in esquemas:
            pesos = pesos_por_esquema[esquema]
            numerador = np.zeros(len(genes_painel), dtype=np.float64)
            denominador = np.zeros(len(genes_painel), dtype=np.float64)
            for doador, peso in pesos.items():
                if peso <= 0.0:
                    continue
                marca = valido[doador]
                numerador[marca] += peso * medias[doador][marca]
                denominador[marca] += peso
            agrupada = np.where(denominador > 0, numerador / np.maximum(denominador, 1e-12), 0.0)
            comparavel = marca_real & (denominador > 0)
            if int(comparavel.sum()) < MIN_GENES_COMPARAVEIS:
                linha.append(f"{'—':>13s}")
                resultados[esquema].append({"retido": retido, "ganho": float("nan")})
                continue
            coeficiente, ganho, correlacao = ganho_sem_intercepto(
                agrupada[comparavel], real[comparavel]
            )
            linha.append(f"{ganho:+13.4f}")
            resultados[esquema].append(
                {
                    "retido": retido,
                    "genes": float(comparavel.sum()),
                    "coeficiente": coeficiente,
                    "correlacao": correlacao,
                    "ganho": ganho,
                    "doadores_com_peso": float(sum(1 for p in pesos.values() if p > 0)),
                }
            )
        print(f"{retido:30s} {'  '.join(linha)}")

    print("\n== resumo por esquema ==")
    placar: dict[str, dict[str, float]] = {}
    for esquema in esquemas:
        ganhos = np.array(
            [float(x["ganho"]) for x in resultados[esquema] if np.isfinite(float(x["ganho"]))]
        )
        if not ganhos.size:
            continue
        placar[esquema] = {
            "ganho_mediano": float(np.median(ganhos)),
            "ganho_medio": float(ganhos.mean()),
            "minimo": float(ganhos.min()),
            "maximo": float(ganhos.max()),
            "contextos_com_ganho_positivo": float((ganhos > 0).sum()),
            "contextos": float(ganhos.size),
        }
    for esquema, v in sorted(placar.items(), key=lambda kv: -kv[1]["ganho_mediano"]):
        print(
            f"  {esquema:14s} mediana {v['ganho_mediano']:+.4f}  media {v['ganho_medio']:+.4f}  "
            f"faixa [{v['minimo']:+.4f}, {v['maximo']:+.4f}]  "
            f"positivos {int(v['contextos_com_ganho_positivo'])}/{int(v['contextos'])}"
        )

    vencedor = max(placar, key=lambda e: placar[e]["ganho_mediano"])
    atual = placar.get("atual", {}).get("ganho_mediano", float("nan"))
    print(
        f"\nvencedor: {vencedor} (mediana {placar[vencedor]['ganho_mediano']:+.4f} "
        f"contra {atual:+.4f} do esquema atual)"
    )
    sufixo = f"_{argumentos.rotulo}" if argumentos.rotulo else ""
    (SAIDA / f"pool_media{sufixo}.json").write_text(
        json.dumps(
            {
                "min_perts_contexto": min_perts,
                "sem_mesmo_repo": bool(argumentos.sem_mesmo_repo),
                "contextos": meta_ctx,
                "por_esquema": resultados,
                "placar": placar,
                "vencedor": vencedor,
                "o_que_cada_esquema_e": {
                    "atual": "profundidade do repo, com o gate do desvio aplicado a media",
                    "profundidade": "profundidade do repo, sem o gate",
                    "concordancia": "ganho mediano medido do doador, recalculado por exclusao",
                    "uniforme": "peso igual para todo contexto",
                    "amostra": "numero de perturbacoes do contexto",
                },
                "limite_honesto": (
                    "Mede agrupamento entre contextos que EXISTEM, nao para D/E/F. O que "
                    "transfere e a ordem entre esquemas, nao o nivel do ganho. E o teste "
                    "e circular num ponto declarado: o esquema `concordancia` usa a mesma "
                    "familia de medida (ganho como doador) que a metrica do teste, so que "
                    "em pares que nao incluem o contexto retido. A exclusao remove o "
                    "vazamento direto; nao remove a afinidade de forma entre as duas."
                ),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(f"escrito em {SAIDA / f'pool_media{sufixo}.json'}")


if __name__ == "__main__":
    main()
