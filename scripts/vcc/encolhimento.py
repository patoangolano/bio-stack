"""Fator de inflacao POR GENE para a parcela de desvio, medido contra verdade retida.

## O buraco que isto ataca

A varredura de `k_desvio` achou o otimo escalar em 4,0 e a tabela dos seis
membros explicou por que ele e um compromisso, e nao um otimo compartilhado:

    membro        k=1     k=2     k=3    [k=4]    k=5     k=6
    fidelity     0,062   0,452   0,781   0,834   0,823   0,802
    jaccard      0,084   0,401   0,532   0,554   0,551   0,541
    reach        0,343   0,620   0,715   0,756   0,779   0,790
    lfc_nmae     0,101   0,238   0,351   0,425   0,439   0,401
    expr_mse     0,232   0,385   0,447   0,395   0,179   0,000

`expr_mse` pica em 3, `lfc_nmae` em 5, `reach` sobe ate 6. Um escalar nao pode
servir os tres, porque um escalar tem um grau de liberdade e a discordancia tem
dois. Mas a discordancia nao e sobre o NIVEL: e sobre quais genes. Os membros de
direcao querem que o gene sub-movido ande; `expr_mse` quer que o gene bem
determinado NAO seja esticado. Esses sao conjuntos diferentes de genes, e um
fator por gene os separa.

## O que foi medido, e por que da direito de aplicar

Atenuacao por gene na direcao doador->alvo, contra a resposta REAL do alvo:

    a_g = sd_p(previsto[.,g]) / sd_p(real[.,g])        k_g = 1 / a_g

onde `previsto` e o preditor de producao (0,5*diagonal + 0,5*programas, ver
`transfer_fit.ajustar`). Medido nas direcoes com alvo profundo:

1. **A atenuacao GLOBAL cai em 5,1 a 6,2.** Isso corrobora o `k_desvio = 4` por
   um instrumento independente — o escalar veio de varredura no avaliador, este
   numero vem de verdade retida — e explica a diferenca sem apelar a erro:
   casar variancia pede 5 a 6, que e onde `reach` (6) e `lfc_nmae` (5) picam.
   `expr_mse` pede 3 porque erro quadratico e minimizado pela media a
   posteriori, que e ENCOLHIDA, nao calibrada. Os dois instrumentos concordam;
   o que difere e o funcional, nao a medida.

2. **O fator por gene varia 6 a 7 vezes** entre p05 e p95 (0,5x a 3,5x o
   global). Um escalar joga tudo isso fora.

3. **E transfere.** Spearman de `k_g` entre direcoes TOTALMENTE DISJUNTAS — sem
   doador em comum, sem alvo em comum — nos nove pares da sondagem: +0,425 a
   +0,638, mediana +0,56. Com doador em comum sobe a +0,835. Nao e ruido de par.

4. **Nao e artefato de cobertura.** `rho(k_g, n_valido) = -0,092`, contra
   `rho(k_g, var_prevista) = -0,806`. O que dirige o fator e quanto o preditor
   MOVE o gene: ele concentra movimento em poucos genes e sub-move o resto. Se
   o driver fosse a contagem de entradas validas, o fator seria sobre a
   qualidade da medida e nao sobre a predicao, e nao teria por que transferir.

## Por que so alvo profundo entra na medida

A variancia REAL do alvo e o denominador. Num alvo raso ela e dominada por
ruido de amostragem, o que a INFLA e deflaciona `k_g` — mediria a profundidade
do ensaio e nao a atenuacao do preditor. O gate e o mesmo de
`transfer_fit.politica_de_beta`, pela mesma razao: alvo em (replogle,
arc_vcc_25) com pelo menos `MIN_PERTS_ALVO` perturbacoes.

## Onde o peso e aplicado, e onde NAO e

Aqui so se mede e se grava a FORMA sobre os genes, com media geometrica 1. O
nivel continua sendo `k_desvio`, ja calibrado. Quem aplica e `varredura_k.py`:

    k_g = k_desvio * c(alpha) * w_g**alpha

`alpha` e o unico grau de liberdade novo, `alpha = 0` reproduz o braco escalar
bit a bit, e `c(alpha)` preserva a variancia TOTAL da parcela de desvio no
contexto — sem ele, `alpha` mexeria na forma e no nivel ao mesmo tempo, e a
varredura mediria a soma dos dois.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/encolhimento.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import transfer_fit as tf  # noqa: E402

VCC = Path("data/raw/vcc/controls-final")
SAIDA = Path("runs/vcc/transfer")

# Repos cuja resposta por perturbacao e profunda o bastante para servir de
# DENOMINADOR. Ver "Por que so alvo profundo entra na medida".
REPOS_ALVO_PROFUNDO = ("replogle", "arc_vcc_25")
# Perturbacoes em comum minimas para a direcao entrar. 400 e o mesmo corte de
# `politica_de_beta`; abaixo dele a variancia por gene do alvo tem erro padrao
# da ordem do proprio efeito que se quer medir.
MIN_PERTS_ALVO = 400
# Entradas validas minimas para o gene ganhar medida numa direcao. Com menos, a
# razao de duas variancias amostrais e ruido multiplicativo de dois digitos.
MIN_ENTRADAS_POR_GENE = 50
# Encolhimento do fator por gene para o global, em perturbacoes equivalentes.
# Mesmo idioma e mesmo valor de `transfer_fit.K_ENCOLHIMENTO`: com milhares de
# entradas acumuladas o encolhimento e leve, e so morde o gene medido numa
# direcao so e rasa — que e exatamente onde tem de morder.
K_ENCOLHIMENTO = 50.0
# Teto e piso da RAZAO CRUA `k_g / k_global` de uma direcao, em multiplos. O p95
# medido numa direcao e 3,5x e o p05 e 0,49x; 8x e 1/8 ficam folgados o bastante
# para nao recortar a cauda real e apertados o bastante para barrar um gene cujo
# preditor quase nao se move e cujo k_g estoura por divisao por quase zero.
#
# O teto APLICADO e menor que este, porque o recorte vem ANTES do encolhimento e
# os dois compoem: `TETO_RELATIVO ** (n/(n+K_ENCOLHIMENTO))`, que com n = 1.000
# da 7,25 e nao 8. A ordem e deliberada — o recorte limita a razao medida, o
# encolhimento limita a confianca nela, e inverter faria o recorte quase inerte
# porque o encolhimento ja teria domado o valor.
TETO_RELATIVO = 8.0


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman como Pearson dos postos — mesmo idioma de `transfer_fit.avaliar`.

    `stats.spearmanr` devolve um tipo que o verificador nao estreita, e
    `rankdata` ja resolve empate por media.
    """
    posto_a = np.asarray(stats.rankdata(a), dtype=np.float64)
    posto_b = np.asarray(stats.rankdata(b), dtype=np.float64)
    if np.ptp(posto_a) == 0 or np.ptp(posto_b) == 0:
        return float("nan")
    return float(np.corrcoef(posto_a, posto_b)[0, 1])


@dataclass(frozen=True)
class MedidaDeDirecao:
    """Atenuacao por gene numa direcao doador -> alvo."""

    rotulo: str
    genes: list[str]
    entradas: np.ndarray  # (n_genes,) entradas validas que entraram
    sq_previsto: np.ndarray  # (n_genes,) soma de quadrados do previsto
    sq_real: np.ndarray  # (n_genes,) soma de quadrados do real
    sq_previsto_diagonal: np.ndarray  # idem, so a parte diagonal
    k_global: float


def direcoes_com_alvo_profundo() -> list[tuple[str, str, str, str, list[str], list[str]]]:
    """Direcoes (doador, alvo) elegiveis, uma entrada por sentido.

    `pares_disponiveis` devolve pares nao ordenados; cada par da duas direcoes,
    e so a que tem alvo profundo entra.
    """
    saida: list[tuple[str, str, str, str, list[str], list[str]]] = []
    for repo_a, ctx_a, repo_b, ctx_b, perts, genes in tf.pares_disponiveis(tf.MINIMO_PERTS_PAR):
        if len(perts) < MIN_PERTS_ALVO:
            continue
        sentidos = (((repo_a, ctx_a), (repo_b, ctx_b)), ((repo_b, ctx_b), (repo_a, ctx_a)))
        for (rd, cd), (ra, ca) in sentidos:
            if ra in REPOS_ALVO_PROFUNDO:
                saida.append((rd, cd, ra, ca, perts, genes))
    return saida


def medir_direcao(
    repo_doador: str,
    ctx_doador: str,
    repo_alvo: str,
    ctx_alvo: str,
    perts: list[str],
    genes: list[str],
) -> MedidaDeDirecao:
    """Ajusta a transferencia nesta direcao e mede a atenuacao por gene.

    O previsto reproduz o preditor de producao: 0,5 diagonal + 0,5 programas.
    A parte diagonal pura viaja junto porque, em D/E/F, os genes fora da base de
    programas recebem SO a diagonal — e se as duas formas discordassem, aplicar
    uma medida com a outra seria erro silencioso. A discordancia e reportada.
    """
    doador = tf.construir(repo_doador, ctx_doador, perts, genes)
    alvo = tf.construir(repo_alvo, ctx_alvo, perts, genes)
    mascara = doador.valido & alvo.valido

    beta, _, _ = tf.coeficiente_por_gene(doador, alvo)
    base = tf.programas_de_resposta(doador, tf.K_PROGRAMAS)
    gama, carga = tf.coeficiente_por_programa(doador, alvo, base)
    diagonal = beta[None, :] * doador.log2r
    programas = (carga * gama[None, :]) @ base.T
    previsto = 0.5 * diagonal + 0.5 * programas

    entradas = mascara.sum(axis=0).astype(np.float64)
    sq_previsto = (np.square(previsto) * mascara).sum(axis=0)
    sq_real = (np.square(alvo.log2r) * mascara).sum(axis=0)
    sq_diagonal = (np.square(diagonal) * mascara).sum(axis=0)
    util = (entradas >= MIN_ENTRADAS_POR_GENE) & (sq_real > 0) & (sq_previsto > 0)
    k_global = float(np.sqrt(sq_real[util].sum() / sq_previsto[util].sum()))
    rotulo = f"{repo_doador}.{ctx_doador}->{repo_alvo}.{ctx_alvo}"
    print(f"  -- {rotulo}: {int(util.sum())} genes medidos, k global {k_global:.3f} --")
    return MedidaDeDirecao(
        rotulo=rotulo,
        genes=genes,
        entradas=np.where(util, entradas, 0.0),
        sq_previsto=np.where(util, sq_previsto, 0.0),
        sq_real=np.where(util, sq_real, 0.0),
        sq_previsto_diagonal=np.where(util, sq_diagonal, 0.0),
        k_global=k_global,
    )


def agrupar(
    medidas: list[MedidaDeDirecao], genes_painel: list[str]
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Forma por gene no eixo do painel, agrupando as direcoes em log.

    `k_g` e uma razao e o erro dela e multiplicativo, entao a media e em log. O
    peso e a contagem de entradas validas, e o encolhimento leva o gene pouco
    medido ao global DA PROPRIA DIRECAO — nao ao 1,0 absoluto — porque o global
    varia de 5,1 a 6,2 entre direcoes e misturar niveis contaminaria a forma.

    Devolve (peso com media geometrica 1, entradas acumuladas, diagnostico).
    """
    indice = {g: i for i, g in enumerate(genes_painel)}
    n = len(genes_painel)
    soma_log = np.zeros(n, dtype=np.float64)
    soma_peso = np.zeros(n, dtype=np.float64)
    globais: list[float] = []

    for medida in medidas:
        globais.append(medida.k_global)
        linhas = np.array([i for i, g in enumerate(medida.genes) if g in indice], dtype=np.int64)
        if linhas.size == 0:
            continue
        colunas = np.array([indice[medida.genes[i]] for i in linhas], dtype=np.int64)
        entradas = medida.entradas[linhas]
        util = entradas > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            k_gene = np.sqrt(medida.sq_real[linhas] / np.maximum(medida.sq_previsto[linhas], 1e-300))
        # Relativo ao global DESTA direcao, e so entao recortado: o recorte vale
        # sobre a forma, nao sobre o nivel.
        relativo = np.where(util, k_gene / medida.k_global, 1.0)
        relativo = np.clip(relativo, 1.0 / TETO_RELATIVO, TETO_RELATIVO)
        # Encolhimento para 1,0 (= o global da direcao) por contagem de entradas.
        fator = entradas / (entradas + K_ENCOLHIMENTO)
        log_relativo = np.where(util, np.log(relativo) * fator, 0.0)
        np.add.at(soma_log, colunas, np.where(util, log_relativo * entradas, 0.0))
        np.add.at(soma_peso, colunas, np.where(util, entradas, 0.0))

    medido = soma_peso > 0
    log_peso = np.zeros(n, dtype=np.float64)
    log_peso[medido] = soma_log[medido] / soma_peso[medido]
    # Media geometrica exatamente 1 sobre os genes medidos: a forma nao carrega
    # nivel. Gene sem medida fica em 1,0 e herda o escalar.
    log_peso[medido] -= log_peso[medido].mean()
    peso = np.exp(log_peso)

    diagnostico = {
        "k_global_mediano_entre_direcoes": float(np.median(globais)),
        "k_global_min": float(np.min(globais)),
        "k_global_max": float(np.max(globais)),
        "genes_medidos": float(medido.sum()),
        "genes_no_painel": float(n),
        "peso_p05": float(np.percentile(peso[medido], 5)),
        "peso_p50": float(np.percentile(peso[medido], 50)),
        "peso_p95": float(np.percentile(peso[medido], 95)),
    }
    return peso, soma_peso, diagnostico


def concordancia_entre_direcoes(medidas: list[MedidaDeDirecao]) -> dict[str, object]:
    """Spearman de `k_g` entre direcoes, separando pares disjuntos dos demais.

    Par disjunto — sem doador nem alvo em comum — e a unica evidencia de que a
    forma e propriedade do GENE. Num par que compartilha um lado, a correlacao
    pode vir do lado compartilhado.
    """
    por_direcao: dict[str, dict[str, float]] = {}
    for medida in medidas:
        util = medida.entradas > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            k = np.sqrt(medida.sq_real / np.maximum(medida.sq_previsto, 1e-300))
        por_direcao[medida.rotulo] = {
            g: float(v) for g, v, u in zip(medida.genes, k, util, strict=True) if u
        }

    disjuntos: list[float] = []
    compartilhados: list[float] = []
    detalhe: list[dict[str, object]] = []
    rotulos = sorted(por_direcao)
    for i, ra in enumerate(rotulos):
        for rb in rotulos[i + 1 :]:
            comuns = sorted(set(por_direcao[ra]) & set(por_direcao[rb]))
            if len(comuns) < 500:
                continue
            ka = np.array([por_direcao[ra][g] for g in comuns])
            kb = np.array([por_direcao[rb][g] for g in comuns])
            rho = spearman(ka, kb)
            lados_a = set(ra.split("->"))
            lados_b = set(rb.split("->"))
            disjunto = not (lados_a & lados_b)
            (disjuntos if disjunto else compartilhados).append(rho)
            detalhe.append(
                {"a": ra, "b": rb, "genes": len(comuns), "rho": rho, "disjunto": disjunto}
            )

    return {
        "pares_disjuntos": len(disjuntos),
        "rho_disjunto_mediano": float(np.median(disjuntos)) if disjuntos else float("nan"),
        "rho_disjunto_min": float(np.min(disjuntos)) if disjuntos else float("nan"),
        "rho_disjunto_max": float(np.max(disjuntos)) if disjuntos else float("nan"),
        "rho_compartilhado_mediano": float(np.median(compartilhados))
        if compartilhados
        else float("nan"),
        "por_par": detalhe,
    }


def discordancia_da_diagonal(medidas: list[MedidaDeDirecao]) -> dict[str, float]:
    """A forma medida com a mistura de producao serve para o gene so-diagonal?

    Em D/E/F os genes fora da base de programas recebem apenas a diagonal. Se a
    forma medida com 0,5/0,5 discordasse da medida com diagonal pura, aplicar
    uma no lugar da outra seria erro silencioso em milhares de genes.
    """
    rhos: list[float] = []
    for medida in medidas:
        util = (medida.entradas > 0) & (medida.sq_previsto_diagonal > 0)
        if util.sum() < 500:
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            mistura = np.sqrt(medida.sq_real[util] / medida.sq_previsto[util])
            diagonal = np.sqrt(medida.sq_real[util] / medida.sq_previsto_diagonal[util])
        rhos.append(spearman(mistura, diagonal))
    return {
        "direcoes": float(len(rhos)),
        "rho_mediano": float(np.median(rhos)) if rhos else float("nan"),
        "rho_min": float(np.min(rhos)) if rhos else float("nan"),
    }


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    genes_painel = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
    direcoes = direcoes_com_alvo_profundo()
    print(
        f"direcoes com alvo profundo ({'/'.join(REPOS_ALVO_PROFUNDO)}) e >= "
        f"{MIN_PERTS_ALVO} perts: {len(direcoes)}"
    )
    if not direcoes:
        raise ValueError("nenhuma direcao elegivel: sem alvo profundo nao ha denominador confiavel")

    medidas = [medir_direcao(*d) for d in direcoes]
    peso, entradas, diagnostico = agrupar(medidas, genes_painel)
    concordancia = concordancia_entre_direcoes(medidas)
    diagonal = discordancia_da_diagonal(medidas)

    print("\n== forma por gene, no eixo do painel ==")
    print(
        f"  {int(diagnostico['genes_medidos'])} de {int(diagnostico['genes_no_painel'])} "
        f"genes medidos; peso p05 {diagnostico['peso_p05']:.3f}  "
        f"p50 {diagnostico['peso_p50']:.3f}  p95 {diagnostico['peso_p95']:.3f}"
    )
    print(
        f"  k global entre direcoes: mediana {diagnostico['k_global_mediano_entre_direcoes']:.3f} "
        f"(min {diagnostico['k_global_min']:.3f}, max {diagnostico['k_global_max']:.3f})"
    )
    print(
        f"  concordancia em {concordancia['pares_disjuntos']} pares DISJUNTOS: rho mediano "
        f"{concordancia['rho_disjunto_mediano']:+.3f} "
        f"(min {concordancia['rho_disjunto_min']:+.3f}, "
        f"max {concordancia['rho_disjunto_max']:+.3f}); "
        f"pares com um lado em comum: {concordancia['rho_compartilhado_mediano']:+.3f}"
    )
    print(
        f"  forma de mistura vs forma so-diagonal: rho mediano {diagonal['rho_mediano']:+.3f} "
        f"em {int(diagonal['direcoes'])} direcoes (min {diagonal['rho_min']:+.3f})"
    )

    destino = SAIDA / "pesos_por_gene.npz"
    np.savez(
        destino,
        genes=np.array(genes_painel, dtype=str),
        peso=peso,
        entradas=entradas,
        k_global_medido=np.float64(diagnostico["k_global_mediano_entre_direcoes"]),
    )
    # O `.npz` nao entra no git — a politica do projeto mantem array binario de
    # `runs/vcc/` fora dele, e este tem 1,1 MB. Entao o que entra e o DIGEST, para
    # que "regere com encolhimento.py" seja conferivel e nao um pedido de fe: a
    # medida e determinista (a SVD randomizada usa `default_rng(0)`), logo uma
    # regeneracao correta reproduz o hash.
    digesto = hashlib.sha256(destino.read_bytes()).hexdigest()
    print(f"  sha256 do npz: {digesto}")
    (SAIDA / "encolhimento.json").write_text(
        json.dumps(
            {
                "direcoes": [m.rotulo for m in medidas],
                "gate_do_alvo": {
                    "repos": list(REPOS_ALVO_PROFUNDO),
                    "min_perts": MIN_PERTS_ALVO,
                    "por_que": (
                        "A variancia real do alvo e o denominador. Num alvo raso ela e "
                        "dominada por ruido de amostragem, o que a infla e deflaciona k_g: "
                        "mediria a profundidade do ensaio, nao a atenuacao do preditor."
                    ),
                },
                "min_entradas_por_gene": MIN_ENTRADAS_POR_GENE,
                "k_encolhimento": K_ENCOLHIMENTO,
                "teto_relativo": TETO_RELATIVO,
                "sha256_do_pesos_por_gene_npz": digesto,
                "diagnostico": diagnostico,
                "concordancia_entre_direcoes": concordancia,
                "forma_mistura_vs_diagonal": diagonal,
                "onde_o_peso_e_aplicado": (
                    "Aqui so se mede a FORMA, com media geometrica 1. O nivel continua "
                    "sendo k_desvio. varredura_k.py aplica k_g = k_desvio * c(alpha) * "
                    "peso**alpha, com c(alpha) preservando a variancia total da parcela "
                    "de desvio no contexto; alpha=0 reproduz o braco escalar bit a bit."
                ),
                "leitura": (
                    "rho disjunto e a evidencia que importa: sem doador nem alvo em comum, "
                    "uma correlacao positiva so pode vir do gene. Se ela caisse a zero, a "
                    "forma seria ruido de par e o escalar ficaria como esta."
                ),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(f"\nescrito em {destino} e {SAIDA / 'encolhimento.json'}")


if __name__ == "__main__":
    main()
