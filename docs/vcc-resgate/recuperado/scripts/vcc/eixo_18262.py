"""
Etapa 1 de docs/VCC-2026-ROTA.md: eixo canonico de genes, cobertura das fontes
de conhecimento, mapas de indice para cada dataset PIE e gate de expressao dos
300 alvos.

Tres cruzamentos, que sao coisas diferentes e vivem em eixos diferentes:

1. EIXO DE SAIDA (gene). Os 18.533 genes do painel VCC contra as 31.684 chaves
   de `gene_text`. Os ausentes nao saem da submissao — a submissao tem de
   devolver os 18.533. Perdem apenas a feature de texto, que e imputada pela
   media da coluna. O `datamodule` da Arc nao aplica alias a `gene_text`, por
   isso aqui tambem nao se aplica: o numero de cobertura tem de ser o mesmo que
   o PIE veria.

2. EIXO DE PERTURBACAO (alvo). Os 300 alvos contra `esm2`, `depmap_gene_effect`,
   `string_space` e `perturbation_text`. Esses quatro tem `aliases.yaml` e o
   alias conta, porque e assim que a fonte e consultada.

3. EIXO CRUZADO (gene e alvo em cada dataset PIE). Para cada repo preprocessado,
   o indice do nosso gene no eixo dele e o indice do nosso alvo no `pert_to_id`.
   Sem isso a etapa 3 nao sabe de que coluna ler a evidencia. Orion mistura
   simbolo HGNC e ENSG no eixo de genes; so o casamento por simbolo vale.

O gate de expressao usa a saida de caracterizar_contextos.py. CRISPRi reprime
transcricao: alvo que o contexto nao transcreve nao tem o que derrubar, e a
predicao ali deve tender ao controle em vez de inventar efeito.

Saidas em runs/vcc/eixo/. Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/eixo_18262.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

VCC = Path("data/raw/vcc/controls-final")
FONTES = Path("data/raw/pie/PIE_sources")
DATASETS = Path("data/raw/pie")
CARACT = Path("runs/vcc/caracterizacao")
SAIDA = Path("runs/vcc/eixo")

# Indexadas por perturbacao: features dos 300 alvos.
FONTES_PERT = ("esm2", "depmap_gene_effect", "string_space", "perturbation_text")
# Repos preprocessados com arrays de resposta. Tahoe fica fora (droga-dose).
REPOS_PIE = (
    "PIE_x_atlas_orion",
    "PIE_replogle_nadig_essential",
    "PIE_arc_vcc_25",
    "PIE_jiang",
)
# Abaixo disso o alvo praticamente nao e transcrito no contexto: nada a reprimir.
FRACAO_GATE = 0.02


def carregar_chaves(nome: str) -> tuple[set[str], dict[str, str]]:
    """Chaves de uma fonte PIE_sources e seu mapa de alias (antigo -> novo)."""
    meta = json.loads((FONTES / nome / "meta.json").read_text(encoding="utf-8"))
    chaves = {str(k) for k in meta["keys"]}
    caminho_alias = FONTES / nome / "aliases.yaml"
    alias: dict[str, str] = {}
    if caminho_alias.exists():
        cru = yaml.safe_load(caminho_alias.read_text(encoding="utf-8")) or {}
        alias = {str(k): str(v) for k, v in cru.items()}
    return chaves, alias


def resolver(simbolo: str, chaves: set[str], alias: dict[str, str]) -> str | None:
    """Nome com que a fonte conhece este simbolo, ou None se nao conhece."""
    if simbolo in chaves:
        return simbolo
    novo = alias.get(simbolo)
    if novo is not None and novo in chaves:
        return novo
    return None


def classificar_ausente(simbolo: str) -> str:
    """Por que este simbolo provavelmente nao esta em gene_text.

    Rotulo descritivo para decidir imputacao com os olhos abertos; nao e usado
    em nenhum calculo.
    """
    if simbolo.startswith(("HIST1", "HIST2", "HIST3", "H1-", "H2A", "H2B", "H3-", "H4-")):
        return "histona (renomeada no HGNC 2019-2021)"
    if simbolo.startswith(("H2BF", "H3F", "H1F")):
        return "histona (simbolo antigo)"
    if "orf" in simbolo:
        return "orf nao caracterizada"
    if simbolo.startswith(("LINC", "MIR", "SNOR", "RNU", "SCARNA")):
        return "nao codificante"
    if "-" in simbolo and simbolo.split("-")[-1].isdigit():
        return "nome de clone / sufixo numerico"
    if simbolo.endswith(("P1", "P2", "P3", "P4")):
        return "pseudogene"
    return "outro"


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
    alvos = pd.read_csv(VCC / "pert_counts.csv")["target_gene"].astype(str).tolist()
    print(f"painel VCC: {len(genes)} genes, {len(alvos)} alvos")

    # --- 1. eixo de saida contra gene_text (sem alias, de proposito) ---
    chaves_gene_text, _ = carregar_chaves("gene_text")
    coberto = [g in chaves_gene_text for g in genes]
    n_coberto = sum(coberto)
    ausentes = [g for g, ok in zip(genes, coberto, strict=True) if not ok]
    print(f"gene_text: {n_coberto} de {len(genes)} cobertos, {len(ausentes)} ausentes")

    tabela_ausentes = pd.DataFrame(
        {
            "gene": ausentes,
            "idx_vcc": [genes.index(g) for g in ausentes],
            "motivo_provavel": [classificar_ausente(g) for g in ausentes],
        }
    )
    tabela_ausentes.to_csv(SAIDA / "ausentes_gene_text.tsv", sep="\t", index=False)
    contagem_motivos = tabela_ausentes["motivo_provavel"].value_counts().to_dict()

    # --- 2. eixo de perturbacao: as quatro fontes contra os 300 alvos ---
    cobertura_alvos: dict[str, dict[str, str | None]] = {}
    for fonte in FONTES_PERT:
        chaves, alias = carregar_chaves(fonte)
        cobertura_alvos[fonte] = {a: resolver(a, chaves, alias) for a in alvos}
        n = sum(1 for v in cobertura_alvos[fonte].values() if v is not None)
        n_por_alias = sum(
            1 for a, v in cobertura_alvos[fonte].items() if v is not None and v != a
        )
        print(f"{fonte:20s}: {n:3d} de {len(alvos)} alvos" + (f" ({n_por_alias} por alias)" if n_por_alias else ""))

    # --- 3. eixo cruzado: indices em cada dataset PIE ---
    idx_gene_por_repo: dict[str, dict[str, int]] = {}
    idx_pert_por_repo: dict[str, dict[str, int]] = {}
    resumo_repos: list[dict[str, str | int]] = []
    for repo in REPOS_PIE:
        meta = json.loads((DATASETS / repo / "meta.json").read_text(encoding="utf-8"))
        eixo_repo = {str(g): i for i, g in enumerate(meta["genes"])}
        perts_repo = {str(p): int(i) for p, i in meta["pert_to_id"].items()}
        idx_gene_por_repo[repo] = {g: eixo_repo[g] for g in genes if g in eixo_repo}
        idx_pert_por_repo[repo] = {a: perts_repo[a] for a in alvos if a in perts_repo}
        resumo_repos.append(
            {
                "repo": repo,
                "dataset": str(meta["dataset"]),
                "genes_no_repo": int(meta["num_genes"]),
                "genes_vcc_mapeados": len(idx_gene_por_repo[repo]),
                "contextos": int(meta["num_contexts"]),
                "perts_no_repo": int(meta["num_perts"]),
                "alvos_vcc_mapeados": len(idx_pert_por_repo[repo]),
                "linhas": int(meta["num_rows"]),
            }
        )
        print(
            f"{repo:30s}: {len(idx_gene_por_repo[repo]):5d}/{len(genes)} genes, "
            f"{len(idx_pert_por_repo[repo]):3d}/{len(alvos)} alvos, {meta['num_rows']} linhas"
        )

    # --- eixo canonico, uma linha por gene do painel ---
    eixo = pd.DataFrame({"gene": genes})
    eixo["idx_vcc"] = range(len(genes))
    eixo["em_gene_text"] = coberto
    for repo in REPOS_PIE:
        eixo[f"idx_{repo.removeprefix('PIE_')}"] = [
            idx_gene_por_repo[repo].get(g, -1) for g in genes
        ]
    eixo["n_repos_com_gene"] = sum(
        (eixo[f"idx_{r.removeprefix('PIE_')}"] >= 0).astype(int) for r in REPOS_PIE
    )
    eixo.to_csv(SAIDA / "eixo_canonico.tsv", sep="\t", index=False)
    eixo.loc[eixo["em_gene_text"], ["gene", "idx_vcc"]].to_csv(
        SAIDA / "eixo_18262.tsv", sep="\t", index=False
    )

    # --- gate de expressao dos 300 alvos, por contexto ---
    expressao = pd.read_csv(CARACT / "alvos_por_contexto.tsv", sep="\t")
    expressao["gate"] = [
        "sem_alvo_transcrito" if f < FRACAO_GATE else "ativo"
        for f in expressao["fracao_expressa"]
    ]
    for fonte in FONTES_PERT:
        expressao[f"em_{fonte}"] = [
            cobertura_alvos[fonte][a] is not None for a in expressao["alvo"]
        ]
    for repo in REPOS_PIE:
        expressao[f"evid_{repo.removeprefix('PIE_')}"] = [
            a in idx_pert_por_repo[repo] for a in expressao["alvo"]
        ]
    expressao.to_csv(SAIDA / "gate_alvos.tsv", sep="\t", index=False)
    gate_por_contexto = (
        expressao.groupby(["contexto", "gate"]).size().unstack(fill_value=0).to_dict()
    )
    print(f"\ngate de expressao (fracao < {FRACAO_GATE}): {gate_por_contexto}")

    # --- evidencia doadora: duas contagens, e elas divergem de proposito ---
    #
    # A contagem CRUA olha `pert_to_id` do repo: tudo que foi medido, inclusive
    # linhas que a Arc reservou para validacao e teste. A contagem de TREINO
    # olha PIE_splits/replogle_xdataset/train.json. A segunda e a que vale para
    # a etapa 3: treinar numa linha de validacao e vazamento, e e exatamente a
    # diferenca entre as duas (5 alvos do ARC_H1) que fazia a conta dar 161.
    com_evidencia_crua = sorted({a for repo in REPOS_PIE for a in idx_pert_por_repo[repo]})
    split = json.loads(
        (DATASETS / "PIE_splits" / "replogle_xdataset" / "train.json").read_text(encoding="utf-8")
    )
    perts_treino_por_dataset: dict[str, set[str]] = {}
    for chave, lista in split.items():
        dataset = str(chave).split(".", 1)[0]
        perts_treino_por_dataset.setdefault(dataset, set()).update(str(p) for p in lista)

    conjunto_alvos = set(alvos)
    treino_por_repo: dict[str, list[str]] = {}
    for item in resumo_repos:
        repo = str(item["repo"])
        dataset = str(item["dataset"])
        no_treino = sorted(conjunto_alvos & perts_treino_por_dataset.get(dataset, set()))
        treino_por_repo[repo] = no_treino
        item["alvos_vcc_no_split_treino"] = len(no_treino)

    com_evidencia_treino = sorted({a for lista in treino_por_repo.values() for a in lista})
    sem_evidencia = sorted(conjunto_alvos - set(com_evidencia_treino))
    reservados = sorted(set(com_evidencia_crua) - set(com_evidencia_treino))
    print(
        f"evidencia doadora: {len(com_evidencia_crua)} medidos, "
        f"{len(com_evidencia_treino)} no split de treino, {len(sem_evidencia)} sem nenhuma"
    )
    detalhe = ", ".join(f"{r.removeprefix('PIE_')}={len(v)}" for r, v in treino_por_repo.items())
    print(f"  por repo no treino: {detalhe}")
    print(f"  medidos mas reservados para validacao/teste ({len(reservados)}): {', '.join(reservados)}")

    resumo = {
        "painel": {"genes": len(genes), "alvos": len(alvos)},
        "gene_text": {
            "chaves_na_fonte": len(chaves_gene_text),
            "cobertos": n_coberto,
            "ausentes": len(ausentes),
            "motivos": contagem_motivos,
            "decisao_imputacao": (
                "Os ausentes PERMANECEM no eixo de saida: a submissao exige os 18.533. "
                "Perdem so a feature de texto, imputada pela media da coluna de treino. "
                "Nao sao gateados a zero — expressao de controle existe para todos, e a "
                "evidencia cruzada de outro dataset nao depende de gene_text. Alias nao "
                "e aplicado aqui porque o datamodule da Arc tambem nao aplica."
            ),
        },
        "fontes_por_alvo": {
            fonte: {
                "cobertos": sum(1 for v in cobertura_alvos[fonte].values() if v is not None),
                "por_alias": sum(
                    1 for a, v in cobertura_alvos[fonte].items() if v is not None and v != a
                ),
            }
            for fonte in FONTES_PERT
        },
        "datasets_pie": resumo_repos,
        "evidencia_doadora": {
            "medidos_em_algum_repo": len(com_evidencia_crua),
            "no_split_de_treino": len(com_evidencia_treino),
            "por_repo_no_treino": {r: len(v) for r, v in treino_por_repo.items()},
            "reservados_para_validacao": reservados,
            "sem_nenhuma": len(sem_evidencia),
            "alvos_sem_evidencia": sem_evidencia,
            "qual_usar": (
                "no_split_de_treino. Treinar numa linha que a Arc reservou para "
                "validacao e vazamento; os reservados entram so como checagem."
            ),
        },
        "gate": {"fracao_minima": FRACAO_GATE, "por_contexto": gate_por_contexto},
    }
    (SAIDA / "resumo.json").write_text(
        json.dumps(resumo, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    print(f"\nescrito em {SAIDA}")


if __name__ == "__main__":
    main()
