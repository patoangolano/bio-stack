#!/usr/bin/env python
"""Baixa os arrays PIE do Hugging Face que o pipeline de VCC 2026 consome.

Subconjunto deliberado, nao o repositorio inteiro. Conferido contra os 4,47 GB
registrados da execucao anterior:

  * Somente `preprocessed/`. Os `de/*.parquet` ficam de fora porque a variavel
    transferida e `log2(fold_changes)` mascarada por `tested` -- razao
    adimensional aplicada ao nosso proprio `b_ctrl`. Em PIE_x_atlas_orion
    sozinho os parquets de DE somam 5,2 GB que nada le.
  * `lfc_true.npy` EXCLUIDO: 2,36 GB so no orion, e e o desfecho do PIE, nao
    entrada nossa.
  * PIE_tahoe100m EXCLUIDO inteiro.
  * PIE_splits baixa completo (1,7 MB) -- `replogle_xdataset/train.json` e o
    que define os 157 alvos com evidencia doadora, contra os 161 que aparecem
    em algum `pert_to_id`.
  * PIE_sources: SO os metadados (1,27 MB de um repositorio de 44 GB). O
    `eixo_18262.py` le a lista de chaves (`meta.json`) e os aliases de cinco
    fontes para saber quais alvos cada uma cobre; os embeddings nao entram em
    nada daqui -- os 143 alvos sem evidencia deram R2 negativo com eles.

Os quatro repositorios sao publicos (HTTP 200 sem token, medido 2026-10-23).

Uso:
    uv run python scripts/vcc/baixar_pie.py
    uv run python scripts/vcc/baixar_pie.py --conferir   # so lista, nao baixa
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from huggingface_hub import HfApi, snapshot_download

RAIZ = Path(__file__).resolve().parents[2]
DESTINO = RAIZ / "data" / "raw" / "pie" / "arrays"

# Os quatro doadores. A profundidade medida (mediana de genes com fdr<0,05 por
# perturbacao) esta no comentario porque e ela, nao a sobreposicao de alvo, que
# decide qual doador presta: a referencia real do desafio tem ~2.739.
DOADORES = {
    "PIE_x_atlas_orion": 2,  # 148 dos 157 alvos, e mediana de DOIS genes
    "PIE_arc_vcc_25": 2016,  # o unico no regime certo -- e o ensaio da propria competicao
    "PIE_jiang": 13,
    "PIE_replogle_nadig_essential": 46,
}

PADROES_ACEITOS = ["preprocessed/*"]
PADROES_RECUSADOS = ["*lfc_true*"]

# As cinco fontes cujas chaves o `eixo_18262.py` consulta.
FONTES = ("gene_text", "esm2", "depmap_gene_effect", "string_space", "perturbation_text")
PADROES_FONTES = [f"{f}/{arq}" for f in FONTES for arq in ("meta.json", "aliases.yaml")]


def tamanho_previsto(api: HfApi) -> int:
    """Soma o que seria baixado, para conferir contra os 4,47 GB registrados."""
    total = 0
    for nome in [*DOADORES, "PIE_splits"]:
        # dataset_info, nao repo_info: o segundo devolve uniao com ModelInfo/SpaceInfo
        # e o verificador de tipo nao consegue estreitar para o ramo com `siblings`.
        info = api.dataset_info(f"arcinstitute/{nome}", files_metadata=True)
        for irmao in info.siblings or []:
            caminho = irmao.rfilename
            if nome != "PIE_splits":
                if not caminho.startswith("preprocessed/"):
                    continue
                if "lfc_true" in caminho:
                    continue
            tamanho = irmao.size or 0
            total += tamanho
            print(f"  {tamanho / 1e6:9.1f} MB  {nome}/{caminho}")
    return total


def baixar() -> None:
    DESTINO.mkdir(parents=True, exist_ok=True)
    for nome, profundidade in DOADORES.items():
        print(f"--- {nome} (profundidade mediana {profundidade})", flush=True)
        snapshot_download(
            repo_id=f"arcinstitute/{nome}",
            repo_type="dataset",
            local_dir=DESTINO / nome,
            allow_patterns=PADROES_ACEITOS,
            ignore_patterns=PADROES_RECUSADOS,
        )
    print("--- PIE_splits (completo)", flush=True)
    snapshot_download(
        repo_id="arcinstitute/PIE_splits",
        repo_type="dataset",
        local_dir=DESTINO.parent / "PIE_splits",
    )
    baixar_fontes()


def baixar_fontes() -> None:
    print("--- PIE_sources (so metadados)", flush=True)
    snapshot_download(
        repo_id="arcinstitute/PIE_sources",
        repo_type="dataset",
        local_dir=DESTINO.parent / "PIE_sources",
        allow_patterns=PADROES_FONTES,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--conferir",
        action="store_true",
        help="lista o que seria baixado e soma o tamanho, sem transferir nada",
    )
    ap.add_argument(
        "--so-fontes",
        action="store_true",
        help="baixa so os metadados do PIE_sources, sem tocar nos 4,2 GB de arrays",
    )
    args = ap.parse_args()
    if args.so_fontes:
        baixar_fontes()
        presentes = sorted((DESTINO.parent / "PIE_sources").rglob("*.*"))
        print(f"{len([p for p in presentes if '.cache' not in p.parts])} arquivos de metadados")
        return 0

    api = HfApi()
    if args.conferir:
        total = tamanho_previsto(api)
        print(f"\ntotal previsto: {total / 1e9:.2f} GB (registrado da execucao anterior: 4,47 GB)")
        return 0

    baixar()
    presentes = sorted(DESTINO.rglob("*.npy")) + sorted(DESTINO.rglob("*.json"))
    no_disco = sum(p.stat().st_size for p in presentes)
    print(f"\n{len(presentes)} arquivos, {no_disco / 1e9:.2f} GB em {DESTINO}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
