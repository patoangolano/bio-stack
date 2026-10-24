"""
Confere o pacote controls-final contra a especificacao do VCC 2026.

Esperado: 18 533 genes, 300 perturbacoes, contextos D/E/F. Os .h5ad sao
abertos em modo backed: a matriz fica no disco, so obs/var sobem para a RAM.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

GENES_ESPERADOS = 18_533
PERTS_ESPERADAS = 300


def main() -> None:
    base = Path(sys.argv[1])

    manifesto = json.loads((base / "manifest.json").read_text(encoding="utf-8"))
    print("manifest.json")
    for chave, valor in manifesto.items():
        texto = json.dumps(valor, ensure_ascii=False)
        print(f"  {chave}: {texto[:140]}{'...' if len(texto) > 140 else ''}")

    genes = pd.read_csv(base / "gene_names.csv", header=None)
    marca = "OK" if len(genes) == GENES_ESPERADOS else f"!! esperado {GENES_ESPERADOS}"
    print(f"\ngene_names.csv: {len(genes)} linhas [{marca}]")

    perts = pd.read_csv(base / "pert_counts.csv")
    marca = "OK" if len(perts) == PERTS_ESPERADAS else f"!! esperado {PERTS_ESPERADAS}"
    print(f"pert_counts.csv: {len(perts)} linhas, colunas {list(perts.columns)} [{marca}]")
    print(f"  primeiras: {perts.iloc[:3, 0].tolist()}")

    for arquivo in sorted(base.glob("context_*.h5ad")):
        adata = ad.read_h5ad(arquivo, backed="r")
        print(f"\n{arquivo.name}: {adata.n_obs:,} celulas x {adata.n_vars:,} genes")
        print(f"  obs: {list(adata.obs.columns)}")
        for coluna in adata.obs.columns:
            serie = adata.obs[coluna]
            if serie.dtype.name in ("category", "object") and serie.nunique() <= 12:
                print(f"  {coluna}: {serie.value_counts().to_dict()}")
        adata.file.close()


if __name__ == "__main__":
    main()
