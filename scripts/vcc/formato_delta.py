"""O arquivo de delta: a interface entre o preditor, o gerador e o validador numerico.

Um `.npz` com:

    genes       texto [n_genes]      eixo de genes, na ordem do gene_names.csv
    alvos       texto [n_alvos]      alvos, na ordem do pert_counts.csv
    contextos   texto [n_contextos]  ex.: D, E, F
    delta_<c>   float32 [n_alvos, n_genes], um por contexto

`delta_<c>[i, g]` e o deslocamento PRETENDIDO do pseudobulk do alvo i no
contexto c, no espaco `bulk_lognorm` (log1p CP50k, ver pseudobulk.py), relativo
ao pseudobulk dos controles do mesmo contexto no controls-final.

Mesmas chaves do `delta_previsto.npz` antigo, com uma diferenca: os rotulos sao
texto de largura fixa, nao arrays de objeto, entao o arquivo carrega com
`allow_pickle=False`. Pickle desserializado executa codigo; um arquivo de dados
nao tem por que exigir isso.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np


@dataclass(frozen=True)
class Delta:
    genes: list[str]
    alvos: list[str]
    contextos: list[str]
    por_contexto: dict[str, np.ndarray]

    def __post_init__(self) -> None:
        if sorted(self.por_contexto) != sorted(self.contextos):
            raise ValueError(f"contextos {self.contextos} vs chaves {sorted(self.por_contexto)}")
        forma = (len(self.alvos), len(self.genes))
        for c, d in self.por_contexto.items():
            if d.shape != forma:
                raise ValueError(f"delta_{c} tem forma {d.shape}, esperado {forma}")


def salvar_delta(caminho: Path, delta: Delta) -> None:
    arrays: dict[str, np.ndarray] = {
        "genes": np.array(delta.genes, dtype=str),
        "alvos": np.array(delta.alvos, dtype=str),
        "contextos": np.array(delta.contextos, dtype=str),
    }
    arrays.update({f"delta_{c}": np.asarray(d, dtype=np.float32) for c, d in delta.por_contexto.items()})
    # cast: sem ele o verificador supoe que um array pode cair no parametro allow_pickle.
    np.savez(caminho, **cast("dict[str, Any]", arrays))


def ler_delta(caminho: Path) -> Delta:
    try:
        with np.load(caminho, allow_pickle=False) as z:
            contextos = [str(c) for c in z["contextos"]]
            return Delta(
                genes=[str(g) for g in z["genes"]],
                alvos=[str(a) for a in z["alvos"]],
                contextos=contextos,
                por_contexto={c: np.asarray(z[f"delta_{c}"]) for c in contextos},
            )
    except ValueError as erro:
        if "allow_pickle" in str(erro):
            raise ValueError(
                f"{caminho}: os rotulos estao em arrays de objeto e so carregariam com pickle; "
                "regrave com formato_delta.salvar_delta"
            ) from erro
        raise
