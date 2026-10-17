"""O sinal fraco de ODAPH fora da posicao, nas aves, e um locus ou e acaso?

No genoma inteiro (s35_busca_genoma), o melhor hit da proteina ODAPH do jacare
no frango (NC_052546.1:3.022.616, chr15) e no kakapo (NC_044279.2:31.599.253)
tem como melhor hit reciproco a propria ODAPH do jacare — mas a 32-33 bits,
encostado no piso de ruido (22-30). Reciprocidade de um trecho curto pode ser
acaso.

O teste: se os dois trechos estao em regioes SINTENICAS entre si (os mesmos
genes vizinhos no frango e no kakapo), o sinal e o mesmo locus em duas
linhagens que divergiram ha ~90 Ma — nao acaso. Compara os genes nomeados a
+-JANELA de cada hit e, como controle, a vizinhanca de um ponto aleatorio do
mesmo cromossomo do kakapo (gerador com semente fixa).

    uv run python scripts/atavismo/s35_sintenia_odaph.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_spec = importlib.util.spec_from_file_location("s27", Path(__file__).parent / "s27_painel_intervalos.py")
if _spec is None or _spec.loader is None:
    raise SystemExit("nao carreguei s27_painel_intervalos.py")
s27 = importlib.util.module_from_spec(_spec)
sys.modules["s27"] = s27
_spec.loader.exec_module(s27)

JANELA = 400_000
HITS = {"Ggal": ("NC_052546.1", 3_022_616), "Shab": ("NC_044279.2", 31_599_253)}
TAM_SHAB = 103_000_000  # NC_044279.2 tem ~103 Mb; o controle sorteia longe do hit


def vizinhos(acc: str, pos: int) -> set[str]:
    genes = s27.tabela_features(acc, max(1, pos - JANELA), pos + JANELA)
    return {g.sym.upper() for g in genes if not g.sym.startswith("LOC")}


def main() -> int:
    viz = {rot: vizinhos(acc, pos) for rot, (acc, pos) in HITS.items()}
    for rot, v in viz.items():
        acc, pos = HITS[rot]
        print(f"{rot} {acc}:{pos:,} +-{JANELA // 1000} kb: {len(v)} genes nomeados")
    comuns = viz["Ggal"] & viz["Shab"]
    print(f"\nem comum: {len(comuns)}  {', '.join(sorted(comuns)) or '-'}")

    rng = np.random.default_rng(20261017)
    print("\ncontrole — pontos aleatorios do mesmo cromossomo do kakapo:")
    for _ in range(3):
        p = int(rng.integers(JANELA, TAM_SHAB - JANELA))
        if abs(p - HITS["Shab"][1]) < 5 * JANELA:
            continue
        v = vizinhos(HITS["Shab"][0], p)
        print(f"  NC_044279.2:{p:,}  {len(v)} genes, em comum com o frango: {len(v & viz['Ggal'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
