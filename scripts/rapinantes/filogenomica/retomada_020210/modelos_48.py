# modelos_48.py — rebuild the 48-group partition model the 020210 supermatrix ran under, from its own ModelFinder cache.
# The scheme comes from concatenado.best_scheme.nex (written by ModelFinder on 2026-10-18 14:25); each group's model is
# the BIC-best entry of concatenado.model.gz, exactly the criterion ModelFinder applied. The totals must reproduce the
# "Best partition model BIC score: 4467486.907 (LnL: -2230458.528  df:535)" line of the original log, or nothing is written.
import gzip
import math
import re
import sys
from pathlib import Path

LNL_ORIGINAL = -2230458.528
# df of the whole partition model = sum of the per-group df (each carries its own rate scaler) - 1 fixed scaler
# + 111 branch lengths of the shared unrooted tree (2 * 57 taxa - 3), with -p (edge-proportional).
DF_ORIGINAL = 535
MODELO = re.compile(r"^ ([A-Z][A-Za-z0-9.]*(?:\+[A-Za-z0-9.]+)*): (-[0-9.]+) ([0-9]+) [0-9.e+-]+$")


def ler_cache(caminho: Path) -> dict[str, dict[str, tuple[float, int]]]:
    """Key -> {model: (lnL, df)} from a ModelFinder .model.gz; parameter, tree and best_* lines are skipped."""
    cache: dict[str, dict[str, tuple[float, int]]] = {}
    chave = None
    with gzip.open(caminho, "rt") as fh:
        for linha in fh:
            if linha.startswith("#"):
                continue
            if not linha.startswith(" ") and linha.rstrip().endswith(":"):
                chave = linha.rstrip()[:-1]
                cache[chave] = {}
            elif chave is not None and (m := MODELO.match(linha)):
                cache[chave][m[1]] = (float(m[2]), int(m[3]))
    return cache


def main(td: Path, saida: Path) -> None:
    cache = ler_cache(td / "concatenado.model.gz")
    grupos = re.findall(r"charset (\S+) = ([^;]+);", (td / "concatenado.best_scheme.nex").read_text())
    if len(grupos) != 48:
        sys.exit(f"best_scheme.nex com {len(grupos)} grupos, esperados 48")

    linhas_charset, linhas_modelo = [], []
    soma_lnl, soma_df, n_total = 0.0, 0, 0
    for nome, faixas in grupos:
        n_sitios = sum(
            int(fim or ini) - int(ini) + 1 for ini, _, fim in (f.partition("-") for f in faixas.split())
        )
        # In memory (and in the cache) ModelFinder joins merged partitions with '+'; best_scheme.nex uses '_'.
        modelos = cache[nome.replace("_", "+")]
        melhor = min(modelos, key=lambda m: -2 * modelos[m][0] + modelos[m][1] * math.log(n_sitios))
        lnl, df = modelos[melhor]
        soma_lnl += lnl
        soma_df += df
        n_total += n_sitios
        linhas_charset.append(f"  charset {nome} = {faixas};")
        linhas_modelo.append(f"{melhor}: {nome}")
        print(f"{len(modelos):4d} modelos  {n_sitios:6d} sitios  {melhor:22s} {nome[:70]}")

    df_total = soma_df - 1 + 111
    print(f"\nsoma lnL {soma_lnl:.3f} | df {df_total} | sitios {n_total}")
    if round(soma_lnl, 3) != LNL_ORIGINAL or df_total != DF_ORIGINAL:
        sys.exit(f"nao reproduz o original (LnL {LNL_ORIGINAL}, df {DF_ORIGINAL}); nada gravado")

    saida.write_text(
        "#nexus\nbegin sets;\n"
        + "\n".join(linhas_charset)
        + "\n  charpartition mymodels =\n    "
        + ",\n    ".join(linhas_modelo)
        + ";\nend;\n"
    )
    print(f"reproduz o original; gravado {saida}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
