# /// script
# requires-python = ">=3.11"
# dependencies = ["requests>=2.31", "pandas>=2.0", "xlrd>=2.0.1"]
# ///
"""Taxa de recombinacao do chrZ da galinha na janela da reliquia de esmalte.

A reliquia (ENAM/AMBN, s25/s26) esta em NC_052572.1:55,15-55,17 Mb (GRCg7b),
no intervalo NRG1..LPL. Se a janela tiver recombinacao baixa, a persistencia
da reliquia se explica por ligacao e deriva, sem selecao — o confundidor do E1.

Fonte: Elferink et al. 2010, BMC Genetics 11:11 (PMID 20141624, PMC2834597),
Additional file 5 (1471-2156-11-11-S5.XLS): bins de ~500 kb, cM e cM/Mb por
sexo e por populacao, coordenadas WASHUC2 (galGal3). O Z aparece como "GGZ" e
so tem taxa de MACHO — e o que se espera: o Z so recombina em machos ZZ.
Baixado pela API de suplementos do Europe PMC (o site da editora deu 403).

A janela em galGal3 NAO e assumida: o script consulta NRG1 e LPL no refGene
do galGal3 pela API do UCSC e usa as coordenadas devolvidas.

    uv run --script scripts/atavismo/s28_recombinacao_z.py
"""

from __future__ import annotations

import io
import sys
import zipfile
from pathlib import Path

import pandas as pd
import requests

PMCID = "PMC2834597"
ARQ = "1471-2156-11-11-S5.XLS"
OUT = Path(__file__).resolve().parents[2] / "runs" / "atavismo" / "recombinacao"
UCSC = "https://api.genome.ucsc.edu/getData/track"

# colunas do S5.XLS (linha 0 = grupo, linha 1 = cM / cM/Mb)
COL = {"crom": 0, "bin": 1, "de": 4, "ate": 5, "tam_mb": 6,
       "macho_cm": 14, "macho_cmmb": 15, "pop1_macho_cmmb": 24, "pop2_macho_cmmb": 33}


def baixa() -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    destino = OUT / ARQ
    if not destino.exists():
        r = requests.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles",
                         timeout=300)
        r.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
            destino.write_bytes(zf.read(ARQ))
    return destino


def gene_galgal3(simbolo: str) -> tuple[int, int] | None:
    """Menor inicio e maior fim do simbolo no refGene do galGal3, chrZ 50-56 Mb."""
    r = requests.get(UCSC, params={"genome": "galGal3", "track": "refGene", "chrom": "chrZ",
                                   "start": 50_000_000, "end": 56_000_000}, timeout=120)
    r.raise_for_status()
    linhas = [g for g in r.json().get("refGene", []) if g.get("name2", "").upper() == simbolo]
    if not linhas:
        return None
    return min(int(g["txStart"]) for g in linhas), max(int(g["txEnd"]) for g in linhas)


def tabela() -> pd.DataFrame:
    bruto = pd.read_excel(baixa(), header=None, engine="xlrd").iloc[2:]
    df = pd.DataFrame({k: bruto[c] for k, c in COL.items()})
    for k in df.columns:
        if k != "crom":
            df[k] = pd.to_numeric(df[k], errors="coerce")
    df["crom"] = df["crom"].astype(str).str.strip()
    return df.dropna(subset=["de", "ate", "macho_cmmb"])


def percentil(valor: float, serie: pd.Series) -> float:
    v = [float(x) for x in serie]
    return 100.0 * sum(x <= valor for x in v) / len(v)


def main() -> int:
    df = tabela()
    z = pd.DataFrame(df.loc[df["crom"] == "GGZ"])
    auto = pd.DataFrame(df.loc[df["crom"].str.match(r"^GGA\d+$")])

    nrg1, lpl = gene_galgal3("NRG1"), gene_galgal3("LPL")
    print(f"galGal3 (UCSC refGene): NRG1 {nrg1}  LPL {lpl}")
    if not (nrg1 and lpl):
        print("NRG1/LPL nao achados no refGene do galGal3 — sem janela confiavel")
        return 1
    ini, fim = min(nrg1[0], lpl[0]), max(nrg1[1], lpl[1])
    print(f"janela NRG1..LPL em galGal3: {ini:,}-{fim:,} ({(fim - ini) / 1e6:.2f} Mb)")

    na_janela = pd.DataFrame(z.loc[(z["ate"] >= ini) & (z["de"] <= fim)])
    print("\nbins de ~500 kb que tocam a janela (taxa de MACHO, cM/Mb):")
    print(na_janela[["bin", "de", "ate", "macho_cm", "macho_cmmb", "pop1_macho_cmmb",
                     "pop2_macho_cmmb"]].to_string(index=False))

    # taxa da janela ponderada pelo tamanho dos bins (cM somados / Mb somados)
    taxa = float(na_janela["macho_cm"].sum()) / float((na_janela["ate"] - na_janela["de"]).sum() / 1e6)
    z_media = float(z["macho_cm"].sum()) / float((z["ate"] - z["de"]).sum() / 1e6)
    a_media = float(auto["macho_cm"].sum()) / float((auto["ate"] - auto["de"]).sum() / 1e6)
    zeros = float((z["macho_cm"] == 0).mean()) * 100
    print(f"\ntaxa na janela (macho, ponderada): {taxa:.2f} cM/Mb")
    print(f"media do Z (macho): {z_media:.2f} cM/Mb  |  media dos autossomos (macho): {a_media:.2f} cM/Mb")
    print(f"percentil da janela entre os {len(z)} bins do Z: "
          f"{percentil(taxa, pd.Series(z['macho_cmmb'])):.0f}  |  entre os {len(auto)} bins autossomicos: "
          f"{percentil(taxa, pd.Series(auto['macho_cmmb'])):.0f}")
    print(f"bins do Z com 0 cM: {zeros:.0f}%  (resolucao: 500 kb e numero finito de meioses —"
          " um 0 isolado nao e deserto de recombinacao)")

    saida = OUT / "s28_janela_reliquia.tsv"
    na_janela.to_csv(saida, sep="\t", index=False)
    print(f"\ntabela: {saida}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
