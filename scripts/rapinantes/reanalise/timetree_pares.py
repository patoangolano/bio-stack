"""
TimeTree 5 (API REST, timetree.temple.edu/api): idades de divergencia dos pares
de rapinantes pedidos.

    uv run python scripts/rapinantes/reanalise/timetree_pares.py

Para cada par: NCBI taxid por E-utilities (esearch, db=taxonomy, nome
cientifico exato), depois /pairwise/<a>/<b>/json. Grava a resposta bruta e uma
tabela com: idade "estimada" do TimeTree (precomputed/adjusted_age, com o IC
que ele publica), mediana do TimeTree (sum_median_time), e - recalculados aqui
a partir de studies.time_estimates - n de estimativas, mediana, quartis e
amplitude. Somente biblioteca padrao.
"""

from __future__ import annotations

import json
import statistics
import time
import urllib.parse
import urllib.request
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]
OUT = RAIZ / "runs" / "rapinantes" / "reanalise" / "timetree"
TT = "https://timetree.temple.edu/api/pairwise"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

PARES = [
    ("Accipitridae", "Cathartidae"),
    ("Accipitridae", "Sagittariidae"),
    ("Accipitridae", "Pandionidae"),
    ("Falconidae", "Psittaciformes"),
    ("Falconidae", "Cariamidae"),
    ("Tytonidae", "Strigidae"),
    ("Strigiformes", "Accipitriformes"),
    ("Herpetotheres", "Falco"),
    ("Caracara", "Falco"),
]


def get(url: str) -> str:
    for k in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "bio-stack"}),
                                        timeout=120) as r:
                return r.read().decode()
        except Exception as e:  # noqa: BLE001 - tenta de novo; a ultima falha sobe
            if k == 4:
                raise
            print(f"  {url}: {e}; nova tentativa")
            time.sleep(3 * (k + 1))
    raise RuntimeError


def taxid(nome: str) -> int:
    # "AND Aves[Subtree]": Sagittariidae tem homonimo fora das aves no NCBI (2 taxids sem o filtro)
    q = urllib.parse.urlencode({"db": "taxonomy", "term": f"{nome}[Scientific Name] AND Aves[Subtree]",
                                "retmode": "json"})
    d = json.loads(get(f"{EUTILS}?{q}"))
    ids = d["esearchresult"]["idlist"]
    if len(ids) != 1:
        raise SystemExit(f"{nome}: {len(ids)} taxids ({ids})")
    time.sleep(0.4)
    return int(ids[0])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ids = {n: taxid(n) for n in sorted({x for p in PARES for x in p})}
    (OUT / "ncbi_taxids.json").write_text(json.dumps(ids, indent=1))
    linhas = []
    for a, b in PARES:
        url = f"{TT}/{ids[a]}/{ids[b]}/json"
        try:
            txt = get(url)
            d = json.loads(txt)
        except Exception as e:  # noqa: BLE001 - documenta a falha na tabela
            linhas.append({"par": f"{a}-{b}", "erro": str(e), "url": url})
            continue
        (OUT / f"pairwise_{a}_{b}.json").write_text(txt)
        st = d.get("studies") or {}
        ests = [float(x) for x in (st.get("time_estimates") or "{}").strip("{}").split(",") if x.strip()]
        q = statistics.quantiles(ests, n=4) if len(ests) >= 2 else [float("nan")] * 3
        linhas.append({
            "par": f"{a}-{b}", "taxid_a": ids[a], "taxid_b": ids[b],
            "no_ncbi": (d.get("ncbi_crown") or {}).get("scientific_name", ""),
            "n_estudos_all_total": d.get("all_total"),
            "n_estimativas": len(ests),
            "idade_estimada_TT_Ma": round(float(d.get("sum_simple_mol_time") or st.get("precomputed_age") or "nan"), 2),
            "IC_TT_baixo": round(float(st.get("precomputed_ci_low") or "nan"), 2),
            "IC_TT_alto": round(float(st.get("precomputed_ci_high") or "nan"), 2),
            "mediana_TT_Ma": d.get("sum_median_time"),
            "mediana_recalc_Ma": round(statistics.median(ests), 2) if ests else "",
            "Q1_Ma": round(q[0], 2), "Q3_Ma": round(q[2], 2),
            "min_Ma": round(min(ests), 2) if ests else "", "max_Ma": round(max(ests), 2) if ests else "",
            "conflitos": len(d.get("conflict_records") or []),
            "versao_TT": d.get("version"), "url": url,
        })
        print(linhas[-1])
        time.sleep(1)
    cab = []
    for li in linhas:
        for k in li:
            if k not in cab:
                cab.append(k)
    with open(OUT / "timetree_pares.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cab) + "\n")
        for li in linhas:
            f.write("\t".join(str(li.get(c, "")) for c in cab) + "\n")


if __name__ == "__main__":
    main()
