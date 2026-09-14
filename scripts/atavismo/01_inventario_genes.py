"""Inventario de presenca/ausencia dos modulos genicos do atavismo em Archosauria.

Pergunta: as aves perderam o *hardware* (matriz de esmalte/dentina) mas
mantiveram o *software* (rede regulatoria odontogenica e de alongamento axial)?

Consulta a API NCBI Datasets v2alpha (sem credencial; usa NCBI_API_KEY se houver)
e grava uma matriz presenca/ausencia por gene x taxon.

ATENCAO METODOLOGICA: ausencia na anotacao RefSeq NAO e prova de ausencia no
genoma. E evidencia de primeira ordem que exige confirmacao por busca em
sequencia (tBLASTn contra o genoma) antes de virar afirmacao no texto.
"""

from __future__ import annotations

import collections
import json
import os
import sys
import time
from pathlib import Path
from urllib.parse import quote

import requests

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
RAIZ = Path(__file__).resolve().parents[2]
SAIDA = RAIZ / "runs" / "atavismo"

# --- Taxons: o eixo filogenetico do problema ---------------------------------
# Archosauria = Aves + Crocodylia. Testudines e grupo-irmao de Archosauria.
# Lepidosauria polariza. Mamiferos sao o controle de toolkit completo.
TAXONS: list[tuple[str, str, str]] = [
    ("Gallus gallus", "Aves/Galliformes", "ave, dente ausente"),
    ("Melopsittacus undulatus", "Aves/Psittaciformes", "ave, dente ausente"),
    ("Taeniopygia guttata", "Aves/Passeriformes", "ave, dente ausente"),
    ("Alligator mississippiensis", "Crocodylia", "dente presente, polifiodonte"),
    ("Crocodylus porosus", "Crocodylia", "dente presente, polifiodonte"),
    ("Gavialis gangeticus", "Crocodylia", "dente presente, polifiodonte"),
    ("Chrysemys picta bellii", "Testudines", "CONTROLE: perda independente de dente"),
    ("Anolis carolinensis", "Lepidosauria", "dente presente (grupo externo)"),
    ("Homo sapiens", "Mammalia", "CONTROLE: toolkit completo"),
    ("Mus musculus", "Mammalia", "CONTROLE: toolkit completo"),
]

ORDEM_CLADOS = [
    "Aves/Galliformes",
    "Aves/Psittaciformes",
    "Aves/Passeriformes",
    "Crocodylia",
    "Testudines",
    "Lepidosauria",
    "Mammalia",
]

# --- Modulos genicos ---------------------------------------------------------
MODULOS: dict[str, dict] = {
    "A_esmalte_dentina_HARDWARE": {
        "hipotese": "perdido nas aves (pseudogenizado); intacto em Crocodylia",
        "genes": [
            "ENAM", "AMBN", "AMELX", "DSPP", "MMP20", "AMTN", "ODAM",
            "KLK4", "ACP4", "SLC24A4", "WDR72", "FAM83H",
        ],
    },
    "B_odontogenese_REGULATORIA": {
        "hipotese": "retido e funcional nas aves (pleiotropia: pena, bico, escama)",
        "genes": [
            "BMP4", "BMP2", "SHH", "PITX1", "PITX2", "MSX1", "MSX2",
            "BARX1", "RUNX2", "FGF8", "FGF4", "LEF1", "WNT10A", "WNT10B",
            "EDA", "EDAR", "EDARADD", "PAX9", "DLX1", "DLX2", "SP6",
            "SOSTDC1", "CTNNB1", "NOTCH1", "JAG1", "TFAP2B",
        ],
    },
    "C_eixo_caudal": {
        "hipotese": "presente nas aves; a diferenca e regulatoria/temporal, nao de conteudo",
        "genes": [
            "TBXT", "GDF11", "LIN28A", "LIN28B", "CDX1", "CDX2", "CDX4",
            "HOXA13", "HOXB13", "HOXC13", "HOXD13", "WNT3A", "WNT5A",
            "ALDH1A2", "CYP26A1", "MSGN1", "TBX6", "RSPO3", "SALL4",
            "HES7", "DLL1", "LFNG", "MESP2", "RIPPLY2",
        ],
    },
    "D_ciliogenese_talpid": {
        "hipotese": "intacto; a via talpid2/C2CD3 revelou germe dentario em frango",
        "genes": [
            "C2CD3", "KIAA0586", "IFT88", "OFD1", "GLI3", "SMO", "KIF7",
            "PTCH1", "EVC", "EVC2",
        ],
    },
}


def consulta(simbolos: list[str], taxon: str) -> dict[str, dict]:
    """Consulta simbolos em lote para um taxon. Devolve {SIMBOLO: registro}."""
    url = f"{API}/gene/symbol/{','.join(simbolos)}/taxon/{quote(taxon)}"
    cab = {"Accept": "application/json"}
    chave = os.environ.get("NCBI_API_KEY")
    if chave:
        cab["api-key"] = chave

    for tentativa in range(4):
        try:
            r = requests.get(url, headers=cab, timeout=60)
            if r.status_code == 200:
                achados: dict[str, dict] = {}
                for rep in r.json().get("reports", []):
                    g = rep.get("gene", {})
                    sim = g.get("symbol")
                    if sim:
                        achados.setdefault(sim.upper(), g)
                return achados
            if r.status_code in (429, 500, 502, 503):
                time.sleep(2**tentativa)
                continue
            print(f"  ! HTTP {r.status_code} em {taxon}", file=sys.stderr)
            return {}
        except requests.RequestException as e:
            print(f"  ! {type(e).__name__} em {taxon} ({tentativa + 1}/4)", file=sys.stderr)
            time.sleep(2**tentativa)
    return {}


def classifica(g: dict | None) -> str:
    """Estado compacto de um gene num taxon."""
    if not g:
        return "-"
    tipo = (g.get("type") or "").upper()
    if "PSEUDO" in tipo:
        return "PSI"
    if tipo == "PROTEIN_CODING":
        n = g.get("protein_count") or 0
        return f"P{n}" if n else "P?"
    return tipo[:4] or "?"


def main() -> int:
    SAIDA.mkdir(parents=True, exist_ok=True)
    linhas: list[dict] = []
    bruto: dict[str, dict] = {}

    for modulo, info in MODULOS.items():
        genes: list[str] = info["genes"]
        print(f"\n### {modulo}  ({len(genes)} genes)")
        for taxon, clado, nota in TAXONS:
            achados: dict[str, dict] = {}
            for i in range(0, len(genes), 8):
                achados.update(consulta(genes[i : i + 8], taxon))
                time.sleep(0.15)
            presentes = sum(1 for s in genes if achados.get(s.upper()))
            print(f"  {taxon:<30s} {presentes:2d}/{len(genes)}  ({clado})")
            for s in genes:
                g = achados.get(s.upper())
                linhas.append(
                    {
                        "modulo": modulo,
                        "gene": s,
                        "taxon": taxon,
                        "clado": clado,
                        "nota_taxon": nota,
                        "estado": classifica(g),
                        "gene_id": (g or {}).get("gene_id", ""),
                        "cromossomo": ",".join((g or {}).get("chromosomes", []) or []),
                    }
                )
                if g:
                    bruto[f"{taxon}|{s}"] = g

    campos = ["modulo", "gene", "taxon", "clado", "nota_taxon", "estado",
              "gene_id", "cromossomo"]
    tsv = SAIDA / "inventario_genes.tsv"
    with tsv.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(campos) + "\n")
        for ln in linhas:
            fh.write("\t".join(str(ln[c]) for c in campos) + "\n")

    (SAIDA / "inventario_bruto.json").write_text(
        json.dumps(bruto, indent=1, ensure_ascii=False), encoding="utf-8"
    )

    md = SAIDA / "matriz_presenca.md"
    with md.open("w", encoding="utf-8") as fh:
        fh.write("# Matriz de presenca dos modulos do atavismo\n\n")
        fh.write(
            "Celula = taxons com o gene presente / taxons consultados no clado.\n\n"
            "Ausente da anotacao RefSeq **nao** e o mesmo que ausente do genoma: "
            "e evidencia de primeira ordem, a confirmar por tBLASTn.\n"
        )
        for modulo, info in MODULOS.items():
            fh.write(f"\n## {modulo}\n\n_Hipotese: {info['hipotese']}_\n\n")
            fh.write("| gene | " + " | ".join(ORDEM_CLADOS) + " |\n")
            fh.write("|" + "---|" * (len(ORDEM_CLADOS) + 1) + "\n")
            for gene in info["genes"]:
                cel = []
                for clado in ORDEM_CLADOS:
                    est = [
                        ln["estado"] for ln in linhas
                        if ln["gene"] == gene and ln["clado"] == clado
                        and ln["modulo"] == modulo
                    ]
                    pres = [e for e in est if e != "-"]
                    cel.append(f"{len(pres)}/{len(est)}" if est else "?")
                fh.write(f"| {gene} | " + " | ".join(cel) + " |\n")

            agg: dict[str, list[int]] = collections.defaultdict(lambda: [0, 0])
            for ln in linhas:
                if ln["modulo"] != modulo:
                    continue
                agg[ln["clado"]][1] += 1
                if ln["estado"] != "-":
                    agg[ln["clado"]][0] += 1
            fh.write("\n**Fracao do modulo presente, por clado:**\n\n")
            for clado in ORDEM_CLADOS:
                if clado in agg:
                    p, t = agg[clado]
                    fh.write(f"- {clado}: **{p}/{t}** ({100 * p / t:.0f}%)\n")

    print(f"\nEscrito: {tsv}")
    print(f"Escrito: {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
