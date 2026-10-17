"""Validacao independente por alinhamento de genoma inteiro: multiz100way do hg38.

Tudo o que o s26-s31 mediu vem de busca dirigida (tBLASTn, miniprot, blastn)
dentro de intervalos sintenicos. O multiz100way do UCSC e outro metodo inteiro:
cadeias e redes lastz genoma x genoma, costuradas pelo multiz. Ele inclui, contra
o humano: Alligator (allMis1), tres tartarugas (chrPic2, cheMyd1, pelSin1),
Anolis (anoCar2) e varias aves (galGal4, melGal1, anaPla1, taeGut2, melUnd1,
falChe1, falPer1, colLiv1, ...).

O Table Browser bloqueia acesso automatizado (Cloudflare Turnstile) e o PECAN de
amniotas do Ensembl so traz primatas nesta regiao; resta o MAF do chr4 inteiro
(4,6 GB, s36_multiz.sh), recortado aqui nos exons codificantes.

    uv run python scripts/atavismo/s36_multiz.py bed        (exons de ENAM/AMBN/AMTN)
    wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s36_multiz.sh
    uv run python scripts/atavismo/s36_multiz.py resumo

RESUMO por gene x montagem: fracao do CDS humano com base alinhada, e
frameshifts — corridas de gap (onde o humano tem base) ou de insercao, de
tamanho nao multiplo de 3, DENTRO dos blocos. Controles: allMis1 e anoCar2 tem
o gene funcional e devem ter cobertura alta e frameshift ~0.

MEDIDO em 2026-10-17: a contagem de frameshift REPROVA no controle — o ENAM
funcional do jacare sai com 54 contra o humano (~320 Ma: alinhamento de
nucleotideo mamifero x reptil cheio de gap que nao e lesao). So a COBERTURA e
usada; lesao continua vindo do s29/s31, onde o controle zerou.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import requests

OUT = Path(__file__).resolve().parents[2] / "runs" / "atavismo" / "multiz"
GENES = {"ENAM": "NM_031889.3", "AMBN": "NM_016519.6", "AMTN": "NM_212557.4"}
DENTE = {"allMis1": "COM", "anoCar2": "COM", "chrPic2": "sem", "cheMyd1": "sem", "pelSin1": "sem"}
AVES = {"galGal4", "melGal1", "anaPla1", "taeGut2", "melUnd1", "falChe1", "falPer1", "colLiv1",
        "ficAlb2", "zonAlb1", "geoFor1", "pseHum1", "amaVit1", "araMac1"}


def bed() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    r = requests.get("https://api.genome.ucsc.edu/getData/track",
                     params={"genome": "hg38", "track": "ncbiRefSeqCurated", "chrom": "chr4",
                             "start": 70_500_000, "end": 70_700_000}, timeout=120)
    r.raise_for_status()
    linhas = []
    for g in r.json()["ncbiRefSeqCurated"]:
        gene = next((k for k, v in GENES.items() if g["name"] == v), None)
        if not gene:
            continue
        ini = [int(x) for x in g["exonStarts"].strip(",").split(",")]
        fim = [int(x) for x in g["exonEnds"].strip(",").split(",")]
        for i, (a, b) in enumerate(zip(ini, fim, strict=True)):
            a, b = max(a, int(g["cdsStart"])), min(b, int(g["cdsEnd"]))
            if a < b:
                linhas.append(f"chr4\t{a}\t{b}\t{gene}_cds{i + 1}\n")
    (OUT / "cds_hg38.bed").write_text("".join(linhas), encoding="utf-8", newline="\n")
    print(f"{len(linhas)} exons codificantes em {OUT / 'cds_hg38.bed'}")
    return 0


def blocos(maf: Path):
    """Gera blocos do MAF como dict montagem -> (inicio no humano, texto alinhado)."""
    bloco: dict[str, str] = {}
    ini_hs = 0
    for linha in maf.read_text(encoding="utf-8").splitlines():
        if linha.startswith("a"):
            if bloco:
                yield ini_hs, bloco
            bloco = {}
        elif linha.startswith("s "):
            c = linha.split()
            mont = c[1].split(".")[0]
            if mont == "hg38":
                ini_hs = int(c[2])
            bloco[mont] = c[6].upper()
    if bloco:
        yield ini_hs, bloco


def resumo() -> int:
    exons = []
    for c in csv.reader((OUT / "cds_hg38.bed").open(encoding="utf-8"), delimiter="\t"):
        exons.append((c[3].split("_")[0], int(c[1]), int(c[2])))
    tam_cds = defaultdict(int)
    for g, a, b in exons:
        tam_cds[g] += b - a

    def gene_de(pos: int) -> str | None:
        return next((g for g, a, b in exons if a <= pos < b), None)

    cob: dict[tuple[str, str], int] = defaultdict(int)
    fs: dict[tuple[str, str], int] = defaultdict(int)
    for ini_hs, bl in blocos(OUT / "cds_hg38.maf"):
        hs = bl.get("hg38")
        if not hs:
            continue
        pos = []  # posicao no humano de cada coluna (None em coluna de gap humano)
        p = ini_hs
        for ch in hs:
            pos.append(p if ch != "-" else None)
            if ch != "-":
                p += 1
        for mont, seq in bl.items():
            if mont == "hg38":
                continue
            corrida, tipo, gene_corr = 0, "", None
            for k, ch in enumerate(seq):
                g = gene_de(pos[k]) if pos[k] is not None else gene_corr
                if pos[k] is not None and ch != "-" and g:
                    cob[(g, mont)] += 1
                t = "del" if (pos[k] is not None and ch == "-") else "ins" if (pos[k] is None and ch != "-") else ""
                if t and t == tipo:
                    corrida += 1
                else:
                    if corrida % 3 and gene_corr:
                        fs[(gene_corr, mont)] += 1
                    corrida, tipo = (1, t) if t else (0, "")
                gene_corr = g or gene_corr
            if corrida % 3 and gene_corr:
                fs[(gene_corr, mont)] += 1

    montagens = sorted({m for _, m in cob}, key=lambda m: (m not in DENTE, m not in AVES, m))
    linhas = []
    for g in GENES:
        print(f"\n=== {g} (CDS humano {tam_cds[g]} nt) ===")
        print(f"{'montagem':10s} {'grupo':10s} {'cobertura':>9s} {'frameshifts':>11s}")
        for m in montagens:
            if m not in DENTE and m not in AVES:
                continue
            grupo = ("ave" if m in AVES else "réptil " + DENTE[m])
            c = cob.get((g, m), 0) / tam_cds[g]
            print(f"{m:10s} {grupo:10s} {c:9.0%} {fs.get((g, m), 0):11d}")
            linhas.append({"gene": g, "montagem": m, "grupo": grupo, "cobertura": round(c, 3),
                           "frameshifts": fs.get((g, m), 0)})
    with (OUT / "s36_multiz_resumo.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    return 0


if __name__ == "__main__":
    modo = sys.argv[1] if len(sys.argv) > 1 else ""
    if modo not in ("bed", "resumo"):
        raise SystemExit("uso: s36_multiz.py bed|resumo")
    raise SystemExit(bed() if modo == "bed" else resumo())
