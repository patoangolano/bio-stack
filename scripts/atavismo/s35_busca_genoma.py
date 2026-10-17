"""Le a busca de genoma inteiro do s35_busca_genoma.sh.

HSPs da mesma consulta na mesma sequencia a menos de JUNTA bp viram um LOCUS;
o escore do locus e a soma dos bits. Por consulta e genoma, os 3 melhores loci.

CONTROLE: o melhor locus de ODAPH no jacare tem de cair no LOC106740521
(NC_081825.1:164.672.404-164.676.187). Nas aves, cada locus com e <= 1e-3 traz
o melhor hit reciproco no proteoma do jacare (s35_reciproco.tsv): so conta como
ODAPH se o reciproco for a proteina do jacare XP_019349100.2 ("uncharacterized
protein LOC106740521").

    uv run python scripts/atavismo/s35_busca_genoma.py
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

SAIDA = Path(__file__).resolve().parents[2] / "runs" / "atavismo" / "genoma"
JUNTA = 50_000
ODAPH_JACARE = ("NC_081825.1", 164_672_404, 164_676_187)
REC_ODAPH = "LOC106740521"


def loci(tsv: Path) -> dict[str, list[dict]]:
    por_q: dict[str, list[list]] = defaultdict(list)
    with tsv.open(encoding="utf-8") as f:
        for c in csv.reader(f, delimiter="\t"):
            # com -parse_seqids o sseqid vem "ref|NC_081825.1|": o nome e o pedaco nao vazio
            q, s = c[0], next(p for p in reversed(c[1].split("|")) if p)
            lo, hi = sorted((int(c[6]), int(c[7])))
            por_q[q].append([s, lo, hi, float(c[8]), float(c[9])])
    saida: dict[str, list[dict]] = {}
    for q, hsps in por_q.items():
        hsps.sort()
        grupos: list[dict] = []
        for s, lo, hi, ev, bits in hsps:
            g = grupos[-1] if grupos else None
            if g and g["seq"] == s and lo - g["fim"] <= JUNTA:
                g["fim"] = max(g["fim"], hi)
                g["bits"] += bits
                g["ev"] = min(g["ev"], ev)
                g["n"] += 1
            else:
                grupos.append({"seq": s, "ini": lo, "fim": hi, "bits": bits, "ev": ev, "n": 1})
        saida[q] = sorted(grupos, key=lambda g: -g["bits"])
    return saida


def main() -> int:
    rec: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    if (SAIDA / "s35_reciproco.tsv").exists():
        with (SAIDA / "s35_reciproco.tsv").open(encoding="utf-8") as f:
            for c in csv.reader(f, delimiter="\t"):
                rec[(c[0], c[1], c[2])].append(f"{c[3]}-{c[4]}|{c[7]}")
    linhas = []
    for rot in ("Amis", "Ggal", "Shab"):
        tsv = SAIDA / f"s35_{rot}.tsv"
        if not tsv.exists():
            continue
        lc = loci(tsv)
        print(f"\n===== {rot} =====")
        for q in sorted(lc, key=lambda x: (x.split("_")[0], x)):
            for i, g in enumerate(lc[q][:3]):
                nota = ""
                if rot == "Amis" and q.startswith("ODAPH") and i == 0:
                    s, a, b = ODAPH_JACARE
                    ok = g["seq"] == s and g["ini"] <= b and g["fim"] >= a
                    nota = "CONTROLE OK (LOC106740521)" if ok else "CONTROLE FALHOU"
                elif g["bits"] >= 30 and (rot != "Amis" or q.startswith("ODAM")):
                    vistos = [r for (r_rot, r_q, r_s), lst in rec.items()
                              if r_rot == rot and r_q == q and r_s == g["seq"] for r in lst
                              if g["ini"] <= int(r.split("-")[0]) <= g["fim"]]
                    if vistos:
                        alvo = vistos[0].split("|", 1)[1]
                        nota = ("RECIPROCO = ODAPH do jacare" if REC_ODAPH in alvo else f"reciproco: {alvo[:70]}")
                print(f"  {q:11s} #{i + 1} {g['seq']:>16s}:{g['ini']:>11,}-{g['fim']:<11,} "
                      f"{g['n']:2d} HSP  {g['bits']:7.1f} bits  e={g['ev']:.1e}  {nota}")
                linhas.append({"genoma": rot, "consulta": q, "rank": i + 1, "seq": g["seq"],
                               "ini": g["ini"], "fim": g["fim"], "hsps": g["n"], "bits": round(g["bits"], 1),
                               "evalue": g["ev"], "nota": nota})
    if linhas:
        with (SAIDA / "s35_loci.tsv").open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(linhas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
