"""
f04 — locos finais: CDS RefSeq (f02f) + CDS miniprot (f03c), filtro de ocupacao e sorteio.

Ordem de entrada = permutacao fixa do f02d (rng = np.random.default_rng(0)); ficam
os primeiros N locos com >= 85 % dos 57 taxons (>= 49). Registro miniprot com
comprimento fora de 0,5x-1,5x a mediana do loco sai (fragmento ou fusao).

Saidas:
  data/proc/rapinantes/filogenomica/locos_aa/<loco>.faa   proteina (cabecalho = taxon)
  ~/rapinantes/locos/cds_final/<loco>.fna                 CDS sem stop, mesmo conjunto
  data/proc/rapinantes/filogenomica/genes.tsv             lista de genes com ocupacao
  data/proc/rapinantes/filogenomica/ocupacao_taxon.tsv    locos por taxon

    uv run python scripts/rapinantes/filogenomica/f04_locos_finais.py <N>
"""

from __future__ import annotations

import csv
import hashlib
import statistics
import sys
from pathlib import Path

from Bio.Seq import Seq

RAIZ = Path(__file__).resolve().parents[3]
PROC = RAIZ / "data" / "proc" / "rapinantes" / "filogenomica"
WSL = Path(r"\\wsl.localhost\Ubuntu-24.04\home\manec\rapinantes")


def ler_fasta(caminho: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    nome, partes = None, []
    with caminho.open(encoding="utf-8") as fh:
        for linha in fh:
            if linha.startswith(">"):
                if nome is not None:
                    seqs[nome] = "".join(partes)
                nome, partes = linha[1:].split()[0], []
            else:
                partes.append(linha.strip())
    if nome is not None:
        seqs[nome] = "".join(partes)
    return seqs


def main() -> None:
    n_final = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    taxons = list(csv.DictReader((PROC / "taxons.tsv").open(encoding="utf-8"), delimiter="\t"))
    nomes = [t["taxon"] for t in taxons]
    mp = [t["taxon"] for t in taxons if t["fonte"] == "miniprot"]
    min_taxa = -(-85 * len(nomes) // 100)  # teto de 85 %
    mp_cds = {t: ler_fasta(f) for t in mp if (f := WSL / "miniprot" / "cds_por_taxon" / f"{t}.fna").is_file()}
    pool = list(csv.DictReader((PROC / "genes_pool.tsv").open(encoding="utf-8"), delimiter="\t"))

    saida_aa = PROC / "locos_aa"
    saida_cds = WSL / "locos" / "cds_final"
    for d in (saida_aa, saida_cds):
        d.mkdir(parents=True, exist_ok=True)
        for f in d.glob("*.f*a"):
            f.unlink()
    finais, ocup = [], dict.fromkeys(nomes, 0)
    fora_mp = 0
    for g in pool:
        if len(finais) >= n_final:
            break
        reg = ler_fasta(WSL / "locos" / "refseq_cds" / f"{g['loco']}.fna")
        med = statistics.median(len(s) for s in reg.values())
        n_mp = 0
        for t, d in mp_cds.items():
            s = d.get(g["loco"])
            if s is None:
                continue
            if not 0.5 * med <= len(s) <= 1.5 * med:
                fora_mp += 1
                continue
            reg[t] = s
            n_mp += 1
        if len(reg) < min_taxa:
            continue
        with (saida_cds / f"{g['loco']}.fna").open("w", encoding="utf-8", newline="\n") as fc, \
             (saida_aa / f"{g['loco']}.faa").open("w", encoding="utf-8", newline="\n") as fa:
            for t in nomes:
                if t not in reg:
                    continue
                prot = str(Seq(reg[t]).translate(table="Standard"))
                assert "*" not in prot, (g["loco"], t)
                fc.write(f">{t}\n{reg[t]}\n")
                fa.write(f">{t}\n{prot}\n")
                ocup[t] += 1
        finais.append({**{k: g[k] for k in ("ordem", "loco", "gene_id_galo", "grupo_ortologo", "cromossomo_galo", "cds_galo_pb", "mediana_aa")},
                       "n_taxons": len(reg), "n_refseq": int(g["n_refseq"]), "n_miniprot": n_mp})
    with (PROC / "genes.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(finais[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(finais)
    with (PROC / "ocupacao_taxon.tsv").open("w", encoding="utf-8", newline="") as fh:
        fh.write("taxon\tfonte\tlocos\tfracao\n")
        for t in taxons:
            fh.write(f"{t['taxon']}\t{t['fonte']}\t{ocup[t['taxon']]}\t{ocup[t['taxon']] / len(finais):.3f}\n")
    h = hashlib.sha256()
    for f in sorted(saida_aa.glob("*.faa")):
        h.update(f.read_bytes())
    total = sum(r["n_taxons"] for r in finais)
    print(f"{len(finais)} locos (pedido {n_final}; minimo {min_taxa}/{len(nomes)} taxons; pool {len(pool)}); "
          f"ocupacao {total / (len(finais) * len(nomes)):.3f}; registros miniprot fora por comprimento: {fora_mp}; "
          f"sha256(locos_aa) {h.hexdigest()[:16]}")


if __name__ == "__main__":
    main()
