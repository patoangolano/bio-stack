"""
f03c — CDS por loco dos taxons derivados de miniprot, com auditoria.

Le ~/rapinantes/miniprot/<taxon>.gff (miniprot --gff) e <taxon>.cds.fna (gffread -x)
e, para cada loco (Target da consulta), fica com o alinhamento Rank=1. O registro
sai do loco se:
  * Frameshift > 0 ou StopCodon > 0 (lesao no quadro de leitura: erro de montagem,
    pseudogene ou alinhamento ruim — nos tres casos o codon nao e confiavel);
  * cobertura da consulta < 50 %;
  * identidade com a consulta < 70 % (suspeita de paralogo: as consultas sao de
    parentes da mesma ordem, a 20-60 Ma);
  * o CDS extraido nao traduz limpo (fora de quadro ou stop interno).
Um segundo alinhamento com escore >= 0,99 do primeiro (o --outs padrao) fica
registrado como "ambiguo" e o loco fica com o Rank=1.

Saidas: ~/rapinantes/miniprot/cds_por_taxon/<taxon>.fna (cabecalho = loco) e
data/proc/rapinantes/filogenomica/auditoria_miniprot.tsv (+ resumo por taxon).

    uv run python scripts/rapinantes/filogenomica/f03c_miniprot_cds.py
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq

RAIZ = Path(__file__).resolve().parents[3]
PROC = RAIZ / "data" / "proc" / "rapinantes" / "filogenomica"
WSL = Path(r"\\wsl.localhost\Ubuntu-24.04\home\manec\rapinantes")
MP = WSL / "miniprot"
MIN_COB, MIN_IDENT = 0.50, 0.70


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


def traduz(cds: str) -> tuple[str | None, str]:
    cds = re.sub(r"[^ACGT]", "N", cds.upper())
    if len(cds) % 3:
        return None, "fora de quadro"
    prot = str(Seq(cds).translate(table="Standard"))
    if prot.endswith("*"):
        prot, cds = prot[:-1], cds[:-3]
    if "*" in prot:
        return None, "stop interno no CDS extraido"
    return cds, ""


def main() -> None:
    taxons = [t for t in csv.DictReader((PROC / "taxons.tsv").open(encoding="utf-8"), delimiter="\t") if t["fonte"] == "miniprot"]
    (MP / "cds_por_taxon").mkdir(exist_ok=True)
    auditoria, resumo = [], []
    for t in taxons:
        taxon, consulta = t["taxon"], t["consulta_miniprot"]
        gff, cds_f = MP / f"{taxon}.gff", MP / f"{taxon}.cds.fna"
        if not gff.is_file() or not cds_f.is_file():
            print(f"{taxon}: sem saida do miniprot")
            continue
        qlen = {k: len(v) for k, v in ler_fasta(WSL / "locos" / "consultas" / f"{consulta}.faa").items()}
        cds = ler_fasta(cds_f)
        alns: dict[str, list[dict[str, str]]] = defaultdict(list)
        for linha in gff.open(encoding="utf-8"):
            if linha.startswith("#"):
                continue
            c = linha.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "mRNA":
                continue
            at = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            alvo = at["Target"].split()
            at.update(loco=alvo[0], q_ini=alvo[1], q_fim=alvo[2], contig=c[0], ini=c[3], fim=c[4], fita=c[6])
            alns[at["loco"]].append(at)
        cont: dict[str, int] = defaultdict(int)
        with (MP / "cds_por_taxon" / f"{taxon}.fna").open("w", encoding="utf-8", newline="\n") as fh:
            for loco in sorted(qlen):
                lista = sorted(alns.get(loco, []), key=lambda a: int(a["Rank"]))
                if not lista:
                    cont["sem alinhamento"] += 1
                    auditoria.append({"taxon": taxon, "loco": loco, "status": "sem alinhamento"})
                    continue
                a = lista[0]
                cob = (int(a["q_fim"]) - int(a["q_ini"]) + 1) / qlen[loco]
                ident = float(a["Identity"])
                fs, st = int(a.get("Frameshift", 0)), int(a.get("StopCodon", 0))
                motivo = ""
                if fs:
                    motivo = "frameshift"
                elif st:
                    motivo = "stop interno"
                elif cob < MIN_COB:
                    motivo = "cobertura < 50 %"
                elif ident < MIN_IDENT:
                    motivo = "identidade < 70 %"
                seq = None
                if not motivo:
                    seq, motivo = traduz(cds.get(a["ID"], ""))
                status = motivo or "ok"
                cont[status] += 1
                if len(lista) > 1:
                    cont["ambiguo (2o alinhamento >= 0,99)"] += 1
                auditoria.append({
                    "taxon": taxon, "loco": loco, "status": status, "identidade": f"{ident:.3f}",
                    "cobertura": f"{cob:.3f}", "frameshift": fs, "stop": st, "n_alinhamentos": len(lista),
                    "contig": a["contig"], "inicio": a["ini"], "fim": a["fim"], "fita": a["fita"],
                })
                if seq:
                    fh.write(f">{loco}\n{seq}\n")
        ok = cont.get("ok", 0)
        resumo.append({"taxon": taxon, "acessao": t["acessao"], "consulta": consulta, "locos": len(qlen), "ok": ok,
                       **{k: v for k, v in sorted(cont.items()) if k != "ok"}})
        print(f"{taxon}: {ok}/{len(qlen)} ok; {dict(cont)}")
    campos = ["taxon", "loco", "status", "identidade", "cobertura", "frameshift", "stop", "n_alinhamentos", "contig", "inicio", "fim", "fita"]
    with (PROC / "auditoria_miniprot.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(auditoria)
    chaves = sorted({k for r in resumo for k in r}, key=lambda k: (k not in ("taxon", "acessao", "consulta", "locos", "ok"), k))
    with (PROC / "auditoria_miniprot_resumo.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=chaves, delimiter="\t", lineterminator="\n", restval=0)
        w.writeheader()
        w.writerows(resumo)


if __name__ == "__main__":
    main()
