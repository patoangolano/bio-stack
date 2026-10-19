"""
f02f — locos a partir dos ortologos NCBI (RefSeq) baixados pelo f02e.

Para cada gene candidato do galo, na ordem da permutacao do f02d:
  * junta os genes do mesmo grupo "NCBI Ortholog" nos 46 taxons RefSeq;
  * um registro por taxon: o CDS mais longo do gene (os relatorios de produto
    das aves nao trazem RefSeq Select); taxon com mais de um gene no grupo sai
    do loco (ambiguidade de ortologia);
  * CDS valido: multiplo de 3, sem stop interno (stop terminal removido) e
    com comprimento entre 0,5x e 1,5x a mediana do loco (modelo fundido ou
    fragmento sai do loco);
  * loco aprovado: >= 40 dos 46 taxons RefSeq (85 %), CDS do galo >= 600 pb e
    proteina de consulta disponivel para o miniprot.

Saidas (no WSL, ~/rapinantes/locos/): refseq_cds/<loco>.fna (cabecalho = taxon),
consultas/<taxon_consulta>.faa (proteina do taxon de consulta, cabecalho = loco).
No projeto: data/proc/rapinantes/filogenomica/genes_pool.tsv e auditoria_refseq.tsv.

    uv run python scripts/rapinantes/filogenomica/f02f_locos_refseq.py
"""

from __future__ import annotations

import csv
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq

RAIZ = Path(__file__).resolve().parents[3]
PROC = RAIZ / "data" / "proc" / "rapinantes" / "filogenomica"
WSL = Path(r"\\wsl.localhost\Ubuntu-24.04\home\manec\rapinantes\locos")
MIN_REFSEQ = 40  # 85 % de 46
MIN_CDS_GALO = 600

# taxon de consulta do miniprot -> alternativas, na ordem, se faltar no loco
CONSULTAS = {
    "Aquila_chrysaetos": ["Aquila_chrysaetos", "Haliaeetus_albicilla", "Astur_gentilis", "Buteo_buteo", "Harpia_harpyja"],
    "Falco_peregrinus": ["Falco_peregrinus", "Falco_naumanni"],
    "Gymnogyps_californianus": ["Gymnogyps_californianus", "Aquila_chrysaetos", "Haliaeetus_albicilla"],
    "Strix_aluco": ["Strix_aluco", "Athene_noctua", "Tyto_alba"],
}


def ler_fasta(caminho: Path) -> dict[str, tuple[str, str]]:
    """acessao -> (cabecalho, sequencia)."""
    seqs: dict[str, tuple[str, str]] = {}
    cab, partes = None, []
    with caminho.open(encoding="utf-8") as fh:
        for linha in fh:
            if linha.startswith(">"):
                if cab is not None:
                    seqs[cab.split()[0]] = (cab, "".join(partes))
                cab, partes = linha[1:].strip(), []
            else:
                partes.append(linha.strip())
    if cab is not None:
        seqs[cab.split()[0]] = (cab, "".join(partes))
    return seqs


def traduz(cds: str) -> str | None:
    """Proteina sem o stop terminal; None se o CDS nao for valido."""
    cds = cds.upper()
    if len(cds) % 3 or re.search(r"[^ACGTN]", cds):
        return None
    prot = str(Seq(cds).translate(table="Standard"))
    if prot.endswith("*"):
        prot = prot[:-1]
    if "*" in prot:
        return None
    return prot


def main() -> None:
    taxons = list(csv.DictReader((PROC / "taxons.tsv").open(encoding="utf-8"), delimiter="\t"))
    tax2nome = {t["tax_id"]: t["taxon"] for t in taxons if t["fonte"] == "RefSeq"}
    assert len(tax2nome) == 46, len(tax2nome)

    grupo_galo: dict[str, str] = {}
    for linha in (WSL / "galo_genes.jsonl").open(encoding="utf-8"):
        g = json.loads(linha)
        for gg in g.get("gene_groups", []):
            if gg.get("method") == "NCBI Ortholog":
                grupo_galo[g["gene_id"]] = gg["id"]

    # grupo -> taxon -> lista de (gene_id, cds) ; cds por GeneID (mais longo depois)
    por_grupo: dict[str, dict[str, dict[str, list[str]]]] = defaultdict(lambda: defaultdict(dict))
    for lote in sorted((WSL / "ortologos").glob("lote_*.d")):
        d = lote / "ncbi_dataset" / "data"
        gene_info: dict[str, tuple[str, str]] = {}
        for linha in (d / "data_report.jsonl").open(encoding="utf-8"):
            g = json.loads(linha)
            grp = next((x["id"] for x in g.get("geneGroups", []) if x.get("method") == "NCBI Ortholog"), None)
            if grp:
                gene_info[g["geneId"]] = (grp, g["taxId"])
        for cab, seq in ler_fasta(d / "cds.fna").values():
            m = re.search(r"\[GeneID=(\d+)\]", cab)
            if not m or m[1] not in gene_info:
                continue
            grp, tax = gene_info[m[1]]
            if tax not in tax2nome:
                continue
            por_grupo[grp][tax2nome[tax]].setdefault(m[1], []).append(seq)

    cands = list(csv.DictReader((PROC / "candidatos_permutados.tsv").open(encoding="utf-8"), delimiter="\t"))
    (WSL / "refseq_cds").mkdir(exist_ok=True)
    (WSL / "consultas").mkdir(exist_ok=True)
    consultas: dict[str, list[tuple[str, str]]] = {k: [] for k in CONSULTAS}
    pool, auditoria = [], []
    vistos_simbolo: set[str] = set()
    n_baixados = 0
    for c in cands:
        grp = grupo_galo.get(c["gene_id"])
        if grp is None or grp not in por_grupo:
            continue
        n_baixados += 1
        loco = re.sub(r"[^A-Za-z0-9_]", "_", c["simbolo"])
        registros: dict[str, tuple[str, str]] = {}
        motivo_taxon: dict[str, str] = {}
        for taxon, genes in por_grupo[grp].items():
            if len(genes) > 1:
                motivo_taxon[taxon] = f"{len(genes)} genes no grupo"
                continue
            (cdss,) = genes.values()
            validos = [(s, p) for s in cdss if (p := traduz(s))]
            if not validos:
                motivo_taxon[taxon] = "CDS invalido (stop interno ou fora de quadro)"
                continue
            s, p = max(validos, key=lambda sp: len(sp[0]))
            registros[taxon] = (s[: len(p) * 3], p)
        if registros:
            med = statistics.median(len(p) for _, p in registros.values())
            for taxon in list(registros):
                n = len(registros[taxon][1])
                if not 0.5 * med <= n <= 1.5 * med:
                    motivo_taxon[taxon] = f"comprimento {n} aa vs mediana {med:.0f}"
                    del registros[taxon]
        galo = registros.get("Gallus_gallus")
        consulta_ok = {k: next((t for t in alts if t in registros), None) for k, alts in CONSULTAS.items()}
        motivo = ""
        if len(registros) < MIN_REFSEQ:
            motivo = f"so {len(registros)} taxons RefSeq"
        elif galo is None or len(galo[0]) < MIN_CDS_GALO:
            motivo = "CDS do galo ausente ou < 600 pb"
        elif not all(consulta_ok.values()):
            motivo = "sem proteina de consulta"
        elif loco in vistos_simbolo:
            motivo = "simbolo repetido"
        auditoria.append({
            "ordem": c["ordem"], "loco": loco, "gene_id_galo": c["gene_id"], "grupo_ortologo": grp,
            "n_refseq": len(registros), "aprovado": "nao" if motivo else "sim", "motivo": motivo,
            "taxons_excluidos": "; ".join(f"{t}: {m}" for t, m in sorted(motivo_taxon.items())),
        })
        if motivo:
            continue
        assert galo is not None
        vistos_simbolo.add(loco)
        with (WSL / "refseq_cds" / f"{loco}.fna").open("w", encoding="utf-8", newline="\n") as fh:
            for taxon in sorted(registros):
                fh.write(f">{taxon}\n{registros[taxon][0]}\n")
        for k, t in consulta_ok.items():
            if t is not None:
                consultas[k].append((loco, registros[t][1]))
        pool.append({
            "ordem": c["ordem"], "loco": loco, "gene_id_galo": c["gene_id"], "grupo_ortologo": grp,
            "cromossomo_galo": c["cromossomo_galo"], "n_refseq": len(registros),
            "cds_galo_pb": len(galo[0]), "mediana_aa": int(statistics.median(len(p) for _, p in registros.values())),
            "consultas": ",".join(f"{k}={t}" for k, t in consulta_ok.items() if t != k),
        })
    for k, lista in consultas.items():
        with (WSL / "consultas" / f"{k}.faa").open("w", encoding="utf-8", newline="\n") as fh:
            for loco, p in lista:
                fh.write(f">{loco}\n{p}\n")
    for nome, linhas in (("genes_pool.tsv", pool), ("auditoria_refseq.tsv", auditoria)):
        if not linhas:
            continue
        with (PROC / nome).open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(linhas)
    motivos: dict[str, int] = defaultdict(int)
    for a in auditoria:
        motivos[re.sub(r"\d+", "N", a["motivo"]) or "aprovado"] += 1
    print(f"{n_baixados} candidatos com ortologos baixados; {len(pool)} aprovados; {dict(motivos)}")


if __name__ == "__main__":
    main()
