"""A reliquia de esmalte e de Neornithes, ou e particular de Gallus?

Estende o s24/s25 para o conjunto da secao 2 de docs/PROJETO-ARCHOSAURIA.md.
O desenho deixa de ser um par e passa a ser comparativo, com DUAS origens
independentes de perda de dentição e TRES controles positivos:

    SEM DENTE — perda 1 (Neornithes, ~116 Ma, origem unica)
        Gallus gallus             GCF_016699485.2   ja medido no s24
        Melopsittacus undulatus   GCF_012275295.1
        Strigops habroptila       GCF_004027225.2
    SEM DENTE — perda 2 (Testudines, INDEPENDENTE, Meredith et al. 2014)
        Chrysemys picta           GCF_011386835.1
    COM DENTE — controle positivo
        Alligator mississippiensis  consulta; os genes estao intactos nele
        Gavialis gangeticus         GCA_030020295.1
        Anolis carolinensis         GCF_035594765.1

A tartaruga e a peca que vale mais: ela entrou no projeto como grupo externo
para polarizar, mas por Meredith 2014 ela e um evento de perda independente.
Se a reliquia aparecer nas aves E na tartaruga, sao duas perdas independentes
deixando rastro comparavel. Se aparecer so nas aves, o rastro e de Neornithes.

METODO. Para cada especie: localizar as ancoras NRG1 e LPL, extrair o intervalo
entre elas com folga, e rodar tblastn das proteinas do jacare contra ele. As
ancoras sao as mesmas verificadas em 2026-10-12:

    NRG1 -> [ENAM . AMBN . AMTN] -> LPL -> PSD3   conservado em Amis e Ggal

Tambem registra, por especie, se os nove genes do painel estao ANOTADOS. Isso
e informacao por si: anotado significa que ainda ha gene reconhecivel.

    uv run python scripts/atavismo/s26_relicto_comparativo.py
    wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s26_alinhar.sh
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import requests

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
CM = {"tool": "bio-stack", "email": "limatheus458@gmail.com"}

OUT = Path(__file__).resolve().parents[2] / "runs" / "atavismo" / "lncrna"

S = requests.Session()
S.headers["User-Agent"] = "bio-stack/s26"

# taxid, rotulo curto, estado dentario, origem da perda
ESPECIES = [
    (9031, "Ggal", "Gallus gallus", "sem dente", "Neornithes"),
    (60713, "Mund", "Melopsittacus undulatus", "sem dente", "Neornithes"),
    (2489341, "Shab", "Strigops habroptila", "sem dente", "Neornithes"),
    (8478, "Cpic", "Chrysemys picta", "sem dente", "Testudines (independente)"),
    (8496, "Amis", "Alligator mississippiensis", "COM dente", "—"),
    (94835, "Ggan", "Gavialis gangeticus", "COM dente", "—"),
    (28377, "Acar", "Anolis carolinensis", "COM dente", "—"),
]

PAINEL = ["ENAM", "AMBN", "AMTN", "AMELX", "MMP20", "DSPP", "ODAM", "ODAPH", "ACP4"]

FOLGA = 60_000  # bp de folga de cada lado do intervalo NRG1..LPL


def gene_pos(sym: str, taxon: int) -> dict | None:
    """Primeira localizacao genomica anotada do gene. None se nao anotado."""
    try:
        d = S.get(f"{API}/gene/symbol/{sym}/taxon/{taxon}",
                  params={"returned_content": "COMPLETE"}, timeout=90).json()
    except (requests.RequestException, json.JSONDecodeError):
        return None
    for rep in d.get("reports", []):
        g = rep.get("gene", {})
        if g.get("symbol", "").upper() != sym.upper():
            continue
        for a in g.get("annotations", []):
            for gl in a.get("genomic_locations", []):
                gr = gl.get("genomic_range", {})
                return {
                    "asm": a.get("assembly_name"),
                    "acc": gl.get("genomic_accession_version"),
                    "seq": gl.get("sequence_name"),
                    "ini": int(gr.get("begin", 0)),
                    "fim": int(gr.get("end", 0)),
                }
    return None


def fatia(acesso: str, ini: int, fim: int, rotulo: str) -> str | None:
    r = S.get(f"{EU}/efetch.fcgi",
              params={**CM, "db": "nuccore", "id": acesso, "rettype": "fasta",
                      "retmode": "text", "seq_start": max(1, ini), "seq_stop": fim},
              timeout=300)
    if r.status_code != 200 or not r.text.strip().startswith(">"):
        return None
    corpo = r.text.strip().split("\n", 1)[1]
    return f">{rotulo} {acesso}:{ini}-{fim}\n{corpo}\n"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    resumo: list[dict] = []

    print("=" * 78)
    print("PARTE 1 — o painel de nove genes esta ANOTADO em cada especie?")
    print("=" * 78)
    print(f"{'especie':28s} {'estado':10s} " + " ".join(f"{g:6s}" for g in PAINEL))
    anotacao: dict[str, dict[str, str]] = {}
    for taxon, rot, nome, estado, _origem in ESPECIES:
        linha = {}
        marcas = []
        for g in PAINEL:
            pos = gene_pos(g, taxon)
            linha[g] = f"{pos['seq']}" if pos else "-"
            marcas.append(f"{(pos['seq'] if pos else '-'):6s}")
            time.sleep(0.25)
        anotacao[rot] = linha
        print(f"{nome:28s} {estado:10s} " + " ".join(marcas))
    (OUT / "anotacao_painel.json").write_text(
        json.dumps(anotacao, indent=2, ensure_ascii=False), encoding="utf-8")

    print()
    print("=" * 78)
    print("PARTE 2 — intervalo sintenico NRG1..LPL por especie")
    print("=" * 78)
    for taxon, rot, nome, estado, origem in ESPECIES:
        nrg1 = gene_pos("NRG1", taxon)
        time.sleep(0.3)
        lpl = gene_pos("LPL", taxon)
        time.sleep(0.3)

        if not (nrg1 and lpl):
            falta = "NRG1" if not nrg1 else "LPL"
            print(f"  {nome:28s} ancora {falta} nao anotada — PULADO")
            resumo.append({"rotulo": rot, "nome": nome, "estado": estado,
                           "origem": origem, "situacao": f"sem ancora {falta}"})
            continue

        if nrg1["seq"] != lpl["seq"]:
            print(f"  {nome:28s} ancoras em sequencias diferentes "
                  f"(NRG1={nrg1['seq']}, LPL={lpl['seq']}) — SINTENIA QUEBRADA")
            resumo.append({"rotulo": rot, "nome": nome, "estado": estado,
                           "origem": origem, "situacao": "sintenia quebrada",
                           "nrg1_seq": nrg1["seq"], "lpl_seq": lpl["seq"]})
            continue

        lo = min(nrg1["ini"], nrg1["fim"], lpl["ini"], lpl["fim"])
        hi = max(nrg1["ini"], nrg1["fim"], lpl["ini"], lpl["fim"])
        # vao intergenico entre as duas ancoras
        vao = max(lpl["ini"], lpl["fim"]) - min(nrg1["ini"], nrg1["fim"])
        vao_limpo = (min(lpl["ini"], lpl["fim"]) - max(nrg1["ini"], nrg1["fim"]))

        arq = OUT / f"intervalo_{rot}.fna"
        fasta = fatia(lpl["acc"], lo - FOLGA, hi + FOLGA, f"{rot}_NRG1_LPL")
        if not fasta:
            print(f"  {nome:28s} efetch falhou para {lpl['acc']} — PULADO")
            continue
        arq.write_text(fasta, encoding="utf-8")
        bp = len(fasta.split("\n", 1)[1].replace("\n", ""))

        print(f"  {nome:28s} {estado:10s} seq={lpl['seq']:>5s} "
              f"NRG1..LPL={vao_limpo:>9,} bp   janela={bp:>9,} bp   {arq.name}")
        resumo.append({
            "rotulo": rot, "nome": nome, "estado": estado, "origem": origem,
            "situacao": "ok", "assembly": lpl["asm"], "sequencia": lpl["seq"],
            "acesso": lpl["acc"], "janela_ini": lo - FOLGA, "janela_fim": hi + FOLGA,
            "vao_nrg1_lpl": vao_limpo, "span_total": vao, "arquivo": arq.name,
        })

    (OUT / "s26_resumo.json").write_text(
        json.dumps(resumo, indent=2, ensure_ascii=False), encoding="utf-8")

    prontos = [r for r in resumo if r.get("situacao") == "ok"]
    print(f"\n{len(prontos)} de {len(ESPECIES)} especies com intervalo extraido.")
    print("agora o alinhamento, no WSL:")
    print("  wsl -d Ubuntu-24.04 -- bash -l "
          "/mnt/c/Users/Manec/bio-stack/scripts/atavismo/s26_alinhar.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
