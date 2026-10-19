"""
f02d — genes candidatos a locos a partir do relatorio de genes do Gallus gallus.

Filtro (antes de olhar qualquer ortologo):
  * PROTEIN_CODING, com grupo "NCBI Ortholog";
  * simbolo oficial (nao LOC/uncharacterized), cromossomo nomeado (nao W, MT, Un);
  * fora de familias multigenicas e de genes com retrocopias abundantes
    (padroes em EXCLUIR).
Os candidatos saem numa permutacao aleatoria fixa, rng = np.random.default_rng(0):
o f02f percorre essa ordem e fica com os primeiros que passam nos filtros de
ortologia — equivale a sortear entre os aprovados.

    uv run python scripts/rapinantes/filogenomica/f02d_candidatos.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parents[3]
WSL = Path(r"\\wsl.localhost\Ubuntu-24.04\home\manec\rapinantes\locos")
SAIDA = RAIZ / "data" / "proc" / "rapinantes" / "filogenomica"

# Familias multigenicas / com muitos paralogos ou retrocopias em aves.
EXCLUIR = re.compile(
    r"^(LOC|OR\d|OR[0-9A-Z]+\b|TAS[12]R|BF\d|BL[AB]?\d|MHC|TAP\d|BG\d|"
    r"IG[HKLJ]|TR[ABDG][VJC]|CD1|CD8|CD4\b|"
    r"H1|H2A|H2B|H3|H4|HIST|"
    r"KRT|FK\d|"  # queratinas alfa e beta
    r"ZNF|ZSCAN|ZFP|ZKSCAN|ZBTB|TRIM|BTN|"
    r"AVBD|GAL\d|DEFB|CATH|"  # defensinas / catelicidinas
    r"CHIR|LILR|KIR|NKR|"
    r"CCL|CXCL|XCL|IFN|IL\d|"
    r"PCDH[ABG]|OLF|VMO|"
    r"UGT|CYP\d|SULT|GST[AMPT]|"
    r"RPL\d|RPS\d|MRP[LS]\d|EEF1A|"
    r"MUC|TTN|NEB|"
    r"SCGB|SPRR|LCE|S100|"
    r"TLR\d|NLR|NOD\d|"
    r"SNOR|MIR|RNU|"
    r"ACT[ABG]|TUB[AB]|"
    r"HSP|HSPA|UBB|UBC|"
    r"GAPDH)",
    re.IGNORECASE,
)


def main() -> None:
    linhas = (WSL / "galo_genes.jsonl").read_text(encoding="utf-8").splitlines()
    cand: list[tuple[str, str, str]] = []
    motivos: dict[str, int] = {}

    def fora(m: str) -> None:
        motivos[m] = motivos.get(m, 0) + 1

    for linha in linhas:
        g = json.loads(linha)
        if g.get("type") != "PROTEIN_CODING":
            fora("nao codificante")
            continue
        sym = g.get("symbol", "")
        if not any(gg.get("method") == "NCBI Ortholog" for gg in g.get("gene_groups", [])):
            fora("sem grupo de ortologos")
            continue
        if EXCLUIR.match(sym) or "uncharacterized" in g.get("description", "").lower():
            fora("familia multigenica / LOC")
            continue
        crom = g.get("chromosomes") or []
        if not crom or crom[0] in {"W", "MT", "Un"}:
            fora("cromossomo W/MT/Un")
            continue
        cand.append((g["gene_id"], sym, crom[0]))
    rng = np.random.default_rng(0)
    ordem = rng.permutation(len(cand))
    SAIDA.mkdir(parents=True, exist_ok=True)
    with (SAIDA / "candidatos_permutados.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("ordem\tgene_id\tsimbolo\tcromossomo_galo\n")
        for i, k in enumerate(ordem):
            gid, sym, c = cand[k]
            fh.write(f"{i}\t{gid}\t{sym}\t{c}\n")
    print(f"{len(linhas)} genes no relatorio; {len(cand)} candidatos; fora: {motivos}")


if __name__ == "__main__":
    main()
