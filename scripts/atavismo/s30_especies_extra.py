"""Mais especies no intervalo NRG1..LPL (C4 de docs/METODO-RELIQUIAS.md).

O s29 levou ENAM e AMBN ao nivel 3 (lesao identificavel) no kakapo e na
tartaruga, mas o nivel 4 — lesao COMPARTILHADA por linhagens que divergiram
depois da perda, o que data a perda — nao se testa: no frango a reliquia esta
degradada demais para o miniprot, e o unico par com lesao dos dois lados e ave x
tartaruga, perdas independentes.

Este script baixa o intervalo NRG1..LPL de mais especies, com o mesmo codigo do
s26 (importado, nao copiado):

    PALEOGNATA  Dromaius novaehollandiae — lesao compartilhada paleognata x
                neognata data a perda no ancestral de Neornithes
    NEOGNATAS   Anas platyrhynchos (Galloanserae, irma de Gallus),
                Taeniopygia guttata, Aquila chrysaetos (Neoaves)
    TARTARUGAS  Chelonia mydas, Gopherus evgoodei, Pelodiscus sinensis —
                compartilhamento entre tartarugas data a perda delas
    COM DENTE   Crocodylus porosus, Python bivittatus, Podarcis muralis —
                mais controles negativos de lesao

Escreve intervalo_<rot>.fna em runs/atavismo/lncrna/, onde o s29_lesoes.sh ja
procura. Depois:

    uv run python scripts/atavismo/s30_especies_extra.py
    wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s29_lesoes.sh
    uv run python scripts/atavismo/s29_lesoes.py
"""

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

AQUI = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s26", AQUI / "s26_relicto_comparativo.py")
if _spec is None or _spec.loader is None:
    raise SystemExit("nao carreguei s26_relicto_comparativo.py")
s26 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s26)

# taxid, rotulo, nome, dente, grupo
EXTRA = [
    (8790, "Dnov", "Dromaius novaehollandiae", "sem dente", "Palaeognathae"),
    (8801, "Scam", "Struthio camelus", "sem dente", "Palaeognathae"),
    (308060, "Arow", "Apteryx rowi", "sem dente", "Palaeognathae"),
    (8787, "Ccas", "Casuarius casuarius", "sem dente", "Palaeognathae"),
    (94827, "Tgua", "Tinamus guttatus", "sem dente", "Palaeognathae"),
    (8839, "Apla", "Anas platyrhynchos", "sem dente", "Galloanserae"),
    (59729, "Tgut", "Taeniopygia guttata", "sem dente", "Neoaves"),
    (8962, "Achr", "Aquila chrysaetos", "sem dente", "Neoaves"),
    (8469, "Cmyd", "Chelonia mydas", "sem dente", "Testudines"),
    (1825980, "Gevg", "Gopherus evgoodei", "sem dente", "Testudines"),
    (13735, "Psin", "Pelodiscus sinensis", "sem dente", "Testudines"),
    (8502, "Cpor", "Crocodylus porosus", "COM dente", "Crocodylia"),
    (176946, "Pbiv", "Python bivittatus", "COM dente", "Squamata"),
    (64176, "Pmur", "Podarcis muralis", "COM dente", "Squamata"),
]
TETO = 3_000_000
# Uma ancora so, quando a outra nao esta anotada (a ema nao tem LPL anotado).
# Nas aves o vao NRG1..LPL tem 105-300 kb; +-800 kb em volta da ancora presente
# cobre o vao dos dois lados, seja qual for a fita.
JANELA_UMA = 800_000


def main() -> int:
    resumo = []
    for taxon, rot, nome, dente, grupo in EXTRA:
        nrg1 = s26.gene_pos("NRG1", taxon)
        time.sleep(0.3)
        lpl = s26.gene_pos("LPL", taxon)
        time.sleep(0.3)
        base = {"rotulo": rot, "nome": nome, "dente": dente, "grupo": grupo}
        nota = "par NRG1/LPL"
        if nrg1 and lpl and nrg1["acc"] == lpl["acc"]:
            lo = min(nrg1["ini"], nrg1["fim"], lpl["ini"], lpl["fim"]) - s26.FOLGA
            hi = max(nrg1["ini"], nrg1["fim"], lpl["ini"], lpl["fim"]) + s26.FOLGA
            acc, seq = lpl["acc"], lpl["seq"]
        elif nrg1 or lpl:
            p = lpl or nrg1
            assert p is not None
            meio = (p["ini"] + p["fim"]) // 2
            lo, hi, acc, seq = meio - JANELA_UMA, meio + JANELA_UMA, p["acc"], p["seq"]
            nota = f"uma ancora ({'LPL' if lpl else 'NRG1'})"
        else:
            print(f"  {nome:28s} sem NRG1 nem LPL anotados — pulado")
            resumo.append({**base, "situacao": "sem ancora"})
            continue
        if hi - lo > TETO:
            print(f"  {nome:28s} NRG1..LPL com {hi - lo:,} bp — acima do teto, pulado")
            resumo.append({**base, "situacao": f"intervalo de {hi - lo:,} bp"})
            continue
        arq = s26.OUT / f"intervalo_{rot}.fna"
        if not arq.exists():
            fasta = s26.fatia(acc, max(1, lo), hi, f"{rot}_NRG1_LPL")
            if not fasta:
                resumo.append({**base, "situacao": "efetch falhou"})
                continue
            arq.write_text(fasta, encoding="utf-8")
        bp = len(arq.read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", ""))
        print(f"  {nome:28s} {dente:9s} seq={seq:>6} {bp:>9,} bp  {nota}")
        resumo.append({**base, "situacao": "ok", "ancoras": nota, "acesso": acc, "sequencia": seq,
                       "janela_ini": max(1, lo), "janela_fim": hi, "bp": bp})
    (s26.OUT / "s30_resumo.json").write_text(json.dumps(resumo, indent=2, ensure_ascii=False),
                                             encoding="utf-8")
    print(f"\n{sum(r['situacao'] == 'ok' for r in resumo)} de {len(EXTRA)} especies com intervalo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
