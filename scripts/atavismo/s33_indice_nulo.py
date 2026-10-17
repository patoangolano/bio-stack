"""Indice de retencao do s26 com nulo por janela e segunda consulta (§8 itens 2 e 4
de docs/METODO-RELIQUIAS.md, critica externa itens a e b).

O problema medido: no s26, tres genes a ~25 bits (ruido) divididos por ~170 dao
indice ~0,4 sem reliquia nenhuma. Aqui cada bitscore e comparado com o NULO da
propria janela: 20 versoes embaralhadas de cada proteina (mesma composicao,
mesmo comprimento, ordem destruida) contra o mesmo intervalo, mesmo tBLASTn. So
conta o EXCESSO sobre o percentil 95 do nulo.

Segunda consulta: as proteinas do Anolis para ENAM, AMBN, AMTN e LPL. Com o
jacare como unica consulta, a distancia a consulta e a especie com dente mais
proxima se confundem; com o Anolis, os crocodilianos passam a ser "distantes".

    uv run python scripts/atavismo/s33_indice_nulo.py preparar
    wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s33_alinhar.sh
    uv run python scripts/atavismo/s33_indice_nulo.py analisar
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parents[2]
LNC = RAIZ / "runs" / "atavismo" / "lncrna"
OUT = RAIZ / "runs" / "atavismo" / "indice_nulo"
GENES = ["ENAM", "AMBN", "AMTN", "LPL"]
N_EMBARALHADAS = 20
SEMENTE = 20261017

_spec = importlib.util.spec_from_file_location("s27", Path(__file__).parent / "s27_painel_intervalos.py")
if _spec is None or _spec.loader is None:
    raise SystemExit("nao carreguei s27_painel_intervalos.py")
s27 = importlib.util.module_from_spec(_spec)
sys.modules["s27"] = s27  # as dataclasses do s27 procuram o proprio modulo aqui
_spec.loader.exec_module(s27)


def le_faa(caminho: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    nome = ""
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        if linha.startswith(">"):
            nome = linha[1:].split()[0]
            seqs[nome] = []
        elif nome:
            seqs[nome].append(linha.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def proteinas_anolis() -> dict[str, str]:
    """Primeira proteina de cada gene no Anolis, pela tabela de features (como no s27)."""
    acessos = {}
    for g in GENES:
        p = s27.gene_pos(g, 28377)
        if not p:
            print(f"  Anolis: {g} nao anotado")
            continue
        genes = s27.tabela_features(p["acc"], p["ini"], p["fim"])
        # no Anolis a tabela de features traz o simbolo em minuscula ("enam")
        gene = next((x for x in genes if x.sym.upper() == g and x.proteina), None)
        if gene and gene.proteina:
            acessos[g] = gene.proteina
        else:
            print(f"  Anolis: {g} sem proteina na tabela de features")
    if len(acessos) < len(GENES):
        raise SystemExit("faltou proteina do Anolis — sem ela nao ha segunda consulta")
    texto = s27.proteinas(acessos)  # cabecalho >SIMBOLO_Amis — renomeado abaixo
    return {k.replace("_Amis", "_Acar"): v for k, v in le_faa_texto(texto).items()}


def le_faa_texto(texto: str) -> dict[str, str]:
    tmp = OUT / "_tmp.faa"
    tmp.write_text(texto, encoding="utf-8")
    seqs = le_faa(tmp)
    tmp.unlink()
    return seqs


def preparar() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEMENTE)  # gerador proprio; nunca np.random.seed
    jac = {k: v for k, v in le_faa(LNC / "jacare_proteinas.faa").items() if k.split("_")[0] in GENES}
    ano = proteinas_anolis()
    consultas = {**jac, **ano}
    linhas = []
    for nome, seq in consultas.items():
        linhas.append(f">{nome}\n{seq}\n")
        for i in range(N_EMBARALHADAS):
            emb = "".join(rng.permutation(list(seq)))
            linhas.append(f">{nome}~emb{i:02d}\n{emb}\n")
    (OUT / "consultas.faa").write_text("".join(linhas), encoding="utf-8", newline="\n")
    print(f"consultas: {len(consultas)} reais ({', '.join(sorted(consultas))}) "
          f"+ {N_EMBARALHADAS} embaralhadas de cada (semente {SEMENTE})")
    return 0


def analisar() -> int:
    esp = {"Amis": "COM", "Ggan": "COM", "Acar": "COM", "Cpic": "sem", "Ggal": "sem", "Shab": "sem"}
    extra = LNC / "s30_resumo.json"
    if extra.exists():
        for r in json.loads(extra.read_text(encoding="utf-8")):
            if r.get("situacao") == "ok":
                esp[r["rotulo"]] = "COM" if r["dente"].startswith("COM") else "sem"
    saida = []
    for fonte in ("Amis", "Acar"):
        print(f"\n=== consulta: proteinas do {'jacare' if fonte == 'Amis' else 'Anolis'} ===")
        print(f"{'esp':5s} {'dente':5s} " + " ".join(f"{g + ' obs/n95':>16s}" for g in GENES)
              + f" {'indice_bruto':>12s} {'indice_nulo':>11s}")
        for rot, dente in esp.items():
            tsv = OUT / f"s33_{rot}.tsv"
            if not tsv.exists():
                continue
            melhor: dict[str, float] = {}
            with tsv.open(encoding="utf-8") as f:
                for c in csv.reader(f, delimiter="\t"):
                    melhor[c[0]] = max(melhor.get(c[0], 0.0), float(c[9]))
            obs, n95 = {}, {}
            for g in GENES:
                q = f"{g}_{fonte}"
                obs[g] = melhor.get(q, 0.0)
                nulo = [melhor.get(f"{q}~emb{i:02d}", 0.0) for i in range(N_EMBARALHADAS)]
                n95[g] = float(np.percentile(nulo, 95))
            exc = {g: max(0.0, obs[g] - n95[g]) for g in GENES}
            bruto = sum(obs[g] for g in GENES[:3]) / obs["LPL"] if obs["LPL"] else float("nan")
            corr = sum(exc[g] for g in GENES[:3]) / exc["LPL"] if exc["LPL"] else float("nan")
            print(f"{rot:5s} {dente:5s} " + " ".join(f"{obs[g]:7.1f}/{n95[g]:<8.1f}" for g in GENES)
                  + f" {bruto:12.2f} {corr:11.2f}")
            saida.append({"consulta": fonte, "especie": rot, "dente": dente,
                          **{f"{g}_obs": obs[g] for g in GENES}, **{f"{g}_nulo95": n95[g] for g in GENES},
                          "indice_bruto": round(bruto, 3), "indice_nulo": round(corr, 3)})
    with (OUT / "s33_indice.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(saida[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(saida)
    print(f"\ntabela: {OUT / 's33_indice.tsv'}")
    return 0


if __name__ == "__main__":
    modo = sys.argv[1] if len(sys.argv) > 1 else ""
    if modo not in ("preparar", "analisar"):
        raise SystemExit("uso: s33_indice_nulo.py preparar|analisar")
    raise SystemExit(preparar() if modo == "preparar" else analisar())
