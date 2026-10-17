"""Indice de retencao do s27, com competicao de paralogos.

Le os s27_<ALVO>_<rot>.tsv do s27_alinhar.sh e decide, trecho a trecho do
genoma, quem fica com cada HSP. A regra:

    um HSP do ALVO e descartado se outra consulta do mesmo intervalo tiver um HSP
    que cubra pelo menos METADE do trecho dele com bitscore maior.

E o equivalente local do melhor-hit-reciproco do s25: se o trecho e mais
parecido com uma colagenase do jacare do que com o MMP20 do jacare, ele e
colagenase, nao reliquia de MMP20.

Indices, por especie, com as ancoras do intervalo como denominador:

    indice_melhor = bits do melhor HSP sobrevivente do alvo / media dos melhores
                    HSPs das ancoras                       (o mesmo do s26)
    indice_soma   = soma dos HSPs sobreviventes do alvo que nao se sobrepoem /
                    media da mesma soma nas ancoras        (pendencia do s26:
                    o melhor HSP sozinho subestima gene fragmentado em exons)

Ancora com zero bits na especie invalida a linha: o intervalo esta errado ou o
metodo nao enxerga ali, e nenhum zero do alvo significa perda.

    uv run python scripts/atavismo/s27_indice.py
"""

from __future__ import annotations

import csv
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path

import requests

PASTA = Path(__file__).resolve().parents[2] / "runs" / "atavismo" / "painel"
ESTADO = {"Amis": "COM", "Ggan": "COM", "Acar": "COM", "Cpic": "sem", "Ggal": "sem",
          "Mund": "sem", "Shab": "sem"}
ORDEM = ["Amis", "Ggan", "Acar", "Cpic", "Ggal", "Mund", "Shab"]
SOBREPOSICAO = 0.5


@dataclass
class Hsp:
    q: str
    pident: float
    s_ini: int
    s_fim: int
    evalue: float
    bits: float

    @property
    def tam(self) -> int:
        return self.s_fim - self.s_ini + 1


def le_tsv(caminho: Path) -> list[Hsp]:
    hsps = []
    with caminho.open(encoding="utf-8") as f:
        for c in csv.reader(f, delimiter="\t"):
            ss, se = int(c[6]), int(c[7])
            hsps.append(Hsp(c[0], float(c[2]), min(ss, se), max(ss, se), float(c[8]), float(c[9])))
    return hsps


def sobrepoe(a: Hsp, b: Hsp) -> bool:
    inter = min(a.s_fim, b.s_fim) - max(a.s_ini, b.s_ini) + 1
    return inter >= SOBREPOSICAO * min(a.tam, b.tam)


def soma_sem_sobreposicao(hsps: list[Hsp]) -> float:
    """Soma gulosa: do maior bitscore para baixo, so entra quem nao encosta nos escolhidos."""
    escolhidos: list[Hsp] = []
    for h in sorted(hsps, key=lambda x: -x.bits):
        if not any(min(h.s_fim, e.s_fim) >= max(h.s_ini, e.s_ini) for e in escolhidos):
            escolhidos.append(h)
    return sum(e.bits for e in escolhidos)


def descricoes(acessos: list[str]) -> dict[str, str]:
    """Titulo da proteina no NCBI, para nomear o competidor que venceu."""
    if not acessos:
        return {}
    params = {"db": "protein", "id": ",".join(acessos), "retmode": "json",
              "tool": "bio-stack", "email": "limatheus458@gmail.com"}
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    try:
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                         params=params, timeout=120).json()["result"]
    except (requests.RequestException, KeyError, ValueError):
        return {}
    return {v["accessionversion"]: v["title"].split(" [")[0]
            for k, v in r.items() if k != "uids"}


def main() -> int:
    resumo = json.loads((PASTA / "s27_resumo.json").read_text(encoding="utf-8"))
    linhas_saida = []
    sobreviventes: list[dict] = []  # HSPs do alvo que venceram a competicao -> s27_reciproco.sh
    produtos: dict[str, str] = {}  # alvo -> titulo da proteina no NCBI -> s27_reciproco.sh
    for painel in resumo:
        alvo = painel["alvo"]
        faa = PASTA / f"painel_{alvo}_proteinas.faa"
        if not faa.exists():
            continue
        acesso_de = {}
        for linha in faa.read_text(encoding="utf-8").splitlines():
            if linha.startswith(">"):
                q, *resto = linha[1:].split()
                acesso_de[q] = resto[0] if resto else ""
        desc = descricoes([a for a in acesso_de.values() if a])
        q_alvo = f"{alvo}_Amis"
        # Produto do alvo pelo titulo, nao pelo acesso: a tabela de features cita
        # versoes de acesso (XP_014463267.1) que o proteoma atual ja substituiu.
        produtos[alvo] = re.sub(r" isoform X\d+", "", desc.get(acesso_de.get(q_alvo, ""), ""))

        print(f"\n=== {alvo}  (ancoras {painel['ancora_esq']} / {painel['ancora_dir']}, "
              f"{len(painel['competidores'])} competidores) ===")
        print(f"{'esp':5s} {'dente':5s} {'alvo_melhor':>11s} {'alvo_soma':>9s} {'HSPs':>5s} "
              f"{'perdidos':>8s} {'ancoras(melhor)':>18s} {'idx_melhor':>10s} {'idx_soma':>8s} "
              f"{'LPL(-)':>6s}  observacao")
        intervalos = {i["rotulo"]: i for i in painel["intervalos"]}
        for rot in ORDEM:
            iv = intervalos.get(rot)
            tsv = PASTA / f"s27_{alvo}_{rot}.tsv"
            if not iv or not iv.get("arquivo") or not tsv.exists():
                motivo = iv["situacao"] if iv else "sem intervalo"
                print(f"{rot:5s} {ESTADO[rot]:5s} {'—':>11s}  {motivo}")
                linhas_saida.append({"alvo": alvo, "especie": rot, "dente": ESTADO[rot], "situacao": motivo})
                continue
            hsps = le_tsv(tsv)
            do_alvo = [h for h in hsps if h.q == q_alvo]
            outros = [h for h in hsps if h.q != q_alvo and h.q != "LPL_Amis"]
            vivos, venceram = [], {}
            for h in do_alvo:
                rival = max((o for o in outros if o.bits > h.bits and sobrepoe(h, o)),
                            key=lambda o: o.bits, default=None)
                if rival is None:
                    vivos.append(h)
                else:
                    venceram[rival.q] = venceram.get(rival.q, 0) + 1
            melhor = max((h.bits for h in vivos), default=0.0)
            soma = soma_sem_sobreposicao(vivos)

            anc_melhor, anc_soma = [], []
            for a in iv["ancoras"]:
                ha = [h for h in hsps if h.q == f"{a}_Amis"]
                anc_melhor.append(max((h.bits for h in ha), default=0.0))
                anc_soma.append(soma_sem_sobreposicao(ha))
            den_m = sum(anc_melhor) / len(anc_melhor) if anc_melhor else 0.0
            den_s = sum(anc_soma) / len(anc_soma) if anc_soma else 0.0
            controle_ok = bool(anc_melhor) and min(anc_melhor) > 50
            idx_m = melhor / den_m if den_m else float("nan")
            idx_s = soma / den_s if den_s else float("nan")
            lpl = max((h.bits for h in hsps if h.q == "LPL_Amis"), default=0.0)

            obs = []
            if not controle_ok:
                obs.append("CONTROLE FALHOU (ancora < 50 bits)")
            if iv["situacao"] != "ok":
                obs.append(iv["situacao"])
            for h in vivos:
                sobreviventes.append({"alvo": alvo, "especie": rot, "arquivo": iv["arquivo"],
                                      "s_ini": h.s_ini, "s_fim": h.s_fim, "bits": h.bits,
                                      "evalue": h.evalue, "lpl_bits": lpl})
            if venceram:
                top = max(venceram, key=lambda k: venceram[k])
                nome = desc.get(acesso_de.get(top, ""), "")
                obs.append(f"{sum(venceram.values())} HSP(s) do alvo cedidos; maior rival {top} ({nome})")
            anc_txt = "/".join(f"{x:.0f}" for x in anc_melhor)
            print(f"{rot:5s} {ESTADO[rot]:5s} {melhor:11.1f} {soma:9.1f} {len(do_alvo):5d} "
                  f"{len(do_alvo) - len(vivos):8d} {anc_txt:>18s} {idx_m:10.3f} {idx_s:8.3f} "
                  f"{lpl:6.1f}  {'; '.join(obs)}")
            linhas_saida.append({
                "alvo": alvo, "especie": rot, "dente": ESTADO[rot], "situacao": iv["situacao"],
                "alvo_melhor_bits": melhor, "alvo_soma_bits": round(soma, 1),
                "hsps_alvo": len(do_alvo), "hsps_cedidos": len(do_alvo) - len(vivos),
                "ancoras": "/".join(iv["ancoras"]), "ancoras_melhor_bits": anc_txt,
                "indice_melhor": round(idx_m, 3), "indice_soma": round(idx_s, 3),
                "lpl_negativo_bits": lpl, "controle_ok": controle_ok, "obs": "; ".join(obs),
            })

    campos = sorted({k for linha in linhas_saida for k in linha}, key=lambda k: (k != "alvo", k != "especie", k))
    with (PASTA / "s27_indice.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos, delimiter="\t")
        w.writeheader()
        w.writerows(linhas_saida)
    with (PASTA / "s27_sobreviventes.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sobreviventes[0]) if sobreviventes else ["alvo"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(sobreviventes)
    # newline="\n": quem le e o bash do WSL, e um \r no fim do titulo faz toda
    # comparacao falhar (medido: 66 de 66 vereditos errados).
    (PASTA / "s27_produtos.tsv").write_text(
        "".join(f"{a}\t{t}\n" for a, t in produtos.items()), encoding="utf-8", newline="\n")
    print(f"\ntabela: {PASTA / 's27_indice.tsv'}")
    print(f"HSPs sobreviventes (para o reciproco): {len(sobreviventes)} em {PASTA / 's27_sobreviventes.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
