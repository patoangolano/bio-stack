"""Le os alinhamentos do miniprot (s29_lesoes.sh): lesoes por gene e especie, e
o teste de lesao compartilhada (nivel 4 da escada da skill atavismo).

LESOES, na coordenada da proteina do jacare:
  frameshift  operacao F ou G no cg (1-2 nt sobrando ou faltando no genoma)
  stop        substituicao por codon de parada no cs (*taa, *tag, *tga)

ESCOLHA DO ALINHAMENTO: o de maior AS. No MMP20 so contam alinhamentos que
encostam num trecho RECIPROCO do s27_reciproco.tsv — o miniprot pode pousar nas
colagenases do cluster, e o reciproco ja decidiu quais trechos sao MMP20.

CONTROLE NEGATIVO DE LESAO: as tres especies com dente. Gene funcional nao tem
frameshift nem stop; lesao ali e artefato de alinhamento ou de montagem, e mede
a taxa de lesao falsa do metodo.

LESAO COMPARTILHADA: mesmo tipo, posicao a <= TOL aa de distancia. O nulo e o
numero de coincidencias esperado se as lesoes caissem ao acaso na regiao que as
duas especies alinham em comum:  E = sum_tipo n1 * n2 * (2*TOL + 1) / L.
Par de aves (perda unica, Meredith 2014): compartilhar acima do nulo data a
lesao antes da separacao das duas linhagens. Ave x tartaruga (perdas
independentes): a expectativa e ficar no nulo — e o controle do teste.

    uv run python scripts/atavismo/s29_lesoes.py
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2] / "runs" / "atavismo"
LES = RAIZ / "lesoes"
RECIPROCO = RAIZ / "painel" / "s27_reciproco.tsv"
GENES = ["ENAM", "AMBN", "AMTN", "MMP20", "AMELX", "DSPP", "ACP4"]
# rotulo -> (dente, grupo). As seis do s26/s27, mais as do s30_especies_extra.py.
BASE = {"Amis": ("COM", "Crocodylia"), "Ggan": ("COM", "Crocodylia"), "Acar": ("COM", "Squamata"),
        "Cpic": ("sem", "Testudines"), "Ggal": ("sem", "Galloanserae"), "Shab": ("sem", "Neoaves")}
FILTRA_RECIPROCO = {"MMP20"}
TOL = 2
STOPS = {"taa", "tag", "tga"}


@dataclass
class Aln:
    qlen: int
    qini: int
    qfim: int
    sini: int
    sfim: int
    escore: int
    fs: int
    st: int
    identidade: float
    lesoes: list[tuple[str, int]] = field(default_factory=list)

    @property
    def cobertura(self) -> float:
        return (self.qfim - self.qini) / self.qlen


def le_paf(caminho: Path) -> list[Aln]:
    alns = []
    if not caminho.exists():
        return alns
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        c = linha.split("\t")
        if len(c) < 12:
            continue
        tags = {t[:2]: t[5:] for t in c[12:]}
        cg, cs = tags.get("cg", ""), tags.get("cs", "")
        lesoes: list[tuple[str, int]] = []
        # frameshift: anda pelo cg; M e I consomem a proteina
        pos = int(c[2])
        for n, op in re.findall(r"(\d+)([A-Z])", cg):
            if op in "MI":
                pos += int(n)
            elif op in "FG":
                lesoes.append(("fs", pos))
        # stop e identidade: anda pelo cs
        pos, iguais, trocas = int(c[2]), 0, 0
        for tok in re.findall(r":\d+|\*[a-z]{3}[A-Z*]|\+[A-Z*]+|-[a-z]+|~[a-z]+\d+[a-z]+", cs):
            if tok.startswith(":"):
                iguais += int(tok[1:])
                pos += int(tok[1:])
            elif tok.startswith("*"):
                if tok[1:4] in STOPS:
                    lesoes.append(("stop", pos))
                trocas += 1
                pos += 1
            elif tok.startswith("+"):
                pos += len(tok) - 1
        alns.append(Aln(int(c[1]), int(c[2]), int(c[3]), int(c[7]), int(c[8]), int(tags.get("AS", 0)),
                        int(tags.get("fs", 0)), int(tags.get("st", 0)),
                        iguais / (iguais + trocas) if iguais + trocas else 0.0, lesoes))
    return alns


def reciprocos() -> dict[tuple[str, str], list[tuple[int, int]]]:
    trechos: dict[tuple[str, str], list[tuple[int, int]]] = {}
    if RECIPROCO.exists():
        with RECIPROCO.open(encoding="utf-8") as f:
            for c in csv.reader(f, delimiter="\t"):
                if len(c) >= 6 and c[5] == "RECIPROCO":
                    trechos.setdefault((c[0], c[1]), []).append((int(c[2]), int(c[3])))
    return trechos


def especies() -> dict[str, tuple[str, str]]:
    esp = dict(BASE)
    extra = RAIZ / "lncrna" / "s30_resumo.json"
    if extra.exists():
        for r in json.loads(extra.read_text(encoding="utf-8")):
            if r.get("situacao") == "ok":
                esp[r["rotulo"]] = ("COM" if r["dente"].startswith("COM") else "sem", r["grupo"])
    return esp


AVES = {"Palaeognathae", "Galloanserae", "Neoaves"}


def tipo_par(g1: str, g2: str) -> str:
    """O que um par de especies sem dente testa."""
    if g1 in AVES and g2 in AVES:
        if "Palaeognathae" in (g1, g2) and g1 != g2:
            return "paleognata x neognata"  # compartilhada = antes de Neornithes
        return "aves, mesmo grupo" if g1 == g2 else "aves: Galloanserae x Neoaves"
    if g1 == g2 == "Testudines":
        return "tartarugas"
    return "ave x tartaruga (indep.)"  # perdas independentes: o controle


def cauda_poisson(obs: int, lam: float) -> float:
    """P(X >= obs) para X ~ Poisson(lam)."""
    return 1.0 - sum(math.exp(-lam) * lam**k / math.factorial(k) for k in range(obs))


def main() -> int:
    rec = reciprocos()
    esp = especies()
    DENTE = {k: v[0] for k, v in esp.items()}  # noqa: N806 — tabela, nao variavel
    melhor: dict[tuple[str, str], Aln] = {}
    linhas = []
    print(f"{'gene':6s} {'esp':5s} {'dente':5s} {'aln':>4s} {'proteina':>13s} {'cob':>5s} {'ident':>6s} "
          f"{'fs':>3s} {'stop':>4s}  leitura")
    for g in GENES:
        for e in esp:
            if not (LES / f"mp_{g}_{e}.paf").exists():
                continue
            alns = le_paf(LES / f"mp_{g}_{e}.paf")
            if g in FILTRA_RECIPROCO:
                regioes = rec.get((g, e), [])
                alns = [a for a in alns
                        if any(min(a.sfim, r[1]) >= max(a.sini, r[0]) for r in regioes)]
            if not alns:
                print(f"{g:6s} {e:5s} {DENTE[e]:5s} {0:4d} {'—':>13s}")
                linhas.append({"gene": g, "especie": e, "dente": DENTE[e], "alinhou": False})
                continue
            a = max(alns, key=lambda x: x.escore)
            melhor[(g, e)] = a
            lesoes = a.fs + a.st
            if DENTE[e] == "COM":
                leitura = "intacto" if lesoes == 0 else "LESAO EM ESPECIE COM DENTE"
            else:
                leitura = "pseudogene (nivel 3)" if lesoes else "sem lesao no trecho alinhado"
            print(f"{g:6s} {e:5s} {DENTE[e]:5s} {len(alns):4d} {a.qini + 1:>5d}-{a.qfim:<5d}/{a.qlen:<3d}"
                  f"{a.cobertura:5.0%} {a.identidade:6.2f} {a.fs:3d} {a.st:4d}  {leitura}")
            linhas.append({"gene": g, "especie": e, "dente": DENTE[e], "alinhou": True,
                           "q_ini": a.qini + 1, "q_fim": a.qfim, "q_len": a.qlen,
                           "cobertura": round(a.cobertura, 3), "identidade": round(a.identidade, 3),
                           "frameshifts": a.fs, "stops": a.st, "escore": a.escore,
                           "s_ini": a.sini + 1, "s_fim": a.sfim,
                           "lesoes": ";".join(f"{t}@{p}" for t, p in a.lesoes)})

    print(f"\n=== lesoes compartilhadas (mesmo tipo, <= {TOL} aa) vs nulo de posicao aleatoria ===")
    print(f"{'gene':6s} {'par':10s} {'tipo de par':30s} {'L':>4s} {'n1':>3s} {'n2':>3s} {'obs':>4s} {'esp':>6s}  posicoes")
    pares = []
    sem_dente = [e for e in esp if esp[e][0] == "sem"]
    for g in GENES:
        for e1, e2 in combinations(sem_dente, 2):
            a1, a2 = melhor.get((g, e1)), melhor.get((g, e2))
            if not (a1 and a2):
                continue
            lo, hi = max(a1.qini, a2.qini), min(a1.qfim, a2.qfim)
            if hi - lo < 30:
                continue
            l1 = [x for x in a1.lesoes if lo <= x[1] <= hi]
            l2 = [x for x in a2.lesoes if lo <= x[1] <= hi]
            obs = [(t, p1, p2) for t, p1 in l1 for t2, p2 in l2 if t == t2 and abs(p1 - p2) <= TOL]
            esperado = sum(sum(1 for x in l1 if x[0] == t) * sum(1 for x in l2 if x[0] == t)
                           for t in ("fs", "stop")) * (2 * TOL + 1) / (hi - lo)
            tipo = tipo_par(esp[e1][1], esp[e2][1])
            print(f"{g:6s} {e1}-{e2:5s} {tipo:30s} {hi - lo:4d} {len(l1):3d} {len(l2):3d} {len(obs):4d} "
                  f"{esperado:6.2f}  {', '.join(f'{t}@{p1}/{p2}' for t, p1, p2 in obs) or '-'}")
            pares.append({"gene": g, "par": f"{e1}-{e2}", "tipo": tipo, "L": hi - lo, "n1": len(l1),
                          "n2": len(l2), "observado": len(obs), "esperado": round(esperado, 2),
                          "posicoes": ";".join(f"{t}@{p1}/{p2}" for t, p1, p2 in obs)})

    # Agregado por tipo de par. Descritivo: pares que dividem especie nao sao
    # independentes, e a cauda de Poisson so diz quao longe do acaso esta.
    print("\n=== agregado por tipo de par (descritivo; pares nao independentes) ===")
    for tipo in sorted({p["tipo"] for p in pares}):
        sel = [p for p in pares if p["tipo"] == tipo]
        o, e = sum(p["observado"] for p in sel), sum(p["esperado"] for p in sel)
        print(f"  {tipo:30s} {len(sel):3d} pares  observado {o:3d}  esperado {e:6.2f}  "
              f"razao {o / e if e else float('nan'):5.1f}  P(X>={o}) = {cauda_poisson(o, e):.2g}")

    for nome, dados in (("s29_lesoes.tsv", linhas), ("s29_compartilhadas.tsv", pares)):
        if dados:
            campos = list(dict.fromkeys(k for d in dados for k in d))
            with (LES / nome).open("w", encoding="utf-8", newline="") as f:
                w = csv.DictWriter(f, fieldnames=campos, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(dados)
    print(f"\ntabelas: {LES / 's29_lesoes.tsv'}  |  {LES / 's29_compartilhadas.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
