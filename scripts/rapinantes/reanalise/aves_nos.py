"""
Sintese Aves 1.6 (McTavish et al. 2025, PNAS 122:e2409658122; AvesData,
Tree_versions/Aves_1.6): arvores-fonte que sustentam e que conflitam com os
mesmos nos-chave de opentree_nos.py, lidas do annotations.json da sintese.

    uv run python scripts/rapinantes/reanalise/aves_nos.py

Entrada (baixada do GitHub McTavishLab/AvesData para opentree/aves_1.6/):
labelled_supertree_ottnames.tre (rotulo de no interno = id do no na sintese),
annotations.json (supported_by, conflicts_with, partial_path_of, terminal,
resolves por no) e rank_collection.json (ordem de prioridade das arvores-fonte:
a sintese e "veto reverso" - o no que entra e o da arvore de maior prioridade).
Saida: opentree/aves_1.6/nos_aves16.tsv.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from opentree_nos import NOS  # noqa: E402

RAIZ = Path(__file__).resolve().parents[3]
D = RAIZ / "runs" / "rapinantes" / "reanalise" / "opentree" / "aves_1.6"


class No:
    __slots__ = ("filhos", "nome", "pai")

    def __init__(self) -> None:
        self.filhos: list[No] = []
        self.nome = ""
        self.pai: No | None = None


def ler(s: str) -> No:
    raiz = No()
    pilha = [raiz]
    atual = raiz
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            novo = No()
            novo.pai = atual
            atual.filhos.append(novo)
            pilha.append(novo)
            atual = novo
            i += 1
        elif c == ",":
            pilha.pop()
            novo = No()
            novo.pai = pilha[-1]
            pilha[-1].filhos.append(novo)
            pilha.append(novo)
            atual = novo
            i += 1
        elif c == ")":
            pilha.pop()
            atual = pilha[-1]
            i += 1
        elif c == ";":
            break
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in "(),;":
                j += 1
            i = j
        else:
            j = i
            while j < n and s[j] not in "(),:;":
                j += 1
            atual.nome += s[i:j].strip().strip("'")
            i = j
    return raiz.filhos[0] if len(raiz.filhos) == 1 and not raiz.nome else raiz


def main() -> None:
    raiz = ler((D / "labelled_supertree_ottnames.tre").read_text().strip())
    ann = json.loads((D / "annotations.json").read_text())["nodes"]
    rank = json.loads((D / "rank_collection.json").read_text())["decisions"]
    prioridade = {f"{d['studyID']}@{d['treeID']}": i for i, d in enumerate(rank)}
    ref = {f"{d['studyID']}@{d['treeID']}": d.get("compactRefText", "") for d in rank}
    # Especie -> no. Na sintese a especie pode ser no interno (com subespecies
    # como folhas: "Pandion_haliaetus_ott509837" acima de
    # "Pandion_haliaetus_haliaetus_ott806938") ou existir so pelas subespecies
    # ("Strix_occidentalis_caurina_ott486530"). Prefere o no da especie; senao,
    # a primeira subespecie (qualquer uma serve: todas estao dentro da especie).
    folhas: dict[str, No] = {}
    exato: set[str] = set()
    pilha = [raiz]
    while pilha:
        no = pilha.pop()
        pilha.extend(no.filhos)
        m = re.match(r"(.+)_ott\d+$", no.nome)
        if not m:
            continue
        partes = m.group(1).split("_")
        if len(partes) < 2 or not partes[0][:1].isupper() or not partes[1][:1].islower():
            continue
        esp = f"{partes[0]} {partes[1]}"
        if len(partes) == 2:
            folhas[esp] = no
            exato.add(esp)
        elif esp not in exato and esp not in folhas:
            folhas[esp] = no

    def chave(rot: str) -> str:
        if rot.startswith("mrca"):
            return rot
        m = re.search(r"(ott\d+)$", rot)
        return m.group(1) if m else rot

    def caminho(no: No) -> list[No]:
        out = []
        x: No | None = no
        while x is not None:
            out.append(x)
            x = x.pai
        return out

    def mrca(nomes: list[str]) -> No:
        cams = [caminho(folhas[t]) for t in nomes]
        comum = set(map(id, cams[0]))
        for c in cams[1:]:
            comum &= set(map(id, c))
        for x in cams[0]:
            if id(x) in comum:
                return x
        raise RuntimeError

    linhas = []
    for nome, dentro, fora in NOS:
        d = [t for t in dentro if t in folhas]
        ausentes = [t for t in dentro if t not in folhas]
        no = mrca(d)
        sob = set()
        pilha = [no]
        while pilha:  # todos os nos do clado: a especie pode ser no interno
            x = pilha.pop()
            pilha.extend(x.filhos)
            sob.add(id(x))
        intrusos = [t for t in fora if t in folhas and id(folhas[t]) in sob]
        a = ann.get(chave(no.nome), {})
        sup = a.get("supported_by", {}) or {}
        con = a.get("conflicts_with", {}) or {}
        par = a.get("partial_path_of", {}) or {}
        topo = min(sup, key=lambda k: prioridade.get(k, 10**6)) if sup else ""
        linhas.append([
            nome, no.nome, "sim" if not intrusos else "nao (inclui: " + "; ".join(intrusos[:4]) + ")",
            str(len(sup)), str(len(con)), str(len(par)),
            f"{topo} ({ref.get(topo, '')}; prioridade {prioridade.get(topo, '')})" if topo else "",
            ";".join(sorted(sup)), ";".join(sorted(con)), ";".join(ausentes),
        ])
        print(nome, no.nome, len(sup), len(con), intrusos[:2])
    cab = ["no", "node_id_aves16", "existe_na_sintese", "n_sustentam", "n_conflitam", "n_caminho_parcial",
           "fonte_de_maior_prioridade_que_sustenta", "sustentam", "conflitam", "exemplares_ausentes"]
    with open(D / "nos_aves16.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cab) + "\n")
        for li in linhas:
            f.write("\t".join(li) + "\n")


if __name__ == "__main__":
    main()
