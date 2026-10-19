"""
Idade (Ma) e intervalo dos nos dos pares do TimeTree em duas arvores datadas
publicadas, para comparar com a mediana do TimeTree:

  - B10K, Stiller et al. 2024: 04_timetree/main_alternative_mod.tre (MCMCtree,
    rotulo height_95%_HPD por no; numeros com tabela Translate);
  - Aves 1.6 (McTavish et al. 2025, versao de 2026 no AvesData):
    Clements2025/summary_dated_clements.nex (arvore MCC de 100 arvores datadas
    pelo chronosynth; rotulo age_95%_HPD).

    uv run python scripts/rapinantes/reanalise/idades_arvores_datadas.py <b10k.tre> <aves.nex>

A idade do no e a altura (distancia ate as folhas; as duas arvores sao
ultrametricas). Grava runs/rapinantes/reanalise/timetree/idades_arvores_datadas.tsv.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]
OUT = RAIZ / "runs" / "rapinantes" / "reanalise" / "timetree"

PARES = [
    ("Accipitridae-Cathartidae", "Aquila_chrysaetos", "Cathartes_aura"),
    ("Accipitridae-Sagittariidae", "Aquila_chrysaetos", "Sagittarius_serpentarius"),
    ("Accipitridae-Pandionidae", "Aquila_chrysaetos", "Pandion_haliaetus"),
    ("Falconidae-Psittaciformes", "Falco_peregrinus", "Nestor_notabilis"),
    ("Falconidae-Cariamidae", "Falco_peregrinus", "Cariama_cristata"),
    ("Tytonidae-Strigidae", "Tyto_alba", "Strix_occidentalis"),
    ("Strigiformes-Accipitriformes", "Tyto_alba", "Aquila_chrysaetos"),
    ("Herpetotheres-Falco", "Herpetotheres_cachinnans", "Falco_peregrinus"),
    ("Caracara-Falco", "Caracara_plancus", "Falco_peregrinus"),
]


class No:
    __slots__ = ("filhos", "nome", "comp", "nota", "pai")

    def __init__(self) -> None:
        self.filhos: list[No] = []
        self.nome = ""
        self.comp = 0.0
        self.nota = ""
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
        elif c == "[":
            j = s.index("]", i)
            atual.nota += s[i + 1 : j]
            i = j + 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in "(),;[":
                j += 1
            atual.comp = float(s[i + 1 : j])
            i = j
        elif c == ";":
            break
        else:
            j = i
            while j < n and s[j] not in "(),:;[":
                j += 1
            atual.nome += s[i:j].strip()
            i = j
    return raiz.filhos[0] if len(raiz.filhos) == 1 and not raiz.nome else raiz


def ler_nexus(caminho: Path) -> No:
    t = caminho.read_text()
    trad = {}
    m = re.search(r"translate(.*?);", t, re.IGNORECASE | re.DOTALL)
    if m:
        for linha in m.group(1).split(","):
            p = linha.split()
            if len(p) == 2:
                trad[p[0]] = p[1]
    mt = re.search(r"^\s*tree\s+\S+\s*=\s*(?:\[&[RU]\]\s*)?(.*;)", t, re.IGNORECASE | re.MULTILINE)
    if mt is None:
        raise SystemExit(f"{caminho}: bloco TREE nao encontrado")
    raiz = ler(mt.group(1))
    if trad:
        pilha = [raiz]
        while pilha:
            no = pilha.pop()
            pilha.extend(no.filhos)
            if not no.filhos and no.nome in trad:
                no.nome = trad[no.nome]
    return raiz


def folhas(raiz: No) -> dict[str, No]:
    out = {}
    pilha = [raiz]
    while pilha:
        no = pilha.pop()
        pilha.extend(no.filhos)
        if not no.filhos:
            out[no.nome] = no
    return out


def altura(no: No) -> float:
    h = 0.0
    while no.filhos:
        no = no.filhos[0]
        h += no.comp
    return h


def mrca(a: No, b: No) -> No:
    anc = set()
    x: No | None = a
    while x is not None:
        anc.add(id(x))
        x = x.pai
    y: No | None = b
    while y is not None and id(y) not in anc:
        y = y.pai
    assert y is not None
    return y


def hpd(nota: str) -> str:
    m = re.search(r"(?:height|age)_95%_HPD=\{([^,}]+),([^}]+)\}", nota)
    return f"{float(m.group(1)):.2f}-{float(m.group(2)):.2f}" if m else ""


def main() -> None:
    b10k = ler_nexus(Path(sys.argv[1]))
    aves = ler_nexus(Path(sys.argv[2]))
    fb, fa = folhas(b10k), folhas(aves)
    linhas = ["par\ttaxon_a\ttaxon_b\tB10K_idade_Ma\tB10K_HPD95\tAves1.6_idade_Ma\tAves1.6_HPD95"]
    for par, a, b in PARES:
        cel = [par, a, b]
        for f in (fb, fa):
            if a in f and b in f:
                no = mrca(f[a], f[b])
                cel += [f"{altura(no):.2f}", hpd(no.nota)]
            else:
                cel += ["ausente", ""]
        linhas.append("\t".join(cel))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "idades_arvores_datadas.tsv").write_text("\n".join(linhas) + "\n")
    print("\n".join(linhas))


if __name__ == "__main__":
    main()
