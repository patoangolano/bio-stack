"""
Poda arvores de genes (uma por linha, Newick) para um conjunto de taxa.

    python3 -I podar_arvores.py <entrada.tre> <taxa.txt> <saida.tre> [--min 4]

Por que nao `gotree prune -r`: o gotree 0.5.2 desta camada aborta o arquivo
inteiro quando uma arvore fica com 2 folhas depois da poda (medido no UCE do
B10K). Aqui cada arvore e podada sozinha; a que fica com menos de `--min`
taxa focais sai do arquivo (nao tem quarteto, o ASTRAL nao a usaria) e entra
na contagem impressa no fim.

No no que fica com um filho so depois da poda, os dois ramos viram um:
comprimento somado e suporte = o MAIOR dos dois, porque cada um dos ramos
originais, restrito aos taxa que ficaram, define a mesma biparticao.
Somente biblioteca padrao.
"""

from __future__ import annotations

import argparse
import sys


class No:
    __slots__ = ("filhos", "nome", "comp", "suporte")

    def __init__(self) -> None:
        self.filhos: list[No] = []
        self.nome: str | None = None
        self.comp: float | None = None
        self.suporte: str | None = None


def ler(s: str) -> No:
    i = 0
    n = len(s)
    raiz = No()
    pilha = [raiz]
    atual = raiz
    while i < n:
        c = s[i]
        if c == "(":
            novo = No()
            atual.filhos.append(novo)
            pilha.append(novo)
            atual = novo
            i += 1
        elif c == ",":
            pilha.pop()
            pai = pilha[-1]
            novo = No()
            pai.filhos.append(novo)
            pilha.append(novo)
            atual = novo
            i += 1
        elif c == ")":
            pilha.pop()
            atual = pilha[-1]
            i += 1
        elif c == ";":
            break
        else:
            j = i
            while j < n and s[j] not in "(),:;":
                j += 1
            rot = s[i:j].strip()
            comp = None
            if j < n and s[j] == ":":
                k = j + 1
                while k < n and s[k] not in "(),;":
                    k += 1
                try:
                    comp = float(s[j + 1 : k])
                except ValueError:
                    comp = None
                j = k
            if atual.filhos:
                atual.suporte = rot or None
            else:
                atual.nome = rot
            atual.comp = comp
            i = j
    # a raiz implicita tem um filho so: a arvore de verdade
    return raiz.filhos[0] if len(raiz.filhos) == 1 and raiz.nome is None else raiz


def _max_sup(a: str | None, b: str | None) -> str | None:
    try:
        fa = float(a) if a is not None else None
    except ValueError:
        fa = None
    try:
        fb = float(b) if b is not None else None
    except ValueError:
        fb = None
    if fa is None:
        return b
    if fb is None:
        return a
    return a if fa >= fb else b


def podar(no: No, manter: set[str]) -> No | None:
    if not no.filhos:
        return no if no.nome in manter else None
    filhos = [f for f in (podar(x, manter) for x in no.filhos) if f is not None]
    if not filhos:
        return None
    if len(filhos) == 1:
        f = filhos[0]
        if no.comp is not None or f.comp is not None:
            f.comp = (no.comp or 0.0) + (f.comp or 0.0)
        if f.filhos:
            f.suporte = _max_sup(f.suporte, no.suporte)
        return f
    no.filhos = filhos
    return no


def contar(no: No) -> int:
    return 1 if not no.filhos else sum(contar(f) for f in no.filhos)


def escrever(no: No, raiz: bool = True) -> str:
    if no.filhos:
        s = "(" + ",".join(escrever(f, False) for f in no.filhos) + ")"
        if not raiz and no.suporte is not None:
            s += no.suporte
    else:
        s = no.nome or ""
    if not raiz and no.comp is not None:
        s += f":{no.comp:.10g}"
    return s + (";" if raiz else "")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("entrada")
    ap.add_argument("taxa")
    ap.add_argument("saida")
    ap.add_argument("--min", type=int, default=4)
    a = ap.parse_args()
    manter = {x.strip() for x in open(a.taxa) if x.strip()}
    sys.setrecursionlimit(10000)
    n_in = n_out = n_fora = 0
    hist: dict[int, int] = {}
    with open(a.entrada) as fin, open(a.saida, "w") as fout:
        for linha in fin:
            linha = linha.strip()
            if not linha:
                continue
            n_in += 1
            t = podar(ler(linha), manter)
            k = contar(t) if t is not None else 0
            hist[k] = hist.get(k, 0) + 1
            if t is None or k < a.min:
                n_fora += 1
                continue
            # raiz com 2 filhos: desenraiza juntando um filho interno a raiz,
            # para nao ficar um ramo "fantasma" com suporte.
            if len(t.filhos) == 2:
                f0, f1 = t.filhos
                interno = f0 if f0.filhos else f1
                outro = f1 if interno is f0 else f0
                if interno.filhos:
                    if outro.comp is not None or interno.comp is not None:
                        outro.comp = (outro.comp or 0.0) + (interno.comp or 0.0)
                    if outro.filhos:
                        outro.suporte = _max_sup(outro.suporte, interno.suporte)
                    t.filhos = interno.filhos + [outro]
            fout.write(escrever(t) + "\n")
            n_out += 1
    media = sum(k * v for k, v in hist.items()) / max(n_in, 1)
    print(f"{a.entrada}: {n_in} arvores lidas, {n_out} escritas, {n_fora} com < {a.min} taxa focais; "
          f"media de taxa focais por arvore {media:.1f}")


if __name__ == "__main__":
    main()
