"""
Contagem de votos por arvore-fonte da Open Tree para as resolucoes alternativas
de cada no contestado das aves de rapina (meta-analise de topologias).

    uv run python scripts/rapinantes/reanalise/opentree_votos.py

Le runs/rapinantes/reanalise/opentree/node_info_*.json (gerados por
opentree_nos.py), junta todas as arvores-fonte que sustentam, conflitam ou
resolvem algum no-chave, baixa cada uma (folhas com ott_id) e, por pergunta,
diz qual das tres resolucoes a arvore contem como clado enraizado:

  pergunta            grupos (X, Y, Z) + resto R
  afroaves_base       Accipitrimorphae, Strigiformes, Coraciimorphae
  australaves_base    Cariamiformes, Falconiformes, Psittacopasserae
  psittacopasserae    Falconiformes, Psittaciformes, Passeriformes
  accipitrimorphae    Cathartidae, Sagittariidae, Pandionidae+Accipitridae
  accipitriformes_ss  Sagittariidae, Pandionidae, Accipitridae
  falconidae          Herpetotheriinae, Polyborinae, Falconinae
  strigiformes        Tytonidae x Strigidae (monofilia reciproca)

Voto "X+Y" = existe clado com todas as folhas amostradas de X e Y e nenhuma de
Z nem de R. Sem nenhuma das tres: "nao resolvido" (politomia ou grupo nao
monofiletico). Falta um grupo: "nao informativo". Sem nenhuma folha de R a
resposta depende so da raiz posta pelo curador: "sem externo".
A pertenca de cada folha vem da taxonomia OTT (taxonomy/subtree). Folhas sem
ott_id (OTU nao mapeada) ficam fora. Somente biblioteca padrao.
"""

from __future__ import annotations

import glob
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

API = "https://api.opentreeoflife.org/v3"
RAIZ = Path(__file__).resolve().parents[3]
OUT = RAIZ / "runs" / "rapinantes" / "reanalise" / "opentree"
FONTES = OUT / "arvores_fonte"

GRUPOS_TAX = {
    # Primeiro o que fica FORA da contagem: folha dentro de um grupo de interesse
    # mas sem subgrupo atribuido (genero fossil de Falconidae, familia extra da
    # OTT em Strigiformes) nao pode virar "resto" (R) - seria um falconideo
    # contado como externo. Os grupos abaixo sobrescrevem este rotulo.
    "EXCL": ["Falconidae", "Strigiformes", "Accipitriformes"],
    "Cathartidae": ["Cathartidae"],
    "Sagittariidae": ["Sagittariidae"],
    "Pandionidae": ["Pandionidae"],
    "Accipitridae": ["Accipitridae"],
    "Tytonidae": ["Tytonidae"],
    # a OTT 3.7 tem um no Bubonidae separado de Strigidae
    "Strigidae": ["Strigidae", "Bubonidae"],
    "Coraciimorphae": ["Coliiformes", "Leptosomiformes", "Trogoniformes", "Bucerotiformes",
                       "Coraciiformes", "Piciformes"],
    "Cariamiformes": ["Cariamiformes"],
    "Herpetotheriinae": ["Herpetotheres", "Micrastur"],
    "Polyborinae": ["Caracara", "Milvago", "Phalcoboenus", "Daptrius", "Ibycter", "Spiziapteryx"],
    # a OTT 3.7 poe os quiriquiris (Falco sparverius etc.) em Cerchneis
    "Falconinae": ["Falco", "Cerchneis", "Polihierax", "Microhierax", "Neohierax"],
    "Psittaciformes": ["Psittaciformes"],
    "Passeriformes": ["Passeriformes"],
}
# grupos compostos
COMPOSTOS = {
    "Accipitrimorphae": ["Cathartidae", "Sagittariidae", "Pandionidae", "Accipitridae"],
    "Strigiformes": ["Tytonidae", "Strigidae"],
    "Falconiformes": ["Herpetotheriinae", "Polyborinae", "Falconinae"],
    "Psittacopasserae": ["Psittaciformes", "Passeriformes"],
    "Pandionidae+Accipitridae": ["Pandionidae", "Accipitridae"],
    "Strigiformes+Coraciimorphae": ["Tytonidae", "Strigidae", "Coraciimorphae"],
    "Australaves": ["Cariamiformes", "Herpetotheriinae", "Polyborinae", "Falconinae",
                    "Psittaciformes", "Passeriformes"],
}
PERGUNTAS = {
    "afroaves_base": ("Accipitrimorphae", "Strigiformes", "Coraciimorphae"),
    "australaves_base": ("Cariamiformes", "Falconiformes", "Psittacopasserae"),
    "psittacopasserae": ("Falconiformes", "Psittaciformes", "Passeriformes"),
    "accipitrimorphae": ("Cathartidae", "Sagittariidae", "Pandionidae+Accipitridae"),
    "accipitriformes_ss": ("Sagittariidae", "Pandionidae", "Accipitridae"),
    "falconidae": ("Herpetotheriinae", "Polyborinae", "Falconinae"),
    "strigiformes": ("Tytonidae", "Strigidae", None),
    # monofilia de Afroaves: X+Y = Accipitrimorphae dentro de Afroaves;
    # Y+Z = Accipitrimorphae irma do resto de Telluraves (Prum et al. 2015)
    "afroaves_mono": ("Accipitrimorphae", "Strigiformes+Coraciimorphae", "Australaves"),
}
TAXONOMIA = {"ot_2019", "ot_2770", "ott3.7draft3"}


def req(metodo: str, caminho: str, corpo: dict | None = None, texto: bool = False):
    url = f"{API}/{caminho}"
    for k in range(5):
        try:
            if metodo == "POST":
                r = urllib.request.Request(url, data=json.dumps(corpo).encode(),
                                           headers={"content-type": "application/json"})
            else:
                r = urllib.request.Request(url)
            with urllib.request.urlopen(r, timeout=300) as resp:
                b = resp.read().decode()
                return b if texto else json.loads(b)
        except Exception as e:  # noqa: BLE001 - tenta de novo; a ultima falha sobe
            if k == 4:
                raise
            print(f"  {caminho}: {e}; nova tentativa", file=sys.stderr)
            time.sleep(4 * (k + 1))
    raise RuntimeError


class No:
    __slots__ = ("filhos", "nome")

    def __init__(self) -> None:
        self.filhos: list[No] = []
        self.nome: str | None = None


def ler_newick(s: str) -> No:
    raiz = No()
    pilha = [raiz]
    atual = raiz
    i, n = 0, len(s)
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
            novo = No()
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
        elif c == "'":
            j = s.index("'", i + 1)
            if not atual.filhos:
                atual.nome = s[i + 1 : j]
            i = j + 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in "(),;":
                j += 1
            i = j
        else:
            j = i
            while j < n and s[j] not in "(),:;'":
                j += 1
            if not atual.filhos:
                atual.nome = s[i:j].strip()
            i = j
    return raiz.filhos[0] if len(raiz.filhos) == 1 else raiz


def ott_de_nomes(nomes: list[str]) -> dict[str, int]:
    r = req("POST", "tnrs/match_names", {"names": nomes, "context_name": "Birds"})
    out = {}
    for x in r["results"]:
        if x["matches"]:
            out[x["name"]] = x["matches"][0]["taxon"]["ott_id"]
    return out


def membros_taxonomicos() -> dict[int, str]:
    """ott_id -> grupo base (de GRUPOS_TAX), via taxonomy/subtree."""
    cache = OUT / "taxonomia_grupos.json"
    if cache.exists():
        return {int(k): v for k, v in json.loads(cache.read_text()).items()}
    nomes = sorted({t for v in GRUPOS_TAX.values() for t in v})
    ids = ott_de_nomes(nomes)
    falta = [t for t in nomes if t not in ids]
    if falta:
        raise SystemExit(f"taxa sem ott_id: {falta}")
    mapa: dict[int, str] = {}
    for g, taxa in GRUPOS_TAX.items():
        for t in taxa:
            sub = req("POST", "taxonomy/subtree", {"ott_id": ids[t], "label_format": "id"})
            for m in re.findall(r"ott(\d+)", sub["newick"]):
                mapa[int(m)] = g
            print(f"  taxonomia {t} (ott{ids[t]}) -> {g}")
    cache.write_text(json.dumps(mapa))
    (OUT / "taxonomia_grupos_ott.json").write_text(json.dumps(ids, indent=1))
    return mapa


def contagens(raiz: No, grupo_de_folha) -> list[dict[str, int]]:
    """Lista de contagens por grupo para cada no interno (pos-ordem, iterativa)."""
    res = []
    pilha: list[tuple[No, bool]] = [(raiz, False)]
    cont: dict[int, dict[str, int]] = {}
    while pilha:
        no, visto = pilha.pop()
        if not no.filhos:
            g = grupo_de_folha(no.nome)
            cont[id(no)] = {g: 1} if g else {}
            continue
        if not visto:
            pilha.append((no, True))
            for f in no.filhos:
                pilha.append((f, False))
            continue
        c: dict[str, int] = {}
        for f in no.filhos:
            for k, v in cont.pop(id(f)).items():
                c[k] = c.get(k, 0) + v
        cont[id(no)] = c
        res.append(c)
    return res


def expandir(g: str) -> list[str]:
    return COMPOSTOS.get(g, [g])


def avaliar(cs: list[dict[str, int]], tot: dict[str, int], pergunta) -> str:
    x, y, z = pergunta
    def n(g: str, c: dict[str, int]) -> int:
        return sum(c.get(b, 0) for b in expandir(g))
    tx, ty = n(x, tot), n(y, tot)
    if z is None:  # monofilia reciproca de x e y
        if tx == 0 or ty == 0:
            return "nao informativo"
        resto_tot = sum(tot.values()) - tx - ty
        mx = any(n(x, c) == tx and n(y, c) == 0 and sum(c.values()) == tx for c in cs) or tx == 1
        my = any(n(y, c) == ty and n(x, c) == 0 and sum(c.values()) == ty for c in cs) or ty == 1
        xy = any(n(x, c) == tx and n(y, c) == ty and sum(c.values()) == tx + ty for c in cs)
        if resto_tot == 0:
            return "sem externo"
        if xy and mx and my:
            return f"{x}|{y} (ambos monofileticos)"
        return f"nao monofileticos (x={mx}, y={my}, x+y={xy})"
    tz = n(z, tot)
    if min(tx, ty, tz) == 0:
        return "nao informativo"
    resto = sum(tot.values()) - tx - ty - tz
    alvo = {f"{x}+{y}": (x, y), f"{y}+{z}": (y, z), f"{x}+{z}": (x, z)}
    achou = []
    for rot, (a, b) in alvo.items():
        ta, tb = n(a, tot), n(b, tot)
        for c in cs:
            if n(a, c) == ta and n(b, c) == tb and sum(c.values()) == ta + tb:
                achou.append(rot)
                break
    if resto == 0:
        return "sem externo" + (f" ({achou[0]})" if achou else "")
    if len(achou) == 1:
        return achou[0]
    if not achou:
        return "nao resolvido"
    return "ambiguo: " + ",".join(achou)


def indexar(raiz: No, grupo_de_folha):
    """pai e profundidade por no (ids inteiros) e folhas por grupo base."""
    pai: dict[int, int] = {}
    prof: dict[int, int] = {id(raiz): 0}
    folhas: dict[str, list[int]] = {}
    pilha = [raiz]
    while pilha:
        no = pilha.pop()
        for f in no.filhos:
            pai[id(f)] = id(no)
            prof[id(f)] = prof[id(no)] + 1
            pilha.append(f)
        if not no.filhos:
            g = grupo_de_folha(no.nome)
            if g:
                folhas.setdefault(g, []).append(id(no))
    return pai, prof, folhas


def _dist(u: int, v: int, pai: dict[int, int], prof: dict[int, int]) -> int:
    a, b, d = u, v, 0
    while prof[a] > prof[b]:
        a = pai[a]
        d += 1
    while prof[b] > prof[a]:
        b = pai[b]
        d += 1
    while a != b:
        a, b = pai[a], pai[b]
        d += 2
    return d


def quartetos(pai, prof, folhas, pergunta, n: int = 2000, semente: int = 12345) -> str:
    """Fracao de quartetos (x, y, z, r) que sustenta cada resolucao.

    Robusto a um taxon fora do lugar, que derruba o criterio de clado estrito.
    Distancia topologica: ab|cd vale quando d(a,b)+d(c,d) e o menor dos tres
    pareamentos; empate = quarteto nao resolvido (politomia na arvore-fonte).
    """
    import random

    rng = random.Random(semente)
    x, y, z = pergunta

    def membros(g: str) -> list[int]:
        return [f for b in expandir(g) for f in folhas.get(b, [])]

    lx, ly = membros(x), membros(y)
    usados = set(expandir(x)) | set(expandir(y))
    if z is None:
        lr = [f for g, fs in folhas.items() if g not in usados for f in fs]
        if not lx or len(ly) < 2:
            return "nao informativo"
        if not lr:
            return "sem externo"
        c = {"mono": 0, "nao": 0, "nr": 0}
        for _ in range(n):
            t, r = rng.choice(lx), rng.choice(lr)
            s1, s2 = rng.sample(ly, 2)
            d = lambda a, b: _dist(a, b, pai, prof)  # noqa: E731
            v = (d(s1, s2) + d(t, r), d(t, s1) + d(s2, r), d(t, s2) + d(s1, r))
            m = min(v)
            if v.count(m) > 1:
                c["nr"] += 1
            elif v[0] == m:
                c["mono"] += 1
            else:
                c["nao"] += 1
        return f"{y}_mono={c['mono'] / n:.2f};nao={c['nao'] / n:.2f};nr={c['nr'] / n:.2f}"
    lz = membros(z)
    usados |= set(expandir(z))
    lr = [f for g, fs in folhas.items() if g not in usados for f in fs]
    if not (lx and ly and lz):
        return "nao informativo"
    if not lr:
        return "sem externo"
    c = {f"{x}+{y}": 0, f"{x}+{z}": 0, f"{y}+{z}": 0, "nr": 0}
    for _ in range(n):
        a, b, e, r = rng.choice(lx), rng.choice(ly), rng.choice(lz), rng.choice(lr)
        d = lambda u, v: _dist(u, v, pai, prof)  # noqa: E731
        v = (d(a, b) + d(e, r), d(a, e) + d(b, r), d(b, e) + d(a, r))
        m = min(v)
        if v.count(m) > 1:
            c["nr"] += 1
        elif v[0] == m:
            c[f"{x}+{y}"] += 1
        elif v[1] == m:
            c[f"{x}+{z}"] += 1
        else:
            c[f"{y}+{z}"] += 1
    return ";".join(f"{k}={v / n:.2f}" for k, v in c.items())


RAIZ_TELLURAVES = {"A": "Accipitrimorphae", "S": "Strigiformes", "C": "Coraciimorphae", "U": "Australaves"}
NOMES_RAIZ = {
    "(((A,S),C),U)": "Afroaves com Hieraves (Stiller 2024)",
    "(((C,S),A),U)": "Afroaves, Accipitrimorphae na base (Hackett 2008, Jarvis 2014)",
    "(((C,S),U),A)": "Accipitrimorphae irma de Telluraves (Prum 2015)",
    "((A,S),(C,U))": "Hieraves irma de Coraciimorphae+Australaves",
}


def raiz_telluraves(pai, prof, folhas, n: int = 2000, semente: int = 12345) -> str:
    """Topologia enraizada dos quatro grupos de Telluraves (A, S, C, U), por amostra.

    Cada amostra sorteia uma folha de cada grupo e uma externa (r). Cada um dos
    quatro tripletos sai do quarteto (x, y, z, r): xy|zr = x e y irmaos. Dois pares
    irmaos disjuntos = arvore balanceada; um par = pente, e o terceiro grupo sai
    do tripleto que sobra. Tripletos incoerentes ou empatados = "nr".
    """
    import random

    rng = random.Random(semente)
    lis = {k: [f for b in expandir(g) for f in folhas.get(b, [])] for k, g in RAIZ_TELLURAVES.items()}
    usados = {b for g in RAIZ_TELLURAVES.values() for b in expandir(g)}
    lr = [f for g, fs in folhas.items() if g not in usados for f in fs]
    if not all(lis.values()):
        return "nao informativo"
    if not lr:
        return "sem externo"
    d = lambda u, v: _dist(u, v, pai, prof)  # noqa: E731
    cont: dict[str, int] = {}
    for _ in range(n):
        f = {k: rng.choice(v) for k, v in lis.items()}
        r = rng.choice(lr)
        irmaos: dict[frozenset[str], int] = {}
        ok = True
        for tri in (("A", "S", "C"), ("A", "S", "U"), ("A", "C", "U"), ("S", "C", "U")):
            x, y, z = tri
            v = (d(f[x], f[y]) + d(f[z], r), d(f[x], f[z]) + d(f[y], r), d(f[y], f[z]) + d(f[x], r))
            m = min(v)
            if v.count(m) > 1:
                ok = False
                break
            par = (x, y) if v[0] == m else (x, z) if v[1] == m else (y, z)
            irmaos[frozenset(par)] = irmaos.get(frozenset(par), 0) + 1
        if not ok:
            cont["nr"] = cont.get("nr", 0) + 1
            continue
        # par irmao "cereja" = aparece como irmao nos dois tripletos que o contem
        cerejas = [p for p, k in irmaos.items() if k == 2]
        if len(cerejas) == 2 and not (cerejas[0] & cerejas[1]):
            a, b = sorted(",".join(sorted(p)) for p in cerejas)
            topo = f"(({a}),({b}))"
        elif len(cerejas) == 1:
            x, y = sorted(cerejas[0])
            resto = [k for k in "ASCU" if k not in cerejas[0]]
            # o terceiro: o que forma par com x ou y no tripleto com o ultimo
            p3 = [p for p in irmaos if p & cerejas[0] and not p <= cerejas[0]]
            terc = [k for k in resto if any(k in p for p in p3)]
            if len(terc) != 1:
                cont["nr"] = cont.get("nr", 0) + 1
                continue
            ult = [k for k in resto if k != terc[0]][0]
            topo = f"((({x},{y}),{terc[0]}),{ult})"
        else:
            cont["nr"] = cont.get("nr", 0) + 1
            continue
        cont[topo] = cont.get(topo, 0) + 1
    return ";".join(f"{k}={v / n:.2f}" for k, v in sorted(cont.items(), key=lambda kv: -kv[1]))


def main() -> None:
    FONTES.mkdir(parents=True, exist_ok=True)
    arvores = set()
    for f in glob.glob(str(OUT / "node_info_*.json")):
        ni = json.loads(Path(f).read_text())["node_info"]
        for k in ("supported_by", "conflicts_with", "partial_path_of", "terminal", "resolves"):
            arvores.update((ni.get(k) or {}).keys())
    arvores = sorted(a for a in arvores if "@" in a and a.split("@")[0] not in TAXONOMIA)
    estudos = json.loads((OUT / "estudos.json").read_text())
    tax = membros_taxonomicos()

    def grupo(nome: str | None) -> str | None:
        if not nome:
            return None
        m = re.fullmatch(r"(?:ott)?(\d+)", nome.strip())
        if not m:
            return None
        g = tax.get(int(m.group(1)), "R")
        return None if g == "EXCL" else g

    linhas = []
    for a in arvores:
        sid, tid = a.split("@")
        arq = FONTES / f"{sid}@{tid}.tre"
        if not arq.exists():
            txt = req("GET", f"study/{sid}/tree/{tid}.tre?tip_label=ot:ottId", texto=True)
            arq.write_text(txt)
        raiz = ler_newick(arq.read_text().strip())
        cs = contagens(raiz, grupo)
        tot: dict[str, int] = {}
        nfolhas = 0
        pilha = [raiz]
        while pilha:
            no = pilha.pop()
            if no.filhos:
                pilha.extend(no.filhos)
            else:
                nfolhas += 1
                g = grupo(no.nome)
                if g:
                    tot[g] = tot.get(g, 0) + 1
        m = estudos.get(sid, {})
        ref = (m.get("ot:studyPublicationReference") or "").split("(")[0].split(",")[0].strip()
        li = {"arvore": a, "ano": m.get("ot:studyYear", ""), "primeiro_autor": ref,
              "doi": m.get("ot:studyPublication", ""), "folhas": nfolhas,
              "folhas_mapeadas": sum(tot.values())}
        pai, prof, folhas = indexar(raiz, grupo)
        for p, q in PERGUNTAS.items():
            li[p] = avaliar(cs, tot, q)
            li[p + "_quartetos"] = quartetos(pai, prof, folhas, q)
        li["raiz_telluraves_A_S_C_U"] = raiz_telluraves(pai, prof, folhas)
        linhas.append(li)
        print(a, li["ano"], li["primeiro_autor"], {p: li[p] for p in PERGUNTAS})
    cab = list(linhas[0].keys())
    with open(OUT / "votos_por_arvore.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cab) + "\n")
        for li in linhas:
            f.write("\t".join(str(li[c]) for c in cab) + "\n")
    # resumo
    resumo = {}
    for p in PERGUNTAS:
        cnt: dict[str, int] = {}
        for li in linhas:
            cnt[li[p]] = cnt.get(li[p], 0) + 1
        resumo[p] = dict(sorted(cnt.items(), key=lambda kv: -kv[1]))
    (OUT / "votos_resumo.json").write_text(json.dumps(resumo, indent=1, ensure_ascii=False))
    print(json.dumps(resumo, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
