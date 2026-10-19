"""
Le as saidas do ASTER (freqQuad.csv), do gCF (IQ-TREE) e as arvores de especies
publicadas do B10K e monta, por ramo contestado das aves de rapina, as
frequencias de quarteto das tres resolucoes, a probabilidade local, o teste de
politomia e o teste de simetria das duas alternativas.

    uv run python scripts/rapinantes/reanalise/sintese_b10k.py

Cada pergunta e uma quadriparticao (X, Y, Z, R = resto). Num freqQuad.csv o
ramo procurado e o que tem exatamente esses quatro grupos em volta, e cada uma
das tres linhas t1-t3 diz qual par (X+Y, X+Z ou Y+Z) aquela topologia junta;
t1 e a resolucao da arvore pontuada. Isso funciona com a arvore que tiver
qualquer uma das tres resolucoes, desde que X, Y e Z sejam clados nela.

Testes (sobre as contagens efetivas f1-f3 e EN que o ASTER grava):
  - politomia (Sayyari & Mirarab 2018): qui-quadrado com 2 g.l. de f1, f2, f3
    contra EN/3 cada; p alto = nao se rejeita a politomia;
  - simetria das alternativas: qui-quadrado com 1 g.l. de f2 contra f3; sob ILS
    pura as duas alternativas tem a mesma frequencia esperada. p baixo = uma
    alternativa domina (introgressao, paralogia, atracao de ramos longos,
    erro sistematico de estimativa das arvores de genes).
Saidas em runs/rapinantes/reanalise/b10k/: quartetos_ramos.tsv,
arvores_publicadas_ramos.tsv, gcf_ramos.tsv.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]
B = RAIZ / "runs" / "rapinantes" / "reanalise" / "b10k"
ASTER = B / "aster"


# ---------------------------------------------------------------- arvores
class No:
    __slots__ = ("filhos", "nome", "rot")

    def __init__(self) -> None:
        self.filhos: list[No] = []
        self.nome = ""
        self.rot = ""


def ler(s: str) -> No:
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
        elif c == "'":
            j = s.index("'", i + 1)
            atual.rot += s[i + 1 : j]
            i = j + 1
        elif c == "[":
            j = s.index("]", i)
            i = j + 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in "(),;[":
                j += 1
            i = j
        elif c == ";":
            break
        else:
            j = i
            while j < n and s[j] not in "(),:;'[":
                j += 1
            tok = s[i:j].strip()
            if atual.filhos:
                atual.rot += tok
            else:
                atual.nome += tok
            i = j
    return raiz.filhos[0] if len(raiz.filhos) == 1 and not raiz.nome else raiz


def clados(raiz: No) -> list[tuple[frozenset[str], str]]:
    """(folhas abaixo, rotulo) para cada no interno."""
    out = []

    def rec(no: No) -> frozenset[str]:
        if not no.filhos:
            return frozenset([no.nome])
        s = frozenset().union(*(rec(f) for f in no.filhos))
        out.append((s, no.rot))
        return s

    rec(raiz)
    return out


def folhas(raiz: No) -> frozenset[str]:
    if not raiz.filhos:
        return frozenset([raiz.nome])
    return frozenset().union(*(folhas(f) for f in raiz.filhos))


def mrca_conj(cl: list[tuple[frozenset[str], str]], taxa: list[str]) -> frozenset[str]:
    alvo = set(taxa)
    melhor = None
    for s, _ in cl:
        if alvo <= s and (melhor is None or len(s) < len(melhor)):
            melhor = s
    assert melhor is not None
    return melhor


# ---------------------------------------------------------------- grupos
REF = ler((ASTER / "63K_topologia.tre").read_text().strip())
CL = clados(REF)
TODOS = folhas(REF)


def g(*taxa: str) -> frozenset[str]:
    return mrca_conj(CL, list(taxa)) if len(taxa) > 1 else frozenset(taxa)


GRUPOS = {
    "Accipitrimorphae": g("Cathartes_aura", "Aquila_chrysaetos"),
    "Strigiformes": g("Tyto_alba", "Strix_occidentalis"),
    "Coraciimorphae": g("Colius_striatus", "Picoides_pubescens"),
    "Cariamiformes": g("Cariama_cristata", "Chunga_burmeisteri"),
    "Falconidae": g("Falco_peregrinus", "Herpetotheres_cachinnans"),
    "Psittacopasserae": g("Nestor_notabilis", "Corvus_cornix"),
    "Psittaciformes": g("Nestor_notabilis", "Melopsittacus_undulatus"),
    "Passeriformes": g("Acanthisitta_chloris", "Corvus_cornix"),
    "Cathartidae": g("Cathartes_aura"),
    "Sagittariidae": g("Sagittarius_serpentarius"),
    "Pandionidae": g("Pandion_haliaetus"),
    "Accipitridae": g("Circaetus_pectoralis", "Aquila_chrysaetos"),
    "Tytonidae": g("Tyto_alba"),
    "Glaucidium": g("Glaucidium_brasilianum"),
    "Strix+Ciccaba": g("Strix_occidentalis", "Ciccaba_nigrolineata"),
}
GRUPOS["Pandionidae+Accipitridae"] = GRUPOS["Pandionidae"] | GRUPOS["Accipitridae"]
GRUPOS["Accipitrimorphae+Strigiformes"] = GRUPOS["Accipitrimorphae"] | GRUPOS["Strigiformes"]
GRUPOS["Strigiformes+Coraciimorphae"] = GRUPOS["Strigiformes"] | GRUPOS["Coraciimorphae"]
GRUPOS["Australaves"] = GRUPOS["Cariamiformes"] | GRUPOS["Falconidae"] | GRUPOS["Psittacopasserae"]

PERGUNTAS = {
    "afroaves_base": ("Accipitrimorphae", "Strigiformes", "Coraciimorphae"),
    "australaves_base": ("Cariamiformes", "Falconidae", "Psittacopasserae"),
    "psittacopasserae": ("Falconidae", "Psittaciformes", "Passeriformes"),
    "accipitrimorphae": ("Cathartidae", "Sagittariidae", "Pandionidae+Accipitridae"),
    "accipitriformes_ss": ("Sagittariidae", "Pandionidae", "Accipitridae"),
    "strigidae": ("Glaucidium", "Strix+Ciccaba", "Tytonidae"),
    # Monofilia de Afroaves. Em quartetos, "Accipitrimorphae irma do resto de
    # Afroaves" e "Accipitrimorphae irma de todo o resto de Telluraves" (Prum et al.
    # 2015; wASTRAL de Sapci et al. 2026) caem na MESMA resolucao Str+Cor da
    # pergunta afroaves_base; o que as separa e este ramo, que so existe na arvore
    # que tiver o clado do meio. Por isso duas versoes:
    "afroaves_mono_se_AS": ("Accipitrimorphae+Strigiformes", "Coraciimorphae", "Australaves"),
    "afroaves_mono_se_SC": ("Accipitrimorphae", "Strigiformes+Coraciimorphae", "Australaves"),
}


def alvo(p: str, taxa: frozenset[str]) -> tuple[dict[str, frozenset[str]], frozenset[frozenset[str]]]:
    x, y, z = PERGUNTAS[p]
    gs = {k: GRUPOS[k] & taxa for k in (x, y, z)}
    gs["R"] = taxa - gs[x] - gs[y] - gs[z]
    return gs, frozenset(gs.values())


# ---------------------------------------------------------------- testes
def p_qui2(x: float, gl: int) -> float:
    if gl == 2:
        return math.exp(-x / 2)
    return math.erfc(math.sqrt(x / 2))


# ---------------------------------------------------------------- freqQuad
QUAD = re.compile(r"\{([^}]*)\}\|\{([^}]*)\}#\{([^}]*)\}\|\{([^}]*)\}")


def ler_freqquad(caminho: Path) -> dict[str, list[tuple[str, list[frozenset[str]], float, float, float]]]:
    ramos: dict[str, list] = {}
    for linha in caminho.read_text().splitlines():
        p = linha.split("\t")
        if len(p) < 6:
            continue
        m = QUAD.fullmatch(p[2])
        if not m:
            continue
        sets = [frozenset(x.split(",")) for x in m.groups()]
        ramos.setdefault(p[0], []).append((p[1], sets, float(p[3]), float(p[4]), float(p[5])))
    return ramos


def nome_par(gs: dict[str, frozenset[str]], a: frozenset[str], b: frozenset[str],
             c: frozenset[str], d: frozenset[str]) -> str:
    inv = {v: k for k, v in gs.items()}
    par1 = {inv[a], inv[b]}
    par2 = {inv[c], inv[d]}
    junto = par2 if "R" in par1 else par1
    x, y = sorted(junto, key=list(gs).index)
    return f"{x}+{y}"


def quartetos(analise: str, caminho: Path, taxa: frozenset[str]) -> list[dict]:
    ramos = ler_freqquad(caminho)
    linhas = []
    for p in PERGUNTAS:
        gs, quad = alvo(p, taxa)
        x, y, z = PERGUNTAS[p]
        alt = [f"{x}+{y}", f"{x}+{z}", f"{y}+{z}"]
        achou = None
        for nid, ts in ramos.items():
            if frozenset(ts[0][1]) == quad:
                achou = (nid, ts)
                break
        li: dict = {"analise": analise, "pergunta": p}
        if achou is None:
            li["resolucao_na_arvore"] = "ramo ausente (X, Y ou Z nao e clado na arvore pontuada)"
            linhas.append(li)
            continue
        nid, ts = achou
        dados = {}
        en = ts[0][4]
        for t, sets, pp, f, _en in ts:
            dados[nome_par(gs, *sets)] = (t, pp, f)
        principal = next(k for k, v in dados.items() if v[0] == "t1")
        f = [dados[a][2] for a in alt]
        e = en / 3
        chi_pol = sum((fi - e) ** 2 / e for fi in f) if en > 0 else float("nan")
        # ordem pela frequencia, nao pela arvore: com localPP ~0 a resolucao da
        # arvore (t1) nao e a mais frequente
        ordem = sorted(alt, key=lambda a: -dados[a][2])
        f1, f2, f3 = (dados[a][2] for a in ordem)
        chi_top = (f1 - f2) ** 2 / (f1 + f2) if (f1 + f2) > 0 else 0.0
        chi_as = (f2 - f3) ** 2 / (f2 + f3) if (f2 + f3) > 0 else 0.0
        li.update({
            "resolucao_na_arvore": principal,
            "localPP_da_resolucao": round(dados[principal][1], 4),
            "EN": round(en, 1),
            **{f"q[{a}]": round(dados[a][2] / en, 4) for a in alt},
            **{f"pp[{a}]": round(dados[a][1], 4) for a in alt},
            "mais_frequente": ordem[0],
            "p_politomia": f"{p_qui2(chi_pol, 2):.3g}",
            "p_top2 (1a vs 2a)": f"{p_qui2(chi_top, 1):.3g}",
            "alternativas_menores": " vs ".join(ordem[1:]),
            "p_simetria_menores": f"{p_qui2(chi_as, 1):.3g}",
            "ramo": nid,
        })
        linhas.append(li)
    return linhas


# ---------------------------------------------------------------- arvores publicadas
def resolucao_publicada(arv: No, p: str) -> tuple[str, str]:
    taxa = folhas(arv)
    gs, _ = alvo(p, taxa)
    x, y, z = PERGUNTAS[p]
    cl = clados(arv)
    # biparticao: clado == par, ou clado == complemento do par (arvore enraizada em outro lugar)
    for a, b in ((x, y), (x, z), (y, z)):
        alvo_s = gs[a] | gs[b]
        comp = taxa - alvo_s
        for s, rot in cl:
            if s == alvo_s or s == comp:
                return f"{a}+{b}", rot
    return "nenhuma (X, Y ou Z nao monofiletico)", ""


def raiz_telluraves(arv: No) -> str:
    """Topologia enraizada dos 4 grupos de Telluraves (A, S, C, U) numa arvore.

    Testa as 15 topologias enraizadas: cada clado delas (como biparticao,
    clado ou complemento) tem de existir na arvore. Serve para arvore com
    raiz em qualquer lugar, desde que fora de Telluraves.
    """
    import itertools

    taxa = folhas(arv)
    g = {"A": GRUPOS["Accipitrimorphae"] & taxa, "S": GRUPOS["Strigiformes"] & taxa,
         "C": GRUPOS["Coraciimorphae"] & taxa, "U": GRUPOS["Australaves"] & taxa}
    splits = set()
    for s, _ in clados(arv):
        splits.add(s)
        splits.add(taxa - s)

    def tem(*ks: str) -> bool:
        return frozenset().union(*(g[k] for k in ks)) in splits

    falhos = [k for k in g if not (tem(k) or len(g[k]) == 1)]
    if falhos:
        # diagnostico: que unioes de grupos (e subgrupos de U) sao clados mesmo assim
        extra = {"K": GRUPOS["Cariamiformes"] & taxa, "F": GRUPOS["Falconidae"] & taxa,
                 "P": GRUPOS["Psittacopasserae"] & taxa}
        g2 = g | extra
        pres = []
        for r in (2, 3):
            for comb in itertools.combinations("ASCKFP", r):
                if frozenset().union(*(g2[k] for k in comb)) in splits:
                    pres.append("+".join(comb))
        return f"nao monofiletico: {','.join(falhos)}; clados: {' '.join(pres)}"
    achadas = []
    for x, y, z, w in itertools.permutations("ASCU"):
        if x < y and tem(x, y) and tem(x, y, z):
            achadas.append(f"((({x},{y}),{z}),{w})")
        if x < y and z < w and x < z and tem(x, y) and tem(z, w):
            achadas.append(f"(({x},{y}),({z},{w}))")
    achadas = sorted(set(achadas))
    return achadas[0] if len(achadas) == 1 else ("nao resolvida" if not achadas else "ambigua: " + ",".join(achadas))


PUBLICADAS = {
    "63K ASTRAL (principal)": "03_species_trees_63K.tre",
    "63K concatenacao RAxML-NG": "03_species_trees_63K_concatenation.tre",
    "teste de politomia 63K (p)": "polytomy_test_63K.tre",
}
PUBLICADAS_WSL = {
    "63K GOBS (bootstrap de locos)": "03_species_trees_63K-GOBS.tre",
    "exons C12 ASTRAL": "03_species_trees_exon_c12.tre",
    "introns ASTRAL": "03_species_trees_intron.tre",
    "UCE ASTRAL": "03_species_trees_UCE.tre",
    "94K ASTRAL": "03_species_trees_94K.tre",
    "128K ASTRAL": "03_species_trees_128K.tre",
    "153K ASTRAL": "03_species_trees_153K.tre",
    "159K ASTRAL": "03_species_trees_159K.tre",
}


def main() -> None:
    # 1. quartetos
    focal = frozenset(x.strip() for x in (B / "taxa_focais.txt").read_text().splitlines() if x.strip())
    linhas = []
    for tipo in ("intergenico63K", "exon_c12", "intron", "uce"):
        for modo, desc, taxa in (
            ("focal.score63K", "42 taxa, ASTRAL-IV, pontua topologia 63K", focal),
            ("focal.astral", "42 taxa, ASTRAL-IV, busca livre", focal),
            ("focal.wastral", "42 taxa, wASTRAL hibrido (aLRT), busca livre", focal),
            ("full363.score63K", "363 taxa, ASTRAL-IV, pontua topologia 63K", TODOS),
        ):
            arq = ASTER / f"{tipo}.{modo}.freqQuad.csv"
            if arq.exists():
                linhas += quartetos(f"{tipo} | {desc}", arq, taxa)
    cab = []
    for li in linhas:
        for k in li:
            if k not in cab:
                cab.append(k)
    with open(B / "quartetos_ramos.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cab) + "\n")
        for li in linhas:
            f.write("\t".join(str(li.get(c, "")) for c in cab) + "\n")
    print(f"quartetos_ramos.tsv: {len(linhas)} linhas")

    # 2. arvores publicadas
    out = ["arvore\tpergunta\tresolucao\trotulo_do_ramo(suporte ou p)"]
    pub = B / "publicados"
    todas = dict(PUBLICADAS) | dict(PUBLICADAS_WSL)
    for nome, arq in todas.items():
        cam = pub / arq
        if not cam.exists():
            continue
        arv = ler(cam.read_text().strip())
        for p in PERGUNTAS:
            r, rot = resolucao_publicada(arv, p)
            out.append(f"{nome}\t{p}\t{r}\t{rot}")
        out.append(f"{nome}\traiz_telluraves_A_S_C_U\t{raiz_telluraves(arv)}\t")
    # arvores de especies desta reanalise (busca livre, 42 taxa)
    for tipo in ("intergenico63K", "exon_c12", "intron", "uce"):
        for modo in ("astral", "wastral"):
            cam = ASTER / f"{tipo}.focal.{modo}.tre"
            if cam.exists() and cam.stat().st_size > 0:
                arv = ler(cam.read_text().strip())
                out.append(f"reanalise {tipo} focal {modo}\traiz_telluraves_A_S_C_U\t{raiz_telluraves(arv)}\t")
    (B / "arvores_publicadas_ramos.tsv").write_text("\n".join(out) + "\n", encoding="utf-8")
    print("arvores_publicadas_ramos.tsv ok")

    # 3. gCF sobre a topologia 63K podada
    gl = ["tipo\tpergunta\tramo_63K\tgCF\tgDF1\tgDF2\tgDFP\tgN"]
    for tipo in ("intergenico63K", "exon_c12", "intron", "uce"):
        br = B / "gcf" / f"gcf_{tipo}.cf.branch"
        st = B / "gcf" / f"gcf_{tipo}.cf.stat"
        if not br.exists():
            continue
        arv = ler(br.read_text().strip())
        taxa = folhas(arv)
        stat = {}
        for linha in st.read_text().splitlines():
            if linha.startswith("#") or linha.startswith("ID"):
                continue
            c = linha.split("\t")
            stat[c[0]] = c
        for p in PERGUNTAS:
            r, rot = resolucao_publicada(arv, p)
            c = stat.get(rot)
            if c:
                gl.append(f"{tipo}\t{p}\t{r}\t{c[1]}\t{c[3]}\t{c[5]}\t{c[7]}\t{c[9]}")
            else:
                gl.append(f"{tipo}\t{p}\t{r}\t\t\t\t\t")
        _ = taxa
    (B / "gcf_ramos.tsv").write_text("\n".join(gl) + "\n", encoding="utf-8")
    print("gcf_ramos.tsv ok")


if __name__ == "__main__":
    main()
