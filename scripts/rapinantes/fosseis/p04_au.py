"""Passo 5 (2c): hipoteses alternativas de posicao para fosseis-chave e teste AU.

  uv run python p04_au.py gerar     # au/<ID>/<fossil>/H_<rotulo>.nwk + lista.tsv
  (WSL) bash -l f04_au.sh           # busca ML restrita por hipotese + -z ... -au
  uv run python p04_au.py ler       # au/au_resultados.tsv

Cada hipotese "fossil irmao de X" vira uma arvore de restricao: o esqueleto
molecular (so atuais), enraizado no grupo externo, com o fossil enxertado no
ramo-tronco de X. Os outros fosseis ficam fora da restricao (livres), entao a
hipotese le-se "depois de podar os demais fosseis, o fossil e irmao de X". A
melhor arvore de cada hipotese, mais a arvore ML com fosseis livres, entram no
mesmo teste AU (10 000 reamostragens RELL).
"""

from __future__ import annotations

import csv
import re
import sys

import clados as cl
import dendropy
from p03_resumo_arvores import raiz_atual

AU = cl.BASE / "au"
COLS = ["logL", "deltaL", "bp-RELL", "p-KH", "p-SH", "p-WKH", "p-WSH", "c-ELW", "p-AU"]

# (matriz, fossil): lista de clados-alvo (nomes de clados.clados ou terminais;
# "A+B" = uniao, que precisa ser clado no esqueleto).
HIPOTESES: dict[tuple[str, str], list[str]] = {
    ("K19", "Pseudasturides_macrocephalus"): [
        "Psittaciformes",
        "Passeriformes",
        "Psittacopasserae",
        "Falconiformes",
        "Eufalconimorphae",
        "Cariamiformes",
        "Australaves",
        "Accipitriformes",
        "Strigiformes",
        "Hieraves",
        "Coraciimorphae",
    ],
    ("K19", "Cyrilavis_colburnorum"): [
        "Psittaciformes",
        "Passeriformes",
        "Psittacopasserae",
        "Falconiformes",
        "Eufalconimorphae",
        "Cariamiformes",
        "Australaves",
        "Accipitriformes",
        "Strigiformes",
        "Hieraves",
        "Coraciimorphae",
    ],
    ("M11", "Messelasturidae"): [
        "Psittaciformes",
        "Passeriformes",
        "Psittacopasserae",
        "Falconidae",
        "Eufalconimorphae",
        "Accipitridae",
        "Accipitriformes",
        "Strigiformes",
        "Hieraves",
        "Coraciimorphae",
    ],
    ("M11", "Halcyornithidae"): [
        "Psittaciformes",
        "Passeriformes",
        "Psittacopasserae",
        "Falconidae",
        "Eufalconimorphae",
        "Accipitridae",
        "Accipitriformes",
        "Strigiformes",
        "Hieraves",
        "Coraciimorphae",
    ],
    ("M05", "Messelastur"): [
        "Strigiformes",
        "Accipitridae",
        "Accipitriformes",
        "Hieraves",
        "Falconidae",
        "Psittaciformes",
        "Falconidae+Psittacidae",
        "Cariamidae",
        "Australaves_M05",
    ],
    ("M05", "Pseudasturidae"): [
        "Strigiformes",
        "Accipitridae",
        "Accipitriformes",
        "Hieraves",
        "Falconidae",
        "Psittaciformes",
        "Falconidae+Psittacidae",
        "Cariamidae",
        "Australaves_M05",
    ],
    ("VZ26", "Messelastur"): [
        "Psittaciformes",
        "Passeriformes",
        "Psittacopasserae",
        "Falconidae",
        "Eufalconimorphae",
        "Strigiformes",
        "Coraciimorphae",
    ],
    ("VZ26", "Pseudasturides"): [
        "Psittaciformes",
        "Passeriformes",
        "Psittacopasserae",
        "Falconidae",
        "Eufalconimorphae",
        "Strigiformes",
        "Coraciimorphae",
    ],
    ("MK22", "Masillaraptor"): ["Falconidae", "Cariama", "Falconidae+Cariama"],
    ("MK22", "Danielsraptor"): ["Falconidae", "Cariama", "Falconidae+Cariama"],
    ("M16", "Bathornis"): [
        "Cariama",
        "Falconidae",
        "Falconidae+Cariama",
        "Cathartidae",
        "Opisthocomidae",
        "Cathartidae+Falconidae+Cariama",
    ],
    ("D21", "Dryornis_pampeanus"): [
        "Cathartidae",
        "Ciconiidae",
        "Accipitriformes",
        "Accipitridae",
        "Sagittariidae",
        "Falconidae",
        "Pandionidae",
    ],
    ("D21", "Pampagyps_imperator"): [
        "Cathartidae",
        "Ciconiidae",
        "Accipitriformes",
        "Accipitridae",
        "Sagittariidae",
        "Falconidae",
        "Pandionidae",
    ],
    ("M16", "Elaphrocnemus"): [
        "Cariama",
        "Falconidae",
        "Falconidae+Cariama",
        "Cathartidae",
        "Opisthocomidae",
        "Cathartidae+Falconidae+Cariama",
    ],
}


def alvo(id_: str, nome: str, cls: dict[str, frozenset[str]]) -> frozenset[str]:
    if nome == "Australaves_M05":
        return cls["Cariamidae"] | cls["Falconidae"] | cls["Psittaciformes"]
    if "+" in nome:
        out: frozenset[str] = frozenset()
        for p in nome.split("+"):
            out |= cls[p]
        return out
    return cls[nome]


def esqueleto_enraizado(id_: str) -> dendropy.Tree:
    t = dendropy.Tree.get(
        path=str(cl.BASE / "matrizes" / id_ / "esqueleto_stiller.nwk"),
        schema="newick",
        preserve_underscores=True,
    )
    og = t.find_node_with_taxon_label(raiz_atual(id_))
    if og is None:
        raise ValueError(f"{id_}: externo ausente do esqueleto")
    t.reroot_at_edge(og.edge, update_bipartitions=False)
    return t


def enxertar(t: dendropy.Tree, fossil: str, x: frozenset[str]) -> str:
    a = t.clone(depth=2)
    folhas = [lf for lf in a.leaf_node_iter() if lf.taxon.label in x]
    no = a.mrca(taxa=[lf.taxon for lf in folhas]) if len(folhas) > 1 else folhas[0]
    if no is None:
        raise ValueError("mrca vazio")
    got = {lf.taxon.label for lf in no.leaf_iter()}
    if got != set(x):
        raise ValueError(
            f"alvo nao e clado no esqueleto: falta {set(x) - got}, sobra {got - set(x)}"
        )
    pai = no.parent_node
    if pai is None:
        raise ValueError("alvo e a raiz")
    novo = dendropy.Node()
    pai.remove_child(no)
    pai.add_child(novo)
    novo.add_child(no)
    novo.add_child(dendropy.Node(taxon=dendropy.Taxon(label=fossil)))
    return a.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True).strip()


def gerar() -> None:
    linhas = []
    for (id_, fossil), alvos in HIPOTESES.items():
        cls = cl.clados(id_)
        t = esqueleto_enraizado(id_)
        pasta = AU / id_ / fossil
        pasta.mkdir(parents=True, exist_ok=True)
        for nome in alvos:
            x = alvo(id_, nome, cls)
            rot = re.sub(r"[^A-Za-z0-9]+", "_", nome)
            (pasta / f"H_{rot}.nwk").write_text(enxertar(t, fossil, x) + "\n", encoding="utf-8")
            linhas.append([id_, fossil, rot, nome, len(x)])
    with open(AU / "lista.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerows(linhas)
    print(len(linhas), "hipoteses")


# Teste da topologia dos ATUAIS: a morfologia rejeita o esqueleto molecular?
# Hipoteses = esqueleto com Falconidae movido para junto dos gavioes (a uniao
# "Falconiformes" classica), mais a arvore morfologica livre.
ATUAIS_TESTES = {
    "K19": ["Accipitridae", "Accipitriformes", "Hieraves"],
    "M11": ["Accipitridae", "Accipitriformes", "Hieraves"],
    "M05": ["Accipitridae", "Accipitriformes", "Hieraves"],
    "MC03": ["Accipitridae", "Accipitriformes", "Hieraves"],
    "C25": ["Accipitridae", "Accipitriformes", "Hieraves"],
    "LZ07": ["Accipitridae", "Accipitriformes", "Hieraves", "Cathartidae>Ciconiidae"],
    "D21": ["Accipitridae", "Cathartidae>Ciconiidae"],
}
PASTA_ATUAIS = "atuais_Falconidae"


def mover(t: dendropy.Tree, grupo: frozenset[str], alvo_: frozenset[str]) -> str:
    a = t.clone(depth=2)

    def mrca(x: frozenset[str]) -> dendropy.Node:
        fs = [lf for lf in a.leaf_node_iter() if lf.taxon.label in x]
        no = a.mrca(taxa=[lf.taxon for lf in fs]) if len(fs) > 1 else fs[0]
        if no is None or {lf.taxon.label for lf in no.leaf_iter()} != set(x):
            raise ValueError(f"nao e clado: {sorted(x)[:4]}")
        return no

    g = mrca(grupo)
    pai = g.parent_node
    if pai is None:
        raise ValueError("grupo e a raiz")
    pai.remove_child(g)
    a.suppress_unifurcations()
    alvo_no = mrca(alvo_)
    pai2 = alvo_no.parent_node
    if pai2 is None:
        raise ValueError("alvo e a raiz")
    novo = dendropy.Node()
    pai2.remove_child(alvo_no)
    pai2.add_child(novo)
    novo.add_child(alvo_no)
    novo.add_child(g)
    return a.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True).strip()


def gerar_atuais() -> None:
    for id_, alvos in ATUAIS_TESTES.items():
        cls = cl.clados(id_)
        fal = cls.get("Falconidae")
        if not fal:
            continue
        t = esqueleto_enraizado(id_)
        pasta = AU / id_ / PASTA_ATUAIS
        pasta.mkdir(parents=True, exist_ok=True)
        for nome in alvos:
            # "A>B" move o grupo A para junto de B; sem ">", A = Falconidae
            grupo, alvo_nome = nome.split(">") if ">" in nome else ("Falconidae", nome)
            if alvo_nome not in cls or grupo not in cls:
                continue
            nwk = mover(t, cls[grupo], cls[alvo_nome])
            (pasta / f"H_{grupo}_com_{alvo_nome}.nwk").write_text(nwk + "\n", encoding="utf-8")
        livre = cl.BASE / "arvores" / id_ / "livre.treefile"
        if livre.exists():
            (pasta / "H_morfologia_livre.nwk").write_text(
                livre.read_text(encoding="utf-8"), encoding="utf-8"
            )
        print(id_, "ok")


def _ler_iqtree(rel, ordem: list[str]) -> list[list[str]]:
    txt = rel.read_text(encoding="utf-8", errors="replace")
    bloco = txt.split("USER TREES")[1]
    cab: list[str] = []
    out = []
    for lin in bloco.splitlines():
        if lin.strip().startswith("Tree") and "logL" in lin:
            cab = lin.split()[1:]
            continue
        m = re.match(r"\s*(\d+)\s+(-[\d.]+)\s+(.*)", lin)
        if not m or not cab:
            continue
        i = int(m.group(1))
        vals = [m.group(2)] + [c for c in m.group(3).split() if c not in ("+", "-")]
        d = dict(zip(cab, vals, strict=False))
        rot = ordem[i - 1] if i - 1 < len(ordem) else str(i)
        out.append([rot] + [d.get(k, "") for k in COLS])
    return out


def ler() -> None:
    saida = []
    chaves = list(HIPOTESES) + [(i, PASTA_ATUAIS) for i in ATUAIS_TESTES]
    for id_, fossil in chaves:
        pasta = AU / id_ / fossil
        for rel, ord_, rotulo in (
            ("AU.iqtree", "ordem.txt", fossil),
            ("AU_esqueletos.iqtree", "ordem_esqueletos.txt", fossil + "_so_esqueletos"),
        ):
            if not (pasta / rel).exists() or not (pasta / ord_).exists():
                continue
            ordem = (pasta / ord_).read_text(encoding="utf-8").split()
            for lin in _ler_iqtree(pasta / rel, ordem):
                saida.append([id_, rotulo, *lin])
    with open(AU / "au_resultados.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["matriz", "fossil", "hipotese", *COLS])
        w.writerows(saida)
    print(len(saida), "linhas")


def lz07_spr() -> None:
    """LZ07: alternativas por SPR sobre a arvore ML com esqueleto (sem nova busca).
    A arvore alternativa nao e reotimizada, o que joga contra ela no AU."""
    cls = cl.clados("LZ07")
    t = dendropy.Tree.get(
        path=str(cl.BASE / "arvores" / "LZ07" / "esqueleto.treefile"),
        schema="newick",
        preserve_underscores=True,
    )
    og = t.find_node_with_taxon_label(raiz_atual("LZ07"))
    if og is None:
        raise ValueError("raiz ausente")
    t.reroot_at_edge(og.edge, update_bipartitions=False)
    pasta = AU / "LZ07" / PASTA_ATUAIS
    pasta.mkdir(parents=True, exist_ok=True)
    base = t.as_string(
        schema="newick",
        suppress_rooting=True,
        unquoted_underscores=True,
        suppress_edge_lengths=True,
        suppress_internal_node_labels=True,
    ).strip()
    linhas, ordem = [base], ["ML_fosseis_livres"]
    for g, a in (("Falconidae", "Accipitridae"), ("Cathartidae", "Ciconiidae")):
        nwk = mover(t, cls[g], cls[a])
        # o R (p06) le a alternativa falcao-gaviao com este nome
        (pasta / f"H_{g}_com_{a}.treefile").write_text(nwk + "\n", encoding="utf-8")
        linhas.append(nwk)
        ordem.append(f"{g}_com_{a}_SPR")
    (pasta / "esqueletos.trees").write_text("\n".join(linhas) + "\n", encoding="utf-8")
    (pasta / "ordem_esqueletos.txt").write_text("\n".join(ordem) + "\n", encoding="utf-8")
    print("LZ07 SPR ok")


if __name__ == "__main__":
    {"gerar": gerar, "gerar_atuais": gerar_atuais, "ler": ler, "lz07_spr": lz07_spr}[sys.argv[1]]()
