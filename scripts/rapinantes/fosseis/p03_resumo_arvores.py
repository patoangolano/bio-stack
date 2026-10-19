"""Passo 4: leitura das arvores do IQ-TREE.

  2a grupos tradicionais x moleculares: frequencia nas 1000 arvores UFBoot
     (podadas aos atuais) da analise morfologica livre, e presenca na arvore ML.
  2b posicao de cada fossil na analise com esqueleto molecular: para cada
     arvore (ML e UFBoot), poda os outros fosseis, enraiza no grupo externo e
     registra o conjunto S de atuais irmao do fossil. Pan-X = fracao das arvores
     com S contido em X (o fossil cai no grupo total de X).
  distancia  fracao dos splits do esqueleto molecular presentes na arvore ML
     morfologica (so atuais).

Saidas em runs/rapinantes/fosseis/arvores/: grupos.tsv, fosseis_posicao.tsv,
fosseis_pan.tsv, distancias.tsv.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import clados as cl
import dendropy

BASE = cl.BASE
ARV = BASE / "arvores"

EXTERNO = {
    "K19": "Eudromia_elegans",
    "K17": "Eudromia_elegans",
    "M11": "Tinamidae",
    "M05": "Tinamidae",
    "M16": "Anhimidae",
    "MC03": "Apsaravis",
    "LZ07": "Lithornis",
    "VZ26": "Anhimidae",
    "MK22": "Anhimidae",
    "C25": "Gansus_yumenensis",
    "D21": "Meleagris_gallopavo",
}
# Raiz entre os ATUAIS (quando o externo da matriz e fossil).
EXTERNO_ATUAL = {"MC03": "Rheidae", "LZ07": "Struthio", "C25": "Eudromia_elegans"}


def raiz_atual(id_: str) -> str:
    return EXTERNO_ATUAL.get(id_, EXTERNO[id_])


PAN = [
    "Cathartidae",
    "Ciconiidae",
    "Psittaciformes",
    "Passeriformes",
    "Psittacopasserae",
    "Falconiformes",
    "Falconidae",
    "Eufalconimorphae",
    "Cariamiformes",
    "Australaves",
    "Accipitriformes",
    "Accipitridae",
    "Strigiformes",
    "Hieraves",
    "Coraciimorphae",
    "Afroaves",
    "Telluraves",
]


def ler_arvores(caminho: Path, tns: dendropy.TaxonNamespace) -> list[dendropy.Tree]:
    return list(
        dendropy.TreeList.get(
            path=str(caminho),
            schema="newick",
            preserve_underscores=True,
            taxon_namespace=tns,
            rooting="force-unrooted",
        )
    )


def splits(arv: dendropy.Tree, manter: set[str], ref: str) -> set[frozenset[str]]:
    a = arv.extract_tree_with_taxa_labels(manter)
    out: set[frozenset[str]] = set()
    todos = frozenset(manter)
    rot: dict[int, frozenset[str]] = {}
    for no in a.postorder_node_iter():
        if no.is_leaf():
            rot[id(no)] = frozenset([no.taxon.label])
        else:
            rot[id(no)] = frozenset().union(*(rot[id(c)] for c in no.child_node_iter()))
        lbl = rot[id(no)]
        lado = lbl if ref not in lbl else todos - lbl
        if 1 < len(lado) < len(todos) - 1:
            out.add(lado)
    return out


def nome_de(s: frozenset[str], inv: dict[frozenset[str], str]) -> str:
    if s in inv:
        return inv[s]
    if len(s) <= 3:
        return "+".join(sorted(s))
    return f"[{len(s)} atuais: {','.join(sorted(s)[:3])},...]"


def irmao(arv: dendropy.Tree, fossil: str, atuais: set[str], externo: str) -> frozenset[str]:
    a = arv.extract_tree_with_taxa_labels(atuais | {fossil})
    og = a.find_node_with_taxon_label(externo)
    f = a.find_node_with_taxon_label(fossil)
    if og is None or f is None:
        raise ValueError(f"{fossil}/{externo} ausente da arvore")
    a.reroot_at_edge(og.edge, update_bipartitions=False)
    pai = f.parent_node
    if pai is None:
        raise ValueError(f"{fossil} sem pai depois do enraizamento")
    folhas = set()
    for c in pai.child_node_iter():
        if c is not f:
            folhas |= {lf.taxon.label for lf in c.leaf_iter()}
    return frozenset(folhas)


def resumir(id_: str) -> tuple[list, list, list, list]:
    tx = cl.ler_taxa(id_)
    atuais = {r["terminal"] for r in tx if r["tipo"] == "atual"}
    fosseis = [r["terminal"] for r in tx if r["tipo"] == "fossil"]
    cls = cl.clados(id_)
    inv: dict[frozenset[str], str] = {}
    prior = (
        cl.MOLECULARES
        + cl.TRADICIONAIS
        + ["Passeriformes", "Psittaciformes", "Falconiformes", "Cariamiformes"]
    )
    for k in sorted(cls, key=lambda k: (k not in prior, len(k))):
        inv.setdefault(cls[k], k)
    ref = raiz_atual(id_)
    g_rows, d_rows, pos_rows, pan_rows = [], [], [], []
    for modo in ("livre", "esqueleto"):
        tf, ub = ARV / id_ / f"{modo}.treefile", ARV / id_ / f"{modo}.ufboot"
        if not tf.exists() or not ub.exists():
            continue
        tns = dendropy.TaxonNamespace()
        ml = ler_arvores(tf, tns)[0]
        boots = ler_arvores(ub, tns)
        sp_ml = splits(ml, atuais, ref)
        if modo == "livre":
            sp_b = [splits(b, atuais, ref) for b in boots]
            for g in cl.TRADICIONAIS + cl.MOLECULARES:
                if g not in cls or len(cls[g]) < 2 or len(cls[g]) > len(atuais) - 2:
                    continue
                lado = cls[g] if ref not in cls[g] else frozenset(atuais) - cls[g]
                freq = sum(lado in s for s in sp_b) / len(sp_b)
                g_rows.append(
                    [
                        id_,
                        g,
                        "tradicional" if g in cl.TRADICIONAIS else "molecular",
                        len(cls[g]),
                        int(lado in sp_ml),
                        f"{100 * freq:.1f}",
                        len(boots),
                    ]
                )
            esq = BASE / "matrizes" / id_ / "esqueleto_stiller.nwk"
            if esq.exists():
                sp_e = splits(ler_arvores(esq, tns)[0], atuais, ref)
                rec = len(sp_e & sp_ml)
                d_rows.append([id_, len(atuais), len(sp_e), rec, f"{rec / max(len(sp_e), 1):.3f}"])
        for f in fosseis:
            cont = Counter(irmao(b, f, atuais, raiz_atual(id_)) for b in boots)
            s_ml = irmao(ml, f, atuais, raiz_atual(id_))
            for s, n in cont.most_common(3):
                pos_rows.append(
                    [
                        id_,
                        modo,
                        f,
                        nome_de(s_ml, inv),
                        nome_de(s, inv),
                        f"{100 * n / len(boots):.1f}",
                    ]
                )
            for x in PAN:
                if x in cls:
                    n = sum(c for s, c in cont.items() if s <= cls[x])
                    pan_rows.append(
                        [id_, modo, f, x, f"{100 * n / len(boots):.1f}", int(s_ml <= cls[x])]
                    )
    return g_rows, d_rows, pos_rows, pan_rows


def main(ids: list[str]) -> None:
    tabs = {
        "grupos.tsv": [
            "matriz",
            "grupo",
            "tipo",
            "n_terminais",
            "na_arvore_ML",
            "UFBoot_pct",
            "n_boot",
        ],
        "distancias.tsv": [
            "matriz",
            "n_atuais",
            "splits_esqueleto",
            "recuperados_na_ML_morfologica",
            "fracao",
        ],
        "fosseis_posicao.tsv": [
            "matriz",
            "analise",
            "fossil",
            "irmao_na_ML",
            "irmao_UFBoot",
            "UFBoot_pct",
        ],
        "fosseis_pan.tsv": ["matriz", "analise", "fossil", "grupo_total", "UFBoot_pct", "na_ML"],
    }
    acum: dict[str, list] = {k: [] for k in tabs}
    for id_ in ids:
        if not (ARV / id_).exists():
            continue
        g, d, p, pa = resumir(id_)
        acum["grupos.tsv"] += g
        acum["distancias.tsv"] += d
        acum["fosseis_posicao.tsv"] += p
        acum["fosseis_pan.tsv"] += pa
        print(id_, "ok", len(g), len(p))
    for nome, cab in tabs.items():
        with open(ARV / nome, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cab)
            w.writerows(acum[nome])


if __name__ == "__main__":
    main(sys.argv[1:] or ["K19", "M11", "M05", "M16", "MC03", "VZ26", "MK22", "C25", "D21", "LZ07"])
