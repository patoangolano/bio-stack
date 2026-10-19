"""Passo 2: matrizes prontas para o IQ-TREE, particoes por numero de estados,
arvore-esqueleto molecular (Stiller et al. 2024, 63K ASTRAL) e metadados.

Decisoes, todas registradas em matrizes/<ID>/decisoes.json:
- polimorfismo e '-' viram '?': o Mk do IQ-TREE nao modela polimorfismo.
- estados de cada carater sao renumerados 0..k-1 pela ordem; as particoes
  agrupam caracteres por k (2, 3, >=4), cada uma com MK+ASC+G4 (Lewis 2001:
  Mk precisa de k certo; +ASC porque matrizes morfologicas nao trazem os
  caracteres invariantes). Particao com < 3 caracteres e fundida a vizinha
  de k menor (o k da fusao e o maior).
- caracteres sem variacao depois da recodificacao saem (o +ASC exige).
- carater ordenado pelos autores entra nao ordenado (o Mk e nao ordenado);
  quantos eram fica no decisoes.json.
- esqueleto: a arvore de Stiller podada aos representantes; terminal
  morfologico que divide representante com outro vira politomia (familia
  sem resolucao interna). Fosseis ficam fora do esqueleto, livres.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys

import dendropy
import morfo_io as mi

BASE = mi.RAIZ / "runs" / "rapinantes" / "fosseis"
TAXO = BASE / "taxonomia"

ACCIPITRIFORMES = {"Accipitridae", "Pandionidae", "Sagittariidae", "Cathartidae"}
CORACIIMORPHAE = {
    "Coliiformes",
    "Leptosomiformes",
    "Trogoniformes",
    "Bucerotiformes",
    "Coraciiformes",
    "Piciformes",
}
RAPINANTE_FAM = {
    "Accipitridae",
    "Pandionidae",
    "Sagittariidae",
    "Falconidae",
    "Strigidae",
    "Tytonidae",
}
# terminais de nivel familiar/ordinal que nao trazem familia NCBI na linhagem
RAPINANTE_TERM = {"Strigiformes", "Strigidae_slash_Tytonidae"}
BICO_GANCHO_EXTRA = {"Cathartidae", "Cariamidae"}


def clado_molecular(ordem: str, familia: str, terminal: str) -> str:
    if familia in ACCIPITRIFORMES or terminal in {"Accipitridae", "Sagittariidae", "Cathartidae"}:
        return "Accipitriformes"
    if ordem in CORACIIMORPHAE:
        return "Coraciimorphae"
    if terminal in {
        "Pici",
        "Galbulae",
        "Upupidae_slash_Phoeniculidae",
        "Coraciidae_slash_Brachypteraciidae",
    }:
        return "Coraciimorphae"
    if terminal in RAPINANTE_TERM:
        return "Strigiformes"
    if ordem in {
        "Strigiformes",
        "Falconiformes",
        "Cariamiformes",
        "Psittaciformes",
        "Passeriformes",
    }:
        return ordem
    if terminal in {"Psittaciformes", "Passeriformes", "Psittacidae"}:
        return {"Psittacidae": "Psittaciformes"}.get(terminal, terminal)
    return "outros"


def guilda(familia: str, terminal: str) -> tuple[str, str]:
    rap = familia in RAPINANTE_FAM or terminal in RAPINANTE_TERM or terminal in RAPINANTE_FAM
    gancho = (
        rap
        or familia in BICO_GANCHO_EXTRA
        or terminal in BICO_GANCHO_EXTRA
        or terminal == "Cariama"
    )
    return ("rapinante" if rap else "nao_rapinante", "bico_gancho" if gancho else "outro")


def ler_mapa(id_: str) -> list[dict]:
    with open(TAXO / f"mapa_{id_}.tsv", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def preparar(id_: str, m: mi.Matriz) -> None:
    """Escreve em pasta temporaria e so substitui o arquivo que mudou: um
    IQ-TREE rodando em paralelo nunca ve um .phy reescrito pela metade."""
    destino = BASE / "matrizes" / id_
    destino.mkdir(parents=True, exist_ok=True)
    pasta = BASE / "matrizes" / f".{id_}.tmp"
    pasta.mkdir(parents=True, exist_ok=True)
    _preparar(id_, m, pasta)
    for arq in pasta.iterdir():
        alvo = destino / arq.name
        if not alvo.exists() or alvo.read_bytes() != arq.read_bytes():
            shutil.copyfile(arq, alvo)
            print("  atualizado:", alvo.name)
    shutil.rmtree(pasta)


def _preparar(id_: str, m: mi.Matriz, pasta) -> None:
    # --- recodificacao e particoes ---
    mantidos, removidos, recod = [], [], {}
    for j in range(m.nchar):
        est = sorted(m.estados_sem_polimorfismo(j))
        if len(est) < 2:
            removidos.append(j)
            continue
        mantidos.append(j)
        recod[j] = {e: str(i) for i, e in enumerate(est)}
    k_de = {j: len(recod[j]) for j in mantidos}
    classes: dict[str, list[int]] = {"k2": [], "k3": [], "k4mais": []}
    for j in mantidos:
        classes["k2" if k_de[j] == 2 else "k3" if k_de[j] == 3 else "k4mais"].append(j)
    if 0 < len(classes["k4mais"]) < 3:
        classes["k3"] += classes.pop("k4mais")
    if 0 < len(classes["k3"]) < 3:
        classes["k2"] += classes.pop("k3")
    classes = {c: sorted(v) for c, v in classes.items() if v}
    ordem_colunas = [j for c in classes for j in classes[c]]
    # --- alinhamento (PHYLIP relaxado) ---
    seqs = {}
    for t in m.taxa:
        s = []
        for j in ordem_colunas:
            c = m.linhas[t][j]
            s.append(recod[j][next(iter(c))] if c is not None and len(c) == 1 else "?")
        seqs[t] = "".join(s)
    with open(pasta / f"{id_}.phy", "w", encoding="utf-8", newline="\n") as f:
        f.write(f"{len(m.taxa)} {len(ordem_colunas)}\n")
        for t in m.taxa:
            f.write(f"{t}  {seqs[t]}\n")
    # --- particoes ---
    linhas = ["#nexus", "begin sets;"]
    pos, partes = 1, []
    for c, cols in classes.items():
        linhas.append(f"  charset {c} = {pos}-{pos + len(cols) - 1};")
        partes.append(f"MK+ASC+G4:{c}")
        pos += len(cols)
    linhas.append(f"  charpartition morfo = {', '.join(partes)};")
    linhas.append("end;")
    (pasta / "particoes.nex").write_text("\n".join(linhas) + "\n", encoding="utf-8")
    # coluna do alinhamento -> carater original (1-based)
    with open(pasta / "colunas.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["coluna", "carater_original", "k", "particao", "ordenado_na_fonte", "descricao"]
        )
        col = 1
        for c, cols in classes.items():
            for j in cols:
                w.writerow([col, j + 1, k_de[j], c, int(j in m.ordenados), m.descricoes.get(j, "")])
                col += 1
    # --- metadados de taxa + esqueleto ---
    mapa = {r["terminal"]: r for r in ler_mapa(id_)}
    with open(pasta / "taxa.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "terminal",
                "tipo",
                "ordem",
                "familia",
                "rep_stiller",
                "clado",
                "guilda",
                "guilda_bico_gancho",
                "pct_codificado",
            ]
        )
        for t in m.taxa:
            r = mapa[t]
            pct = 100 * sum(ch != "?" for ch in seqs[t]) / len(ordem_colunas)
            if t in m.fosseis:
                w.writerow([t, "fossil", "", "", "", "fossil", "fossil", "fossil", f"{pct:.1f}"])
                continue
            g1, g2 = guilda(r["familia_ncbi"], t)
            w.writerow(
                [
                    t,
                    "atual",
                    r["ordem_ncbi"],
                    r["familia_ncbi"],
                    r["rep_stiller"],
                    clado_molecular(r["ordem_ncbi"], r["familia_ncbi"], t),
                    g1,
                    g2,
                    f"{pct:.1f}",
                ]
            )
    esqueleto(id_, m, mapa, pasta)
    dec = {
        "fonte": m.fonte,
        "ntax": len(m.taxa),
        "nchar_fonte": m.nchar,
        "nchar_usados": len(ordem_colunas),
        "removidos_invariantes_apos_recodificacao": [j + 1 for j in removidos],
        "polimorfismos_para_missing": sum(
            1 for t in m.taxa for c in m.linhas[t] if c is not None and len(c) > 1
        ),
        "ordenados_na_fonte_tratados_como_nao_ordenados": len(m.ordenados & set(mantidos)),
        "particoes": {c: len(v) for c, v in classes.items()},
        "fosseis": sorted(m.fosseis),
        "externo": m.externo,
    }
    (pasta / "decisoes.json").write_text(
        json.dumps(dec, indent=1, ensure_ascii=False), encoding="utf-8"
    )
    print(
        id_, json.dumps({k: dec[k] for k in ("ntax", "nchar_fonte", "nchar_usados", "particoes")})
    )


def esqueleto(id_: str, m: mi.Matriz, mapa: dict, pasta) -> None:
    arv = dendropy.Tree.get(
        path=str(mi.BRUTOS / "stiller_2024" / "63K.tre"), schema="newick", preserve_underscores=True
    )
    por_rep: dict[str, list[str]] = {}
    for t in m.taxa:
        if t in m.fosseis:
            continue
        por_rep.setdefault(mapa[t]["rep_stiller"], []).append(t)
    arv.retain_taxa_with_labels(set(por_rep))
    arv.suppress_unifurcations()
    for no in arv.leaf_node_iter():
        terms = por_rep[no.taxon.label]
        if len(terms) == 1:
            no.taxon.label = terms[0]
        else:
            no.taxon = None
            for t in terms:
                filho = dendropy.Node(taxon=dendropy.Taxon(label=t))
                no.add_child(filho)
    for no in arv.preorder_node_iter():
        no.label = None
        no.edge.length = None
    nwk = arv.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True).strip()
    (pasta / "esqueleto_stiller.nwk").write_text(nwk + "\n", encoding="utf-8")


if __name__ == "__main__":
    ids = sys.argv[1:] or None
    for id_, m in mi.todas().items():
        if ids is None or id_ in ids:
            preparar(id_, m)
