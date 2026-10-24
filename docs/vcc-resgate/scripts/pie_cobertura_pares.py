# pie_cobertura_pares.py
# (a) cobertura dos 300 alvos do VCC 2026 nas fontes de conhecimento do PIE_sources
# (b) perturbacoes medidas em dois contextos nos dados PIE (validacao deixa-um-contexto-de-fora)
# Uso: python -I pie_cobertura_pares.py <raiz do projeto>
import itertools
import json
import sys
from pathlib import Path

import numpy as np

raiz = Path(sys.argv[1])
pie = raiz / "data/raw/pie"
alvos = [
    linha.strip()
    for linha in (raiz / "data/raw/vcc/controls-final/pert_counts.csv").read_text(encoding="utf-8").splitlines()[1:]
    if linha.strip()
]
assert len(alvos) == 300 and len(set(alvos)) == 300, (len(alvos), len(set(alvos)))
alvos_set = set(alvos)


def ler_aliases(caminho: Path) -> dict[str, str]:
    # "CHAVE: VALOR" por linha; sem yaml.load para nao executar nada do arquivo baixado
    mapa = {}
    if not caminho.exists():
        return mapa
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        linha = linha.split("#", 1)[0].strip()
        if not linha:
            continue
        k, v = linha.split(":", 1)
        mapa[k.strip()] = v.strip()
    return mapa


# ---------- (a) ----------
fontes = ["esm2", "string_space", "depmap_gene_effect", "ncbi_text", "perturbation_text", "gene_text"]
cobertura = {}
print("== (a) cobertura por fonte")
for f in fontes:
    meta = json.loads((pie / "PIE_sources" / f / "meta.json").read_text(encoding="utf-8"))
    chaves = set(meta["keys"])
    assert len(chaves) == len(meta["keys"]), f"chaves duplicadas em {f}"
    aliases = ler_aliases(pie / "PIE_sources" / f / "aliases.yaml")
    diretos = [a for a in alvos if a in chaves]
    via_alias = [a for a in alvos if a not in chaves and a in aliases and aliases[a] in chaves]
    faltam = [a for a in alvos if a not in chaves and not (a in aliases and aliases[a] in chaves)]
    # alias que apontaria para um alvo mas cujo destino nao existe
    alias_quebrado = [a for a in alvos if a not in chaves and a in aliases and aliases[a] not in chaves]
    cobertura[f] = set(diretos) | set(via_alias)
    print(f"{f}\tindex={meta['index']}\tn_chaves={len(chaves)}\tn_aliases={len(aliases)}"
          f"\tdiretos={len(diretos)}\tvia_alias={len(via_alias)} {via_alias}\ttotal={len(cobertura[f])}"
          f"\tfaltam={len(faltam)} {sorted(faltam)}\talias_quebrado={alias_quebrado}")

todas = set.intersection(*cobertura.values())
print(f"em TODAS as 6 fontes: {len(todas)}")
print(f"faltam em pelo menos uma: {sorted(alvos_set - todas)}")
cinco = set.intersection(*(cobertura[f] for f in fontes if f != "perturbation_text"))
print(f"em todas exceto perturbation_text: {len(cinco)}")
print(f"faltam em pelo menos uma das 5 (sem perturbation_text): {sorted(alvos_set - cinco)}")

# ---------- (b) ----------
print("\n== (b) pares de contextos")
alias_uniao = {}
for f in ("esm2", "ncbi_text", "depmap_gene_effect"):
    alias_uniao.update(ler_aliases(pie / "PIE_sources" / f / "aliases.yaml"))


def medidas_por_contexto(repo: str) -> dict[str, set[str]]:
    meta = json.loads((pie / repo / "meta.json").read_text(encoding="utf-8"))
    assert meta["controls_only"] is False, repo
    ctx = np.load(pie / repo / "ctx_ids.npy", allow_pickle=False)
    pert = np.load(pie / repo / "pert_ids.npy", allow_pickle=False)
    assert ctx.shape == pert.shape == (meta["num_rows"],), repo
    id2ctx = {v: k for k, v in meta["context_to_id"].items()}
    id2pert = {v: k for k, v in meta["pert_to_id"].items()}
    linhas = set(zip(ctx.tolist(), pert.tolist(), strict=True))
    assert len(linhas) == len(ctx), f"(ctx, pert) repetido em {repo}"
    saida: dict[str, set[str]] = {}
    for c, p in linhas:
        saida.setdefault(id2ctx[c], set()).add(id2pert[p])
    return saida


orion = medidas_por_contexto("PIE_x_atlas_orion")
arc = medidas_por_contexto("PIE_arc_vcc_25")
rep = medidas_por_contexto("PIE_replogle_nadig_essential")
for nome, d in (("orion", orion), ("arc_vcc_25", arc), ("replogle", rep)):
    for c, s in sorted(d.items()):
        print(f"{nome}\t{c}\tn_perts={len(s)}\talvos_vcc={len(s & alvos_set)}")

h1_uniao = arc["ARC_H1"] | arc["ARC_H1_VAL"]


def alvos_em(s: set[str]) -> tuple[list[str], list[str]]:
    # casamento direto; depois por alias (alvo antigo -> simbolo novo) nos dois sentidos
    diretos = s & alvos_set
    inv = {v: k for k, v in alias_uniao.items()}
    via = {a for a in alvos_set - diretos if alias_uniao.get(a) in s or inv.get(a) in s}
    return sorted(diretos), sorted(via)


pares = [
    ("PIE_x_atlas_orion", "hct116", orion["hct116"], "PIE_x_atlas_orion", "hek293t", orion["hek293t"]),
    ("PIE_arc_vcc_25", "ARC_H1", arc["ARC_H1"], "PIE_x_atlas_orion", "hct116", orion["hct116"]),
    ("PIE_arc_vcc_25", "ARC_H1_VAL", arc["ARC_H1_VAL"], "PIE_x_atlas_orion", "hct116", orion["hct116"]),
    ("PIE_arc_vcc_25", "ARC_H1+ARC_H1_VAL", h1_uniao, "PIE_x_atlas_orion", "hct116", orion["hct116"]),
    ("PIE_arc_vcc_25", "ARC_H1+ARC_H1_VAL", h1_uniao, "PIE_x_atlas_orion", "hek293t", orion["hek293t"]),
    ("PIE_arc_vcc_25", "ARC_H1", arc["ARC_H1"], "PIE_arc_vcc_25", "ARC_H1_VAL", arc["ARC_H1_VAL"]),
]
for c1, c2 in itertools.combinations(sorted(rep), 2):
    pares.append(("PIE_replogle_nadig_essential", c1, rep[c1], "PIE_replogle_nadig_essential", c2, rep[c2]))

todas4 = set.intersection(*rep.values())
saida = raiz / "runs/vcc/caracterizacao/pares_validacao.tsv"
cab = ["repo_a", "contexto_a", "n_a", "repo_b", "contexto_b", "n_b", "n_comum",
       "n_comum_alvos_vcc", "n_comum_alvos_vcc_via_alias", "alvos_vcc_comuns"]
linhas_tsv = ["\t".join(cab)]
for ra, ca, sa, rb, cb, sb in pares:
    comum = sa & sb
    d, v = alvos_em(comum)
    linhas_tsv.append("\t".join(map(str, [ra, ca, len(sa), rb, cb, len(sb), len(comum), len(d), len(v),
                                          ",".join(d + [f"{x}(alias)" for x in v])])))
d4, v4 = alvos_em(todas4)
linhas_tsv.append("\t".join(map(str, ["PIE_replogle_nadig_essential", "hepg2&jurkat&k562&rpe1", len(todas4),
                                      "-", "-", "-", len(todas4), len(d4), len(v4), ",".join(d4 + v4)])))
saida.write_text("\n".join(linhas_tsv) + "\n", encoding="utf-8")
for linha_tsv in linhas_tsv:
    print(linha_tsv)
print(f"gravado: {saida}")

# sanidade contra fatos ja medidos: uniao Orion com alvos
print("orion uniao alvos:", len((orion["hct116"] | orion["hek293t"]) & alvos_set),
      "hct116:", len(orion["hct116"] & alvos_set), "hek293t:", len(orion["hek293t"] & alvos_set))
print("H1 uniao alvos:", len(h1_uniao & alvos_set), "replogle uniao alvos:", len(set.union(*rep.values()) & alvos_set))
# alvos que aparecem em PIE so via alias
todos_pie = orion["hct116"] | orion["hek293t"] | h1_uniao | set.union(*rep.values())
print("alvos em algum PIE (direto):", len(todos_pie & alvos_set), "via alias extra:", alvos_em(todos_pie)[1])
