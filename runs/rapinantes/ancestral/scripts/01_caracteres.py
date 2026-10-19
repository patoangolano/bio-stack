# 01_caracteres.py
"""Matriz de caracteres por família (pontas da árvore de Stiller et al. 2024).

Entradas (todas em runs/rapinantes/ancestral/):
  downloads/erda/07_scripts_figures__Fig1__annotations.txt   espécie -> família/ordem (Stiller 2024)
  downloads/avonet/AVONET2_eBird.csv,
      crosswalk_eBird.csv, crosswalk_BirdTree.csv              AVONET (Tobias et al. 2022), figshare 16586228 v7
  downloads/eltontraits/BirdFuncDat.txt                       EltonTraits 1.0 (Wilman et al. 2014), esapubs E095-178
  caracteres/morfologia_codificacao_manual.tsv                bico e pés (codificação manual com fonte)
  caracteres/fosseis_codificacao_manual.tsv                   pontas fósseis

Saídas:
  caracteres/familias_pontas.tsv        família -> espécie representante, famílias Clements agregadas
  caracteres/caracteres_familia.tsv     proporções e codificações, com a fonte de cada célula
  caracteres/matriz_mk.tsv              só as codificações (0, 1, 0&1) por ponta, inclusive fósseis

Execução (do Windows, raiz do projeto):
  uv run python -I runs/rapinantes/ancestral/scripts/01_caracteres.py
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]
DL = BASE / "downloads"
OUT = BASE / "caracteres"

CLADES_INCLUIDOS = {"Telluraves", "Strisores", "Aequornithes", "Galloanseres", "Columbimorphae", "Otidimorphae"}
ORDENS_EXTRA = {"Charadriiformes"}

# Limiares (inferior, superior): p >= sup -> 1; inf <= p < sup -> 0&1; p < inf -> 0
# sinônimos que nenhuma tabela de correspondência do AVONET resolve
SINONIMOS = {"Rhadina_sibilatrix": "Phylloscopus_sibilatrix"}

LIMIARES = {"T50": (0.10, 0.50), "T25": (0.05, 0.25), "T75": (0.25, 0.75)}


def codifica(p: float, inf: float, sup: float) -> str:
    if pd.isna(p):
        return "0&1"
    if p >= sup:
        return "1"
    if p >= inf:
        return "0&1"
    return "0"


def estados(c: str) -> set[int]:
    c = str(c).strip()
    if c in {"?", "0&1", "nan", ""}:
        return {0, 1}
    return {int(x) for x in c.split("&")}


def sindrome(trof: str, bico: str, pes: str) -> str:
    """0 = nenhum componente; 1 = um ou dois; 2 = os três."""
    soma = {sum(t) for t in itertools.product(estados(trof), estados(bico), estados(pes))}
    s = sorted({0 if x == 0 else (2 if x == 3 else 1) for x in soma})
    return "&".join(str(x) for x in s)


def main() -> None:
    ann = pd.read_csv(DL / "erda/07_scripts_figures__Fig1__annotations.txt", sep="\t")
    ann.columns = [c.strip() for c in ann.columns]
    ann = ann[ann["clade_mag7"].isin(CLADES_INCLUIDOS) | ann["order"].isin(ORDENS_EXTRA)].copy()

    eb = pd.read_csv(DL / "avonet/AVONET2_eBird.csv")
    cw_eb = pd.read_csv(DL / "avonet/crosswalk_eBird.csv")
    cw_bt = pd.read_csv(DL / "avonet/crosswalk_BirdTree.csv")
    eb["sp"] = eb["Species2"].str.replace(" ", "_")
    fam_eb = dict(zip(eb["sp"], eb["Family2"], strict=True))
    bl2eb: dict[str, list[str]] = {}
    for a, b in zip(cw_eb["Species1"], cw_eb["Species2"], strict=True):
        if isinstance(a, str) and isinstance(b, str):
            bl2eb.setdefault(a.replace(" ", "_"), []).append(b.replace(" ", "_"))
    bt2bl: dict[str, list[str]] = {}
    for a, b in zip(cw_bt["Species1"], cw_bt["Species3"], strict=True):
        if isinstance(a, str) and isinstance(b, str):
            bt2bl.setdefault(b.replace(" ", "_"), []).append(a.replace(" ", "_"))
    genero_fam = eb.assign(gen=eb["sp"].str.split("_").str[0]).groupby("gen")["Family2"].agg(
        lambda s: s.value_counts().index[0]
    )

    # 1. espécie de Stiller -> família Clements (AVONET eBird)
    linhas = []
    for sp0 in ann["taxa"]:
        sp = SINONIMOS.get(sp0, sp0)
        if sp in fam_eb:
            linhas.append((sp0, fam_eb[sp], "nome exato eBird" if sp == sp0 else "sinônimo manual"))
            continue
        alvo = [e for e in bl2eb.get(sp, []) if e in fam_eb]
        if alvo:
            linhas.append((sp0, fam_eb[alvo[0]], "BirdLife->eBird"))
            continue
        alvo = [e for b in bt2bl.get(sp, []) for e in bl2eb.get(b, []) if e in fam_eb]
        if alvo:
            linhas.append((sp0, fam_eb[alvo[0]], "BirdTree->BirdLife->eBird"))
            continue
        gen = sp.split("_")[0]
        linhas.append((sp0, genero_fam.get(gen, pd.NA), "gênero"))
    mapa = pd.DataFrame(linhas, columns=["taxa", "fam_clements", "metodo"])
    ann = ann.merge(mapa, on="taxa", how="left")
    mapa.to_csv(OUT / "pareamento_especies_stiller_avonet.tsv", sep="\t", index=False)
    if ann["fam_clements"].isna().any():
        raise SystemExit(f"espécies sem família Clements: {ann.loc[ann.fam_clements.isna(), 'taxa'].tolist()}")

    # 2. família de Stiller -> conjunto de famílias Clements; e o inverso
    st2cl = ann.groupby("family")["fam_clements"].agg(lambda s: sorted(set(s)))
    cl2st: dict[str, list[str]] = {}
    for f, cls in st2cl.items():
        for c in cls:
            cl2st.setdefault(c, []).append(f)
    compartilhadas = {c: v for c, v in cl2st.items() if len(v) > 1}

    reps = ann.sort_values("taxa").groupby("family").first().reset_index()
    fams = reps[["family", "order", "clade_mag7", "taxa"]].rename(columns={"taxa": "ponta"})
    fams["especies_stiller"] = fams["family"].map(ann.groupby("family")["taxa"].agg(lambda s: ";".join(sorted(s))))
    fams["familias_clements"] = fams["family"].map(lambda f: ";".join(st2cl[f]))
    fams["clements_compartilhada"] = fams["family"].map(
        lambda f: ";".join(c for c in st2cl[f] if c in compartilhadas)
    )
    fams.to_csv(OUT / "familias_pontas.tsv", sep="\t", index=False)

    # 3. AVONET: proporções por família de Stiller
    eb_tn = eb[["sp", "Family2", "Trophic.Niche"]].copy()
    regs = []
    for f in fams["family"]:
        pool = eb_tn[eb_tn["Family2"].isin(st2cl[f])]
        tn = pool["Trophic.Niche"]
        n = int(tn.notna().sum())
        nv, ns, na = int((tn == "Vertivore").sum()), int((tn == "Scavenger").sum()), int((tn == "Aquatic predator").sum())
        regs.append(
            {
                "family": f,
                "AV_n": n,
                "AV_vertivore": nv,
                "AV_scavenger": ns,
                "AV_aquatic_predator": na,
                "p_VS": (nv + ns) / n if n else float("nan"),
                "p_VSA": (nv + ns + na) / n if n else float("nan"),
            }
        )
    av = pd.DataFrame(regs)

    # 4. EltonTraits: % da dieta em vertebrados terrestres + carniça, por espécie
    et = pd.read_csv(DL / "eltontraits/BirdFuncDat.txt", sep="\t", encoding="latin-1")
    et = et[et["Scientific"].notna()].copy()
    et["sp"] = et["Scientific"].str.strip().str.replace(" ", "_")
    et["pct_vert_scav"] = et[["Diet-Vend", "Diet-Vect", "Diet-Vunk", "Diet-Scav"]].sum(axis=1)

    def fam_et(sp: str) -> str | None:
        cands = [e for b in bt2bl.get(sp, [sp]) for e in bl2eb.get(b, []) if e in fam_eb]
        if not cands and sp in fam_eb:
            cands = [sp]
        return fam_eb[cands[0]] if cands else None

    et["fam_clements"] = et["sp"].map(fam_et)
    et_ok = et[et["fam_clements"].notna()]
    regs = []
    for f in fams["family"]:
        pool = et_ok[et_ok["fam_clements"].isin(st2cl[f])]
        n = len(pool)
        k = int((pool["pct_vert_scav"] >= 50).sum())
        regs.append({"family": f, "ET_n": n, "ET_vert_scav_ge50": k, "p_ET": k / n if n else float("nan")})
    etf = pd.DataFrame(regs)

    tab = fams.merge(av, on="family").merge(etf, on="family")

    # 5. codificações tróficas
    for nome, (inf, sup) in LIMIARES.items():
        tab[f"C1_{nome}"] = [codifica(p, inf, sup) for p in tab["p_VS"]]
    tab["C1_bin50"] = ["1" if p >= 0.5 else "0" for p in tab["p_VS"]]
    tab["C2_aquat_T50"] = [codifica(p, *LIMIARES["T50"]) for p in tab["p_VSA"]]
    tab["C3_elton_T50"] = [codifica(p, *LIMIARES["T50"]) for p in tab["p_ET"]]

    def fonte_av(r: pd.Series, col: str, lim: str) -> str:
        return (
            f"AVONET2_eBird Trophic.Niche [{r['familias_clements']}]: n={r['AV_n']}, "
            f"Vertivore={r['AV_vertivore']}, Scavenger={r['AV_scavenger']}, "
            f"Aquatic predator={r['AV_aquatic_predator']}; p={r[col]:.3f}; limiar {lim}"
        )

    tab["fonte_C1"] = [fonte_av(r, "p_VS", "T50/T25/T75/bin50") for _, r in tab.iterrows()]
    tab["fonte_C2"] = [fonte_av(r, "p_VSA", "T50") for _, r in tab.iterrows()]
    tab["fonte_C3"] = [
        f"EltonTraits 1.0 [{r['familias_clements']}]: n={r['ET_n']}, spp com Vend+Vect+Vunk+Scav>=50%: "
        f"{r['ET_vert_scav_ge50']}; p={r['p_ET']:.3f}; limiar T50"
        for _, r in tab.iterrows()
    ]

    # 6. morfologia
    morf = pd.read_csv(OUT / "morfologia_codificacao_manual.tsv", sep="\t", dtype=str)
    faltam = set(morf["familia"]) - set(tab["family"])
    if faltam:
        raise SystemExit(f"famílias da tabela manual ausentes da árvore: {sorted(faltam)}")
    morf = morf.set_index("familia")
    tab["C4_bico_estrito"] = tab["family"].map(morf["bico_estrito"]).fillna("0")
    tab["C4b_bico_amplo"] = tab["family"].map(morf["bico_amplo"]).fillna("0")
    tab["C5_pes"] = tab["family"].map(morf["pes"]).fillna("0")
    tab["fonte_C4"] = tab["family"].map(morf["fonte_bico"]).fillna("ANALISTA (padrão: sem gancho rapinante)")
    tab["fonte_C5"] = tab["family"].map(morf["fonte_pes"]).fillna("ANALISTA (padrão: sem pé rapinante)")
    tab["nota_morfologia"] = tab["family"].map(morf["nota"]).fillna("")
    tab["C6_sindrome"] = [
        sindrome(a, b, c) for a, b, c in zip(tab["C1_T50"], tab["C4_bico_estrito"], tab["C5_pes"], strict=True)
    ]
    tab["fonte_C6"] = "derivado: C1_T50 + C4_bico_estrito + C5_pes (0 = nenhum, 1 = 1-2, 2 = os 3)"
    tab.to_csv(OUT / "caracteres_familia.tsv", sep="\t", index=False, float_format="%.4f")

    # 7. fósseis e matriz final para o Mk
    fos = pd.read_csv(OUT / "fosseis_codificacao_manual.tsv", sep="\t", dtype=str)
    cols = ["C1_T50", "C1_T25", "C1_T75", "C1_bin50", "C2_aquat_T50", "C3_elton_T50",
            "C4_bico_estrito", "C4b_bico_amplo", "C5_pes", "C6_sindrome"]
    mat = tab[["ponta", "family"] + cols].copy()
    mat["tipo"] = "vivente"
    frows = []
    for _, r in fos.iterrows():
        t = "0&1" if r["trofico"] == "?" else r["trofico"]
        row = {"ponta": r["ponta"], "family": r["grupo"], "tipo": "fossil"}
        for c in cols[:6]:
            row[c] = t
        row["C4_bico_estrito"] = "0&1" if r["bico_estrito"] == "?" else r["bico_estrito"]
        row["C4b_bico_amplo"] = "0&1" if r["bico_amplo"] == "?" else r["bico_amplo"]
        row["C5_pes"] = "0&1" if r["pes"] == "?" else r["pes"]
        row["C6_sindrome"] = sindrome(row["C1_T50"], row["C4_bico_estrito"], row["C5_pes"])
        frows.append(row)
    mat = pd.concat([mat, pd.DataFrame(frows)], ignore_index=True)
    # sensibilidade: dieta dos fósseis desconhecida (0&1) -- separa o efeito da posição/idade
    # dos fósseis do efeito da dieta inferida da morfologia (que é parcialmente circular)
    fosq = mat["tipo"] == "fossil"
    for c in ["C1_T50", "C3_elton_T50"]:
        mat[f"{c}_fosQ"] = mat[c].where(~fosq, "0&1")
    mat["C6_sindrome_fosQ"] = [
        sindrome(t, b, p) if fo else s
        for t, b, p, s, fo in zip(
            mat["C1_T50_fosQ"], mat["C4_bico_estrito"], mat["C5_pes"], mat["C6_sindrome"], fosq, strict=True
        )
    ]
    mat.to_csv(OUT / "matriz_mk.tsv", sep="\t", index=False)

    # resumo no console
    print(f"pontas viventes: {len(tab)}; fósseis: {len(frows)}")
    print("métodos de pareamento de nome:", mapa["metodo"].value_counts().to_dict())
    print("famílias Clements compartilhadas por >1 família de Stiller:", compartilhadas)
    for c in cols:
        print(c, mat[mat.tipo == "vivente"][c].value_counts().to_dict())
    show = tab[(tab["p_VS"] > 0) | (tab["C4_bico_estrito"] != "0") | (tab["C5_pes"] != "0")]
    print(show[["family", "AV_n", "p_VS", "p_VSA", "p_ET", "C1_T50", "C4_bico_estrito", "C5_pes", "C6_sindrome"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
