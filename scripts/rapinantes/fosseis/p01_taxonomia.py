"""Passo 1: nomes para o taxonkit (WSL) e, na segunda chamada, o mapeamento
terminal morfologico -> especie representante na arvore de Stiller et al. 2024.

  uv run python p01_taxonomia.py nomes      # escreve taxonomia/nomes.txt
  (WSL) bash -l f02_taxonkit.sh             # escreve taxonomia/linhagem.tsv
  uv run python p01_taxonomia.py mapear     # escreve taxonomia/mapa_<ID>.tsv

Regra do mapeamento: terminal atual -> especie de Stiller da mesma familia
(NCBI); terminal de nivel ordinal -> da mesma ordem. Os casos ambiguos ficam
fixados a mao em FIXOS, com o motivo.
"""

from __future__ import annotations

import csv
import sys

import dendropy
import morfo_io as mi

SAIDA = mi.RAIZ / "runs" / "rapinantes" / "fosseis" / "taxonomia"
STILLER = mi.BRUTOS / "stiller_2024" / "63K.tre"

# Terminais cujo nome nao e um taxon NCBI direto, ou que agregam familias.
ALIAS = {
    "Strigidae_slash_Tytonidae": "Strigidae",
    "Coraciidae_slash_Brachypteraciidae": "Coraciidae",
    "Apodidae_slash_Hemiprocnidae": "Apodidae",
    "Upupidae_slash_Phoeniculidae": "Upupidae",
    "Cracidae_slash_Phasianidae": "Cracidae",
    "Todidae_slash_Momotidae": "Momotidae",
    "Pici": "Picidae",
    "Galbulae": "Galbulidae",
    "Pseudasturidae": None,
    "Catharthes": "Cathartes",
    "Polyborus": "Caracara",
    "Lophortyx": "Callipepla",
    "Porphyrula": "Porphyrio",
    "Afrotis": "Eupodotis",
    "Semornis": "Semnornis",
    "Megalaima": "Psilopogon",
    "Aulocorhynchus": "Aulacorhynchus",
    "Heteroscelus": "Tringa",
    "Ajaia": "Platalea",
    "Diatryma": None,
    # grafias da matriz de Degrange et al. 2021 (MorphoBank P3896)
    "Mirastur_gilvicollis": "Micrastur gilvicollis",
    "Ichtyophaga": None,  # "Ichthyophaga" no NCBI e um platelminto (homonimo)
    "Acrylium_vulturinum": "Acryllium vulturinum",
    "Ephippiorhinchus_senegalensis": "Ephippiorhynchus senegalensis",
    "Ephippiorhinchus_asiaticus": "Ephippiorhynchus asiaticus",
    "Haliaetus_leucocephalus": "Haliaeetus leucocephalus",
}

# Escolhas explicitas (terminal -> especie de Stiller), com motivo.
FIXOS = {
    # Strigiformes como um terminal: qualquer coruja serve; Strigidae escolhida.
    "Strigiformes": "Strix_occidentalis",
    "Psittaciformes": "Amazona_guildingii",
    "Passeriformes": "Corvus_cornix",
    # NCBI poe Agapornis (Psittaculidae na classificacao atual) em Psittacidae;
    # Amazona e Psittacidae s.s., como Psittacus.
    "Psittacus_erithacus": "Amazona_guildingii",
    "Psittacidae": "Amazona_guildingii",
    # Dacelo e Halcyoninae, como Halcyon; Ceyx e Alcedininae.
    "Dacelo_novaeguineae": "Halcyon_senegalensis",
    # NCBI ainda junta Calcarius/Emberiza em Fringillidae; Loxia e Carduelinae, como Spinus.
    "Spinus_pinus": "Loxia_curvirostra",
    # Buphagus virou Buphagidae; Sturnus e Sturnidae s.s., como Gracula.
    "Gracula_religiosa": "Sturnus_vulgaris",
    # Rhinopomastus (Phoeniculidae) esta em Stiller, mas o nome nao resolve no NCBI local.
    "Phoeniculus": "Rhinopomastus_cyanomelas",
    # Stiller nao tem Haematopodidae; Recurvirostridae e o parente mais proximo amostrado.
    "Haematopus": "Himantopus_himantopus",
    # D21 (Degrange et al. 2021). Accipitridae inteiro numa politomia: Stiller
    # tem so 5 acipitrideos (Aquilinae, Haliaeetinae, Circaetinae) e mapear
    # Elanus/Pernis/Gyps neles forcaria parentescos falsos dentro da familia.
    "Circaetus": "Aquila_chrysaetos",
    "Haliaetus_leucocephalus": "Aquila_chrysaetos",
    "Stephanoetus_coronatus": "Aquila_chrysaetos",
    "Ichtyophaga": "Aquila_chrysaetos",
    # Falconidae: Micrastur com Herpetotheres (Herpetotherinae); caracaras
    # (Polyborinae) no lugar de F. cherrug, irmaos de Falco (F. peregrinus).
    "Mirastur_gilvicollis": "Herpetotheres_cachinnans",
    "Falco_tinnunculus": "Falco_peregrinus",
    "Phalcoboenus_megalopterus": "Falco_cherrug",
    "Caracara_cherwayi": "Falco_cherrug",
    "Milvago_chimango": "Falco_cherrug",
    "Milvago_chimachima": "Falco_cherrug",
}


def nome_ncbi(terminal: str) -> str | None:
    if terminal in ALIAS:
        return ALIAS[terminal]
    partes = terminal.split("_")
    if len(partes) >= 2 and partes[1] and partes[1][0].islower() and partes[1] != "sp":
        return f"{partes[0]} {partes[1]}"
    return partes[0]


def especies_stiller() -> list[str]:
    t = dendropy.Tree.get(path=str(STILLER), schema="newick", preserve_underscores=True)
    return [lf.taxon.label for lf in t.leaf_node_iter()]


def cmd_nomes() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    nomes: set[str] = set()
    for sp in especies_stiller():
        g, e = sp.split("_")[:2]
        nomes.add(f"{g} {e}")
        nomes.add(g)
    for m in mi.todas().values():
        for t in m.taxa:
            if t in m.fosseis:
                continue
            n = nome_ncbi(t)
            if n:
                nomes.add(n)
                nomes.add(n.split()[0])
    (SAIDA / "nomes.txt").write_text("\n".join(sorted(nomes)) + "\n", encoding="utf-8")
    print(len(nomes), "nomes")


def ler_linhagem() -> dict[str, tuple[str, str]]:
    out: dict[str, tuple[str, str]] = {}
    with open(SAIDA / "linhagem.tsv", encoding="utf-8") as f:
        for lin in csv.reader(f, delimiter="\t"):
            if len(lin) < 4 or not lin[1]:
                continue
            nome, _taxid, ordem, familia = lin[0], lin[1], lin[2], lin[3]
            if nome not in out or (not out[nome][1] and familia):
                out[nome] = (ordem, familia)
    return out


def cmd_mapear() -> None:
    lin = ler_linhagem()
    sp = especies_stiller()
    por_fam: dict[str, list[str]] = {}
    por_ord: dict[str, list[str]] = {}
    for s in sp:
        g, e = s.split("_")[:2]
        o, f = lin.get(f"{g} {e}", lin.get(g, ("", "")))
        if f:
            por_fam.setdefault(f, []).append(s)
        if o:
            por_ord.setdefault(o, []).append(s)
    for id_, m in mi.todas().items():
        linhas = []
        for t in m.taxa:
            if t in m.fosseis:
                linhas.append((t, "fossil", "", "", "", ""))
                continue
            n = nome_ncbi(t)
            o, f = ("", "")
            if n:
                o, f = lin.get(n, lin.get(n.split()[0], ("", "")))
                if not f and n.endswith("idae"):
                    f = n
                if not o and n.endswith("iformes"):
                    o = n
            if t in FIXOS:
                if not f:  # nome fora do NCBI local: herda ordem/familia do representante
                    g_, e_ = FIXOS[t].split("_")[:2]
                    o, f = lin.get(f"{g_} {e_}", lin.get(g_, (o, f)))
                linhas.append((t, "atual", o, f, FIXOS[t], "fixo"))
                continue
            exato = n.replace(" ", "_") if n else ""
            mesmo_gen = sorted(s for s in sp if n and s.startswith(n.split()[0] + "_"))
            if exato in sp:
                rep, nivel = exato, "especie"
            elif mesmo_gen:
                rep, nivel = mesmo_gen[0], "genero"
            elif f and f in por_fam:
                cand = sorted(por_fam[f])
                gen = [c for c in cand if n and c.startswith(n.split()[0] + "_")]
                rep, nivel = (gen or cand)[0], ("genero" if gen else "familia")
            elif o and o in por_ord:
                rep, nivel = sorted(por_ord[o])[0], "ordem"
            else:
                rep, nivel = "", "SEM_MAPA"
            linhas.append((t, "atual", o, f, rep, nivel))
        with open(SAIDA / f"mapa_{id_}.tsv", "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["terminal", "tipo", "ordem_ncbi", "familia_ncbi", "rep_stiller", "nivel"])
            w.writerows(linhas)
        sem = [x[0] for x in linhas if x[5] == "SEM_MAPA"]
        ordem = [x[0] for x in linhas if x[5] == "ordem"]
        reps = [x[4] for x in linhas if x[4]]
        dup = sorted({r for r in reps if reps.count(r) > 1})
        print(f"{id_}: sem_mapa={sem} nivel_ordem={ordem} reps_repetidos={dup}")


if __name__ == "__main__":
    {"nomes": cmd_nomes, "mapear": cmd_mapear}[sys.argv[1]]()
