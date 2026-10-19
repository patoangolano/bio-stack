"""
Tabela-sintese: nos-chave das aves de rapina x fontes (B10K quartetos, arvores
publicadas do B10K, Open Tree 16.1, Aves 1.6, arvores-fonte da Open Tree,
TimeTree). Le os TSV das outras etapas; nao recalcula nada.

    uv run python scripts/rapinantes/reanalise/tabela_sintese.py

Saidas: runs/rapinantes/reanalise/tabela_sintese.tsv (todas as colunas) e
tabela_sintese.md (versao compacta para o texto).
"""

from __future__ import annotations

import csv
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]
R = RAIZ / "runs" / "rapinantes" / "reanalise"


def tsv(caminho: Path) -> list[dict[str, str]]:
    if not caminho.exists():
        return []
    return list(csv.DictReader(open(caminho, encoding="utf-8"), delimiter="\t"))


# no-chave -> (pergunta B10K, nos Open Tree, pergunta votos OTL, par TimeTree, rotulo das 3 resolucoes)
NOS = [
    ("Base de Afroaves: Accipitrimorphae+Strigiformes | Strigiformes+Coraciimorphae | Accipitrimorphae+Coraciimorphae",
     "afroaves_base", ["Accipitrimorphae+Strigiformes", "Strigiformes+Coraciimorphae",
                       "Accipitrimorphae+Coraciimorphae"], "afroaves_base", "Strigiformes-Accipitriformes"),
    ("Monofilia de Afroaves | Accipitrimorphae (ou Hieraves) irma do resto de Telluraves",
     "afroaves_mono_se_AS", ["Afroaves (Accipitrimorphae+Strigiformes+Coraciimorphae)"], "afroaves_mono", None),
    ("Base de Australaves: Falconiformes+Psittacopasserae (Eufalconimorphae) | Cariamiformes+Falconiformes",
     "australaves_base", ["Eufalconimorphae (Falconiformes+Psittacopasserae)", "Falconiformes+Cariamiformes"],
     "australaves_base", "Falconidae-Cariamidae"),
    ("Psittacopasserae: Psittaciformes+Passeriformes | Falconiformes+Psittaciformes",
     "psittacopasserae", ["Psittacopasserae", "Falconiformes+Psittaciformes"], "psittacopasserae",
     "Falconidae-Psittaciformes"),
    ("Cathartidae irma de Accipitriformes s.s. (Sagittariidae+Pandionidae+Accipitridae)",
     "accipitrimorphae", ["Accipitriformes s.s. (Sagittariidae+Pandionidae+Accipitridae)"],
     "accipitrimorphae", "Accipitridae-Cathartidae"),
    ("Pandionidae+Accipitridae | Sagittariidae+Accipitridae",
     "accipitriformes_ss", ["Pandionidae+Accipitridae", "Sagittariidae+Accipitridae"], "accipitriformes_ss",
     "Accipitridae-Pandionidae"),
    ("Strigidae monofiletico, Tytonidae irmao (Strigiformes)",
     "strigidae", ["Strigidae", "Tytonidae"], "strigiformes", "Tytonidae-Strigidae"),
    ("Raiz de Telluraves: topologia enraizada de A=Accipitrimorphae, S=Strigiformes, C=Coraciimorphae, U=Australaves",
     None, [], "raiz_telluraves_A_S_C_U", None),
    ("Falconidae: Herpetotheriinae irma de Polyborinae+Falconinae",
     None, ["Herpetotheriinae (Herpetotheres+Micrastur)", "Falconinae+Polyborinae (Falco+Caracara)"],
     "falconidae", "Herpetotheres-Falco"),
]

# Arvores-fonte que nao sao inferencia a partir de dados proprios (pelo titulo do
# estudo ou da arvore na Open Tree): superarvore, arvore compilada para estudo
# comparativo, esqueleto de outra arvore, datacao sobre topologias previas.
DERIVADAS = {
    "ot_809@tree2",      # Jetz 2012, "Hackett backbone"
    "ot_2018@tree8",     # Kimball 2019, superarvore filogenomica
    "ot_2018@tree9",
    "ot_116@tree1",      # Hugall & Stuart-Fox 2012, estudo comparativo
    "ot_116@tree4",
    "ot_1174@tree3",     # D'Emic & Benson 2013, estudo comparativo (osteocitos)
    "ot_534@tree1",      # Ksepka & Phillips 2015, estudo de datacao
    "ot_2879@tree1",     # Claramunt et al. 2025, arvore datada de 9 mil especies (esqueleto)
}

ABREV = {
    "Accipitrimorphae": "Acc", "Strigiformes": "Str", "Coraciimorphae": "Cor", "Cariamiformes": "Car",
    "Falconidae": "Fal", "Falconiformes": "Fal", "Psittacopasserae": "PsP", "Psittaciformes": "Psi",
    "Passeriformes": "Pas", "Cathartidae": "Cat", "Sagittariidae": "Sag", "Pandionidae": "Pan",
    "Accipitridae": "Acd", "Pandionidae+Accipitridae": "Pan+Acd", "Glaucidium": "Gla",
    "Strix+Ciccaba": "Strix", "Tytonidae": "Tyto", "Herpetotheriinae": "Her", "Polyborinae": "Pol",
    "Falconinae": "Fni",
}


def abrev(par: str) -> str:
    # "X+Y" onde X ou Y podem conter "+" (Pandionidae+Accipitridae)
    for k in sorted(ABREV, key=len, reverse=True):
        par = par.replace(k, ABREV[k])
    return par


TAXONOMIAS = ("ott", "ot_2770@", "ot_2019@")  # OTT, Clements 2024, eBird: arvores taxonomicas


def filog(lista: str) -> int:
    """Conta so arvores filogeneticas (tira as taxonomicas da lista do annotations)."""
    return sum(1 for x in lista.split(";") if x and not x.startswith(TAXONOMIAS))


def quart(linha: dict[str, str]) -> str:
    qs = sorted(((k[2:-1], float(v)) for k, v in linha.items() if k.startswith("q[") and v), key=lambda kv: -kv[1])
    return " / ".join(f"{abrev(k)} {v:.3f}" for k, v in qs)


def main() -> None:
    qr = tsv(R / "b10k" / "quartetos_ramos.tsv")
    pub = tsv(R / "b10k" / "arvores_publicadas_ramos.tsv")
    gcf = tsv(R / "b10k" / "gcf_ramos.tsv")
    otl = {x["no"]: x for x in tsv(R / "opentree" / "nos.tsv")}
    aves = {x["no"]: x for x in tsv(R / "opentree" / "aves_1.6" / "nos_aves16.tsv")}
    votos = tsv(R / "opentree" / "votos_por_arvore.tsv")
    tt = {x["par"]: x for x in tsv(R / "timetree" / "timetree_pares.tsv")}
    dat = {x["par"]: x for x in tsv(R / "timetree" / "idades_arvores_datadas.tsv")}

    # escore total de quartetos (ASTRAL-IV -C) das topologias alternativas fixas
    escores: dict[str, dict[str, int]] = {}
    for arq in sorted((R / "b10k" / "topologias").glob("escores_*.tsv")):
        for linha in arq.read_text().splitlines():
            c = linha.split("	")
            if len(c) == 3 and c[2].isdigit():
                escores.setdefault(c[0], {})[c[1]] = int(c[2])
    delta_por_perg = {
        "afroaves_base": ["T4_Afroaves_S+C", "T5_Afroaves_A+C"],
        "australaves_base": ["T2_Australaves_K+F", "T3_Australaves_K+P"],
    }

    linhas = []
    for rot, perg, nos_otl, perg_v, par in NOS:
        li: dict[str, str] = {"no_chave": rot}
        if perg:
            for x in qr:
                if x["pergunta"] != perg or not x.get("EN"):
                    continue
                tipo, desc = x["analise"].split(" | ")
                chave = f"B10K {tipo} {desc.split(',')[0]} {'wASTRAL' if 'wASTRAL' in desc else 'ASTRAL-IV'}"
                chave += " livre" if "livre" in desc else " (pontua 63K)"
                li[chave] = (f"{quart(x)}; localPP={x['localPP_da_resolucao']} ({abrev(x['resolucao_na_arvore'])}); "
                             f"EN={x['EN']}; p_polit={x['p_politomia']}; p_sim={x['p_simetria_menores']}")
            for x in pub:
                if x["pergunta"] == perg:
                    li[f"publicada: {x['arvore']}"] = f"{abrev(x['resolucao'])} [{x['rotulo_do_ramo(suporte ou p)']}]"
            for x in gcf:
                if x["pergunta"] == perg and x.get("gCF"):
                    li[f"gCF {x['tipo']} (ramo 63K {abrev(x['ramo_63K'])})"] = (
                        f"gCF {x['gCF']} gDF1 {x['gDF1']} gDF2 {x['gDF2']} gDFP {x['gDFP']} (gN {x['gN']})")
        for tipo, sc in escores.items():
            base = sc.get("T1_63K")
            alts = delta_por_perg.get(perg or "", [])
            if base and alts:
                li[f"B10K {tipo} 42 taxa: escore(63K) - escore(alternativa), quartetos brutos"] = "; ".join(
                    f"{a.split('_', 1)[1]}: {base - sc[a]:+,}" for a in alts if a in sc)
            if base and perg == "afroaves_base":
                outras = [a for a in ("T7_Acc_irma_de_Telluraves", "T8_AccStr_irma_de_Telluraves") if a in sc]
                if outras:
                    li[f"B10K {tipo} 42 taxa: escore(63K) - escore(sem Afroaves)"] = "; ".join(
                        f"{a.split('_', 1)[1]}: {base - sc[a]:+,}" for a in outras)
        for n in nos_otl:
            o, a = otl.get(n), aves.get(n)
            # no ausente da sintese: as contagens seriam do MRCA (outro no), nao do clado pedido
            if o:
                li[f"OpenTree16.1: {n}"] = (
                    f"sustentam={filog(o['sustentam'])}; conflitam={filog(o['conflitam'])}"
                    if o["existe_na_sintese"] == "sim" else "ausente da sintese")
            if a:
                li[f"Aves1.6: {n}"] = (
                    f"sustentam={filog(a['sustentam'])}; conflitam={filog(a['conflitam'])}; "
                    f"no entra por: {a['fonte_de_maior_prioridade_que_sustenta']}"
                    if a["existe_na_sintese"] == "sim" else "ausente da sintese")
        cont: dict[str, int] = {}
        cont_p: dict[str, int] = {}
        for v in votos:
            q = v.get(perg_v if perg_v.startswith("raiz_") else perg_v + "_quartetos", "")
            if "=" not in q:
                continue
            d = {k: float(x) for k, x in (kv.split("=") for kv in q.split(";"))}
            d.pop("nr", None)
            melhor = max(d, key=lambda k: d[k])
            if d[melhor] >= 0.5:
                cont[melhor] = cont.get(melhor, 0) + 1
                if v["arvore"] not in DERIVADAS:
                    cont_p[melhor] = cont_p.get(melhor, 0) + 1
        li["OTL arvores-fonte, todas (maioria de quartetos >= 0,5)"] = "; ".join(
            f"{abrev(k)}: {n}" for k, n in sorted(cont.items(), key=lambda kv: -kv[1]))
        li["OTL arvores-fonte, so inferencia primaria"] = "; ".join(
            f"{abrev(k)}: {n}" for k, n in sorted(cont_p.items(), key=lambda kv: -kv[1]))
        t = tt.get(par) if par else None
        d2 = dat.get(par) if par else None
        if t:
            li[f"TimeTree {par}"] = (f"mediana {t['mediana_TT_Ma']} Ma (Q1-Q3 {t['Q1_Ma']}-{t['Q3_Ma']}; "
                                     f"amplitude {t['min_Ma']}-{t['max_Ma']}; n={t['n_estimativas']}; "
                                     f"estimada {t['idade_estimada_TT_Ma']} [{t['IC_TT_baixo']}-{t['IC_TT_alto']}])")
        if d2:
            li[f"datadas {par}"] = (f"B10K {d2['B10K_idade_Ma']} [{d2['B10K_HPD95']}]; "
                                    f"Aves1.6 {d2['Aves1.6_idade_Ma']} [{d2['Aves1.6_HPD95']}]")
        linhas.append(li)

    cab: list[str] = []
    for li in linhas:
        for k in li:
            if k not in cab:
                cab.append(k)
    with open(R / "tabela_sintese.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cab) + "\n")
        for li in linhas:
            f.write("\t".join(li.get(c, "") for c in cab) + "\n")
    print(f"tabela_sintese.tsv: {len(linhas)} nos x {len(cab)} colunas")


if __name__ == "__main__":
    main()
