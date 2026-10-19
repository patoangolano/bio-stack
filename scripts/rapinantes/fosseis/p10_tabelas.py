"""Tabelas em Markdown, geradas dos TSV, para o resumo (evita transcrever numero a mao).

Escreve runs/rapinantes/fosseis/tabelas_resumo.md.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict

import clados as cl

B = cl.BASE


def ler(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def md(cab, linhas):
    out = ["| " + " | ".join(cab) + " |", "|" + "---|" * len(cab)]
    out += ["| " + " | ".join(str(x) for x in li) + " |" for li in linhas]
    return "\n".join(out)


def main() -> None:
    partes = []
    # 2a
    g = ler(B / "arvores" / "grupos.tsv")
    por = defaultdict(dict)
    for r in g:
        por[r["grupo"]][r["matriz"]] = f"{r['UFBoot_pct']}{'*' if r['na_arvore_ML'] == '1' else ''}"
    ids = sorted({r["matriz"] for r in g})
    linhas = [
        [k] + [por[k].get(i, "—") for i in ids]
        for k in cl.TRADICIONAIS + cl.MOLECULARES
        if k in por
    ]
    partes.append(
        "### 2a. Grupos (UFBoot %, morfologia sozinha; * = na arvore ML)\n\n"
        + md(["grupo"] + ids, linhas)
    )
    # distancias
    d = ler(B / "arvores" / "distancias.tsv")
    partes.append(
        "### Splits do esqueleto molecular recuperados pela ML morfologica\n\n"
        + md(
            ["matriz", "atuais", "splits esqueleto", "recuperados", "fracao"],
            [
                [
                    r["matriz"],
                    r["n_atuais"],
                    r["splits_esqueleto"],
                    r["recuperados_na_ML_morfologica"],
                    r["fracao"],
                ]
                for r in d
            ],
        )
    )
    # 2b
    pos = ler(B / "arvores" / "fosseis_posicao.tsv")
    chave = {
        "Masillaraptor",
        "Danielsraptor",
        "Messelastur",
        "Tynskya",
        "Messelasturidae",
        "Halcyornithidae",
        "Pseudasturidae",
        "Pseudasturides",
        "Pseudasturides_macrocephalus",
        "Cyrilavis_colburnorum",
        "Bathornis",
        "Elaphrocnemus",
        "Paracrax",
        "Psilopterus",
        "Llallawavis",
        "Dryornis_pampeanus",
        "Pampagyps_imperator",
    }
    top: dict[tuple, list] = defaultdict(list)
    for r in pos:
        if r["fossil"] in chave:
            top[(r["matriz"], r["analise"], r["fossil"])].append(
                f"{r['irmao_UFBoot']} {r['UFBoot_pct']}"
            )
    ml = {(r["matriz"], r["analise"], r["fossil"]): r["irmao_na_ML"] for r in pos}
    partes.append(
        "### 2b. Irmao do fossil (outros fosseis podados): ML e as 3 posicoes mais frequentes no UFBoot\n\n"
        + md(
            ["matriz", "analise", "fossil", "irmao na ML", "UFBoot (top 3, %)"],
            [[k[0], k[1], k[2], ml[k], "; ".join(v)] for k, v in sorted(top.items())],
        )
    )
    pan = ler(B / "arvores" / "fosseis_pan.tsv")
    linhas = [
        [r["matriz"], r["fossil"], r["grupo_total"], r["UFBoot_pct"]]
        for r in pan
        if r["analise"] == "esqueleto"
        and r["fossil"] in chave
        and float(r["UFBoot_pct"]) >= 5
        and r["grupo_total"] not in ("Telluraves",)
    ]
    partes.append(
        "### 2b. Grupo total (fossil dentro ou no tronco de X), analise com esqueleto, UFBoot >= 5%\n\n"
        + md(["matriz", "fossil", "grupo total", "UFBoot %"], linhas)
    )
    # 2c
    au = ler(B / "au" / "au_resultados.tsv")
    linhas = [[r["matriz"], r["fossil"], r["hipotese"], r["deltaL"], r["p-AU"]] for r in au]
    partes.append(
        "### 2c. Teste AU\n\n"
        + md(["matriz", "fossil/teste", "hipotese", "deltaL", "p-AU"], linhas)
    )
    # 2d
    for nome, arq in (
        ("ensemble", "ensemble_CI_RI.tsv"),
        ("permutacao", "permutacao_rapinantes.tsv"),
        ("uniao", "uniao_falcao_gaviao.tsv"),
    ):
        p = B / "homoplasia" / arq
        if p.exists():
            rs = ler(p)
            cab = list(rs[0].keys())
            partes.append(f"### 2d. {nome}\n\n" + md(cab, [[r[c] for c in cab] for r in rs]))
    # 2e
    p = B / "morfoespaco" / "testes_morfoespaco.json"
    if p.exists():
        js = json.loads(p.read_text(encoding="utf-8"))
        linhas = []
        for j in js:
            if "PERMANOVA_clado" not in j:
                continue

            def f(k: str, j: dict = j) -> str:
                return f"R2={j[k]['R2']}, p={j[k]['p']}" if j[k]["R2"] is not None else "—"

            linhas.append(
                [
                    j["matriz"],
                    j["n_atuais_Telluraves"],
                    j["n_rapinantes"],
                    f("PERMANOVA_clado"),
                    f("PERMANOVA_guilda"),
                    f("PERMANOVA_guilda_dado_clado"),
                    f"r={j['Mantel_guilda_dado_clado']['r_parcial']}, p={j['Mantel_guilda_dado_clado']['p']}",
                    f"{j.get('MORD_medio_rapinantes_entre_clados')} x {j.get('MORD_medio_outros_pares_entre_clados')}"
                    f" (p={j.get('p_rapinantes_mais_proximos_que_acaso')})",
                ]
            )
        partes.append(
            "### 2e. Morfoespaco\n\n"
            + md(
                [
                    "matriz",
                    "n Telluraves",
                    "n rapinantes",
                    "PERMANOVA clado",
                    "PERMANOVA guilda",
                    "guilda | clado",
                    "Mantel guilda | clado",
                    "MORD rapinantes entre clados x outros pares",
                ],
                linhas,
            )
        )
    (B / "tabelas_resumo.md").write_text("\n\n".join(partes) + "\n", encoding="utf-8")
    print("ok")


if __name__ == "__main__":
    main()
