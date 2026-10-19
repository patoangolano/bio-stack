"""
Tabelas em Markdown para o resumo (research_notes): quartetos do B10K por ramo
e a tabela-sintese compacta. Le os TSV ja gerados; nao recalcula nada.

    uv run python scripts/rapinantes/reanalise/tabelas_md.py

Grava runs/rapinantes/reanalise/tabelas_md.md.
"""

from __future__ import annotations

import csv
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]
R = RAIZ / "runs" / "rapinantes" / "reanalise"

PERG = {
    "afroaves_base": ("Base de Afroaves", ["Accipitrimorphae+Strigiformes", "Strigiformes+Coraciimorphae",
                                            "Accipitrimorphae+Coraciimorphae"], ["A+S", "S+C", "A+C"]),
    "afroaves_mono_se_AS": ("Monofilia de Afroaves (na topologia 63K)",
                            ["Accipitrimorphae+Strigiformes+Coraciimorphae", "Coraciimorphae+Australaves",
                             "Accipitrimorphae+Strigiformes+Australaves"], ["(A+S)+C", "C+U", "(A+S)+U"]),
    "australaves_base": ("Base de Australaves", ["Falconidae+Psittacopasserae", "Cariamiformes+Falconidae",
                                                  "Cariamiformes+Psittacopasserae"], ["F+PsP", "K+F", "K+PsP"]),
    "psittacopasserae": ("Psittacopasserae", ["Psittaciformes+Passeriformes", "Falconidae+Psittaciformes",
                                              "Falconidae+Passeriformes"], ["Psi+Pas", "F+Psi", "F+Pas"]),
    "accipitrimorphae": ("Cathartidae x Accipitriformes s.s.",
                         ["Sagittariidae+Pandionidae+Accipitridae", "Cathartidae+Sagittariidae",
                          "Cathartidae+Pandionidae+Accipitridae"], ["Sag+(Pan+Acd)", "Cat+Sag", "Cat+(Pan+Acd)"]),
    "accipitriformes_ss": ("Sagittariidae / Pandionidae / Accipitridae",
                           ["Pandionidae+Accipitridae", "Sagittariidae+Pandionidae", "Sagittariidae+Accipitridae"],
                           ["Pan+Acd", "Sag+Pan", "Sag+Acd"]),
    "strigidae": ("Strigidae x Tyto", ["Glaucidium+Strix+Ciccaba", "Glaucidium+Tytonidae",
                                       "Strix+Ciccaba+Tytonidae"], ["Strigidae", "Gla+Tyto", "Strix+Tyto"]),
}
TIPOS = {"intergenico63K": "intergênico 63K", "exon_c12": "éxons C12", "intron": "íntrons", "uce": "UCE"}


def tsv(c: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(open(c, encoding="utf-8"), delimiter="\t")) if c.exists() else []


def dec(s: str) -> str:
    """Decimal com virgula (texto em portugues)."""
    return s.replace(".", ",")


def fmt_p(p: str) -> str:
    try:
        v = float(p)
    except ValueError:
        return p
    if v == 0:
        return "<1e-300"
    return f"{v:.1e}".replace("e-0", "e-") if v < 0.001 else dec(f"{v:.3f}")


def main() -> None:
    qr = tsv(R / "b10k" / "quartetos_ramos.tsv")
    out = []
    for p, (titulo, chaves, curtos) in PERG.items():
        out.append(f"\n**{titulo}** — q de {' / '.join(curtos)}\n")
        out.append("| Dado | Táxons, estimador | EN | " + " | ".join(f"q {c}" for c in curtos)
                   + " | Resolução da árvore (localPP) | p 1ª×2ª | p simetria |")
        out.append("|---|---|---|" + "---|" * len(curtos) + "---|---|---|")
        for tipo, nome in TIPOS.items():
            for x in qr:
                if x["pergunta"] != p or not x.get("EN") or not x["analise"].startswith(tipo + " |"):
                    continue
                desc = x["analise"].split(" | ")[1]
                if "busca livre" in desc and "wASTRAL" not in desc:
                    continue  # mesmos q que a pontuacao quando a topologia coincide; a arvore livre vai em outra tabela
                rot = desc.replace("ASTRAL-IV, pontua topologia 63K", "ASTRAL-IV (63K)").replace(
                    "wASTRAL hibrido (aLRT), busca livre", "wASTRAL (aLRT), livre")
                qs = []
                for k in chaves:
                    v = x.get(f"q[{k}]", "")
                    qs.append(dec(f"{float(v):.3f}") if v else "—")
                res = x["resolucao_na_arvore"]
                try:
                    res_c = curtos[chaves.index(res)]
                except ValueError:
                    res_c = res
                en = f"{float(x['EN']):,.0f}".replace(",", " ")
                pp = dec(f"{float(x['localPP_da_resolucao']):.2f}")
                rot = rot.replace("42 taxa", "42 táxons").replace("363 taxa", "363 táxons")
                out.append(f"| {nome} | {rot} | {en} | " + " | ".join(qs)
                           + f" | {res_c} ({pp}) | {fmt_p(x['p_top2 (1a vs 2a)'])}"
                           f" | {fmt_p(x['p_simetria_menores'])} |")
    (R / "tabelas_md.md").write_text("\n".join(out) + "\n", encoding="utf-8")
    print("\n".join(out))


if __name__ == "__main__":
    main()
