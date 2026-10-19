"""
f08 — tabela dos nos-chave por conjunto de dados e sensibilidade.

Para cada execucao do pipeline (rotulo=run_id) e cada no-chave (clado, lido como
biparticao, entao a raiz nao importa):
  * supermatriz: presente?, SH-aLRT / aBayes / UFBoot, gCF, gDF1, gDF2, gDFP (e N),
    sCF, sDF1, sDF2 (sCFL); qui-quadrado gDF1_N vs gDF2_N e sDF1_N vs sDF2_N
    (H0: as duas alternativas NNI na mesma frequencia, como sob ILS pura;
    Minh, Hahn & Lanfear 2020; Lanfear, blog "concordance factors" 2018);
  * wASTRAL: presente?, probabilidade local (localPP);
  * arvores de genes: quantas tem o clado (restrito aos taxons do gene) com
    SH-aLRT >= 80 e UFBoot >= 95, e quantas tem um ramo nesse criterio que o
    contradiz (biparticao incompativel).
Os ramos do .cf.stat sao casados com a arvore .cf.tree pelo rotulo e pelo
comprimento (os IDs do gCF e do sCFL nao coincidem entre si).

    uv run python scripts/rapinantes/filogenomica/f08_nos_chave.py <saida.tsv> rotulo=run_id [...]
"""

from __future__ import annotations

import csv
import importlib.util
import io
import json
import re
import sys
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from scipy.stats import chisquare

RAIZ = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location("f07", Path(__file__).with_name("f07_hipoteses.py"))
assert _spec and _spec.loader
f07 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(f07)
G: dict[str, set[str]] = dict(f07.G)
G["TELL"] = G["AFRO"] | G["AUST"]
G["STRISORES"] = {"Calypte_anna", "Apus_apus", "Antrostomus_carolinensis"}
G["STRIGIDAE"] = G["STRIG"] - {"Tyto_alba"}
G["ACCIPITRIDAE_SEM_ELANUS"] = G["ACCIPITRIDAE"] - {"Elanus_axillaris"}
G["FALC_SEM_HERPETOTHERES"] = G["FALC"] - {"Herpetotheres_cachinnans"}

NOS: list[tuple[str, str]] = [
    ("Telluraves", "TELL"),
    ("Afroaves", "AFRO"),
    ("Australaves", "AUST"),
    ("Accipitrimorphae", "ACCMORPH"),
    ("Accipitriformes s.s. (sem Cathartidae)", "ACC"),
    ("Cathartidae", "CATH"),
    ("Pandion + Accipitridae", "PAN+ACCIPITRIDAE"),
    ("Sagittarius + Accipitridae", "SAG+ACCIPITRIDAE"),
    ("Sagittarius + Pandion", "SAG+PAN"),
    ("Accipitridae", "ACCIPITRIDAE"),
    ("Accipitridae exceto Elanus", "ACCIPITRIDAE_SEM_ELANUS"),
    ("Strigiformes", "STRIG"),
    ("Strigidae", "STRIGIDAE"),
    ("Coraciimorphae", "CORAC"),
    ("Strigiformes + Coraciimorphae (H1/H3)", "STRIG+CORAC"),
    ("Strigiformes + Accipitrimorphae (H2)", "STRIG+ACCMORPH"),
    ("Accipitrimorphae + Coraciimorphae (H1b)", "ACCMORPH+CORAC"),
    ("Eufalconimorphae", "EUFAL"),
    ("Psittacopasserae", "PSPA"),
    ("Falconiformes", "FALC"),
    ("Falconidae exceto Herpetotheres", "FALC_SEM_HERPETOTHERES"),
    ("Cariamiformes + Falconiformes (H7b)", "CARI+FALC"),
    ("Falconidae + Accipitrimorphae (H4)", "FALC+ACCMORPH"),
    ("Strigiformes + Caprimulgidae (H5)", "STRIG+CAPRI"),
    ("Cathartidae + Ciconiiformes (H6)", "CATH+CICO"),
    ("Strisores", "STRISORES"),
]


def conj(expr: str, taxa: set[str]) -> frozenset[str]:
    s: set[str] = set()
    for p in expr.split("+"):
        s |= G[p]
    return frozenset(s & taxa)


def ler(caminho: Path) -> Tree:
    txt = re.sub(r":nan", "", caminho.read_text(encoding="utf-8").strip())
    # rotulos internos viram nome entre aspas, para o parser nao tentar le-los como numero
    txt = re.sub(r"\)([^:,;()']+)", lambda m: ")'" + m[1] + "'", txt)
    return Phylo.read(io.StringIO(txt), "newick")  # pyright: ignore[reportPrivateImportUsage]


def biparticoes(arv: Tree) -> tuple[set[str], list[tuple[frozenset[str], str, float]]]:
    """(taxa, [(lado sem a 1a folha em ordem alfabetica, rotulo, comprimento)])."""
    taxa = {t.name for t in arv.get_terminals()}
    ref = min(taxa)
    out = []
    for c in arv.get_nonterminals():
        f = frozenset(t.name for t in c.get_terminals())
        if len(f) < 2 or len(taxa) - len(f) < 2:
            continue
        lado = frozenset(taxa - f) if ref in f else f
        out.append((lado, c.name or (str(c.confidence) if c.confidence is not None else ""), c.branch_length or 0.0))
    return taxa, out


def canon(s: frozenset[str], taxa: set[str]) -> frozenset[str]:
    return frozenset(taxa - s) if min(taxa) in s else s


def compativel(a: frozenset[str], b: frozenset[str], taxa: set[str]) -> bool:
    a2, b2 = frozenset(taxa - a), frozenset(taxa - b)
    return not (a & b and a & b2 and a2 & b and a2 & b2)


def ler_stat(caminho: Path) -> list[dict[str, str]]:
    linhas = [ln for ln in caminho.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    return list(csv.DictReader(linhas, delimiter="\t"))


def casar(stat: list[dict[str, str]], biparts: list[tuple[frozenset[str], str, float]]) -> dict[frozenset[str], dict[str, str]]:
    m: dict[frozenset[str], dict[str, str]] = {}
    for row in stat:
        comp = float(row["Length"])
        cands = [b for b in biparts if abs(b[2] - comp) <= 1e-6 * max(1.0, comp) + 5e-7 and b[1].startswith(row["Label"])]
        if len(cands) != 1:
            cands = [b for b in biparts if abs(b[2] - comp) <= 5e-7 + 1e-6 * comp]
        if len(cands) == 1:
            m[cands[0][0]] = row
    return m


def pval(a: str | None, b: str | None) -> str:
    if a is None or b is None:
        return ""
    x, y = float(a), float(b)
    if x + y == 0:
        return ""
    return f"{chisquare([x, y]).pvalue:.3g}"


def main() -> None:
    saida = Path(sys.argv[1])
    linhas = []
    for par in sys.argv[2:]:
        rotulo, run_id = par.split("=", 1)
        run = RAIZ / "runs" / "filogenia" / run_id
        rot_sup = json.loads((run / "resumo.json").read_text(encoding="utf-8")).get("rotulos_suporte") or ["SH-aLRT", "aBayes", "UFBoot"]
        # o resumo.json grava alrt/abayes/ufboot; aqui os nomes de exibicao
        rot_sup = [{"alrt": "SH-aLRT", "abayes": "aBayes", "ufboot": "UFBoot"}.get(x, x) for x in rot_sup]
        taxa, b_concat = biparticoes(ler(run / "05_concatenado" / "concatenado.treefile"))
        concat = {b[0]: b[1] for b in b_concat}
        g_taxa, b_g = biparticoes(ler(run / "06_concordancia" / "concordancia.gcf.cf.tree"))
        s_taxa, b_s = biparticoes(ler(run / "06_concordancia" / "concordancia.scfl.cf.tree"))
        gstat = casar(ler_stat(run / "06_concordancia" / "concordancia.gcf.cf.stat"), b_g)
        sstat = casar(ler_stat(run / "06_concordancia" / "concordancia.scfl.cf.stat"), b_s)
        astral_f = run / "07_arvore_especies" / "especies.tre"
        a_taxa, b_a = biparticoes(ler(astral_f)) if astral_f.is_file() else (set(), [])
        astral = {b[0]: b[1] for b in b_a}
        genes = []
        for gf in sorted((run / "04_arvores_genes").glob("*.treefile")):
            gt, gb = biparticoes(ler(gf))
            genes.append((gt, gb))
        for nome, expr in NOS:
            c = conj(expr, taxa)
            if len(c) < 2 or len(taxa) - len(c) < 2:
                continue
            k = canon(c, taxa)
            r: dict[str, object] = {"conjunto": rotulo, "run_id": run_id, "no": nome, "n_taxons_no": len(c)}
            r["na_supermatriz"] = "sim" if k in concat else "nao"
            if k in concat:
                vals = concat[k].split("/")
                for i, nm in enumerate(rot_sup):
                    r[nm] = vals[i] if i < len(vals) else ""
            g = gstat.get(k)
            if g:
                for col in ("gCF", "gDF1", "gDF2", "gDFP", "gN"):
                    r[col] = g[col]
                r["gCF_N/gDF1_N/gDF2_N"] = f"{g['gCF_N']}/{g['gDF1_N']}/{g['gDF2_N']}"
                r["p_gDF1_vs_gDF2"] = pval(g.get("gDF1_N"), g.get("gDF2_N"))
            s = sstat.get(k)
            if s:
                for col in ("sCF", "sDF1", "sDF2"):
                    r[col] = s[col]
                r["p_sDF1_vs_sDF2"] = pval(s.get("sDF1_N"), s.get("sDF2_N"))
            if a_taxa:
                ka = canon(frozenset(c & a_taxa), a_taxa)
                r["no_wASTRAL"] = "sim" if ka in astral else "nao"
                r["localPP"] = astral.get(ka, "")
            a_favor = contra = informativos = 0
            for gt, gb in genes:
                cg = frozenset(c & gt)
                if len(cg) < 2 or len(gt) - len(cg) < 2:
                    continue
                informativos += 1
                kg = canon(cg, gt)
                fortes = []
                for lado, lab, _ in gb:
                    v = lab.split("/")
                    try:
                        sh, ufb = float(v[rot_sup.index("SH-aLRT")]), float(v[rot_sup.index("UFBoot")])
                    except (ValueError, IndexError):
                        continue
                    if sh >= 80 and ufb >= 95:
                        fortes.append(lado)
                if kg in fortes:
                    a_favor += 1
                elif any(not compativel(kg, f, gt) for f in fortes):
                    contra += 1
            r["genes_informativos"] = informativos
            r["genes_fortes_a_favor"] = a_favor
            r["genes_fortes_contra"] = contra
            linhas.append(r)
    campos: list[str] = []
    for r in linhas:
        for k in r:
            if k not in campos:
                campos.append(k)
    saida.parent.mkdir(parents=True, exist_ok=True)
    with saida.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(linhas)
    print(f"{len(linhas)} linhas -> {saida}")


if __name__ == "__main__":
    main()
