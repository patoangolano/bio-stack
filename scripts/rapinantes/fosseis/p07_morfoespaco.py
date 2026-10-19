"""Passo 8 (2e): morfoespaco discreto (MORD e Gower), PCoA e testes de
agrupamento por clado molecular x por guilda rapinante.

Distancia por par de taxa, so nos caracteres codificados nos dois:
  MORD (Lloyd 2016) = soma das diferencas / soma das maximas possiveis
  Gower             = media das diferencas escaladas a [0, 1]
Diferenca por carater: 0 se os conjuntos de estados (polimorfismo incluso) se
tocam; senao 1 (nao ordenado) ou a menor |a - b| em postos (ordenado pelos
autores, escala k - 1). Pares com < 15 caracteres comparaveis sao faltantes;
taxa sao removidos um a um (o de mais faltantes) ate a matriz ficar completa.

Testes so com ATUAIS de Telluraves (fosseis entram na PCoA, nao nos testes):
  PERMANOVA (Anderson 2001) com clado molecular (Accipitriformes, Strigiformes,
  Falconiformes, Cariamiformes, Psittaciformes, Passeriformes, Coraciimorphae)
  e com guilda (rapinante x nao); termo sequencial guilda | clado; Mantel
  parcial distancia ~ mesma-guilda | mesmo-clado e o inverso. 9999 permutacoes,
  rng = np.random.default_rng(20261017).
"""

from __future__ import annotations

import csv
import json

import matplotlib

matplotlib.use("Agg")
import clados as cl
import matplotlib.pyplot as plt
import morfo_io as mi
import numpy as np

SAIDA = cl.BASE / "morfoespaco"
NPERM = 9999
MIN_COMP = 15
CORES = {
    "Accipitriformes": "#b2182b",
    "Strigiformes": "#762a83",
    "Falconiformes": "#e08214",
    "Cariamiformes": "#8c510a",
    "Psittaciformes": "#1b7837",
    "Passeriformes": "#2166ac",
    "Coraciimorphae": "#5aae61",
    "outros": "#999999",
    "fossil": "#000000",
}


def distancias(m: mi.Matriz) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n, c = len(m.taxa), m.nchar
    postos = []
    for j in range(c):
        est = sorted(m.estados_observados(j))
        postos.append({e: i for i, e in enumerate(est)})
    num = np.zeros((n, n))
    den = np.zeros((n, n))
    gow = np.zeros((n, n))
    ncomp = np.zeros((n, n), dtype=int)
    for j in range(c):
        k = len(postos[j])
        if k < 2:
            continue
        ordd = j in m.ordenados
        esc = (k - 1) if ordd else 1
        cel = [m.linhas[t][j] for t in m.taxa]
        for a in range(n):
            ca = cel[a]
            if ca is None:
                continue
            for b in range(a + 1, n):
                cb = cel[b]
                if cb is None:
                    continue
                if ca & cb:
                    d = 0.0
                elif ordd:
                    d = float(min(abs(postos[j][x] - postos[j][y]) for x in ca for y in cb))
                else:
                    d = 1.0
                num[a, b] += d
                den[a, b] += esc
                gow[a, b] += d / esc
                ncomp[a, b] += 1
    num, den, gow, ncomp = (x + x.T for x in (num, den, gow, ncomp))
    with np.errstate(invalid="ignore", divide="ignore"):
        mord = np.where(ncomp >= MIN_COMP, num / den, np.nan)
        gower = np.where(ncomp >= MIN_COMP, gow / ncomp, np.nan)
    np.fill_diagonal(mord, 0.0)
    np.fill_diagonal(gower, 0.0)
    return mord, gower, ncomp


def aparar(d: np.ndarray, nomes: list[str]) -> list[int]:
    idx = list(range(len(nomes)))
    while True:
        sub = d[np.ix_(idx, idx)]
        falt = np.isnan(sub).sum(axis=1)
        if falt.max() == 0:
            return idx
        idx.pop(int(falt.argmax()))


def pcoa(d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = d.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d**2) @ j
    val, vec = np.linalg.eigh(b)
    ordem = np.argsort(val)[::-1]
    val, vec = val[ordem], vec[:, ordem]
    pos = val > 1e-10
    coords = vec[:, pos] * np.sqrt(val[pos])
    return coords, val


def gower_centrada(d: np.ndarray) -> np.ndarray:
    n = d.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * j @ (d**2) @ j


def ss_residual(g: np.ndarray, x: np.ndarray) -> float:
    h = x @ np.linalg.pinv(x)
    r = np.eye(len(g)) - h
    return float(np.trace(r @ g @ r))


def dummies(rot: list[str]) -> np.ndarray:
    nv = sorted(set(rot))
    x = np.ones((len(rot), 1))
    if len(nv) > 1:
        x = np.column_stack([x] + [[1.0 if r == v else 0.0 for r in rot] for v in nv[1:]])
    return x


def permanova(
    d: np.ndarray, rot: list[str], cond: list[str] | None, rng: np.random.Generator
) -> dict:
    """F de Anderson; com cond, termo sequencial rot | cond (permutando linhas)."""
    g = gower_centrada(d)
    n = len(rot)
    sst = float(np.trace(g))

    def f_stat(rr: list[str], cc: list[str] | None) -> tuple[float, float]:
        x0 = dummies(cc) if cc is not None else np.ones((n, 1))
        x1 = np.column_stack([x0, dummies(rr)[:, 1:]])
        ss0, ss1 = ss_residual(g, x0), ss_residual(g, x1)
        gl_t = np.linalg.matrix_rank(x1) - np.linalg.matrix_rank(x0)
        gl_r = n - np.linalg.matrix_rank(x1)
        if gl_t <= 0 or gl_r <= 0:
            return float("nan"), float("nan")
        return ((ss0 - ss1) / gl_t) / (ss1 / gl_r), (ss0 - ss1) / sst

    f_obs, r2 = f_stat(rot, cond)
    if np.isnan(f_obs):
        return {"F": None, "R2": None, "p": None}
    cont = 0
    for _ in range(NPERM):
        p = rng.permutation(n)
        rr = [rot[i] for i in p]
        cc = [cond[i] for i in p] if cond is not None else None
        f, _ = f_stat(rr, cc)
        cont += f >= f_obs
    return {"F": round(f_obs, 3), "R2": round(r2, 3), "p": (cont + 1) / (NPERM + 1)}


def mantel_parcial(d: np.ndarray, a: np.ndarray, b: np.ndarray, rng: np.random.Generator) -> dict:
    """r parcial entre d e a controlando b; permuta linhas/colunas de d."""
    iu = np.triu_indices(len(d), 1)

    def rparc(dd: np.ndarray) -> float:
        x, y, z = dd[iu], a[iu], b[iu]
        rxy, rxz, ryz = np.corrcoef(x, y)[0, 1], np.corrcoef(x, z)[0, 1], np.corrcoef(y, z)[0, 1]
        return (rxy - rxz * ryz) / np.sqrt((1 - rxz**2) * (1 - ryz**2))

    obs = rparc(d)
    cont = 0
    for _ in range(NPERM):
        p = rng.permutation(len(d))
        cont += rparc(d[np.ix_(p, p)]) >= obs
    return {"r_parcial": round(float(obs), 3), "p": (cont + 1) / (NPERM + 1)}


def analisar(id_: str, m: mi.Matriz) -> dict:
    tx = {r["terminal"]: r for r in cl.ler_taxa(id_)}
    mord, gower, ncomp = distancias(m)
    idx = aparar(mord, m.taxa)
    nomes = [m.taxa[i] for i in idx]
    d = mord[np.ix_(idx, idx)]
    dg = gower[np.ix_(idx, idx)]
    fora = sorted(set(m.taxa) - set(nomes))
    coords, val = pcoa(d)
    pos = val[val > 1e-10]
    expl = pos / pos.sum()
    # --- testes nos atuais de Telluraves ---
    tell = [
        i
        for i, t in enumerate(nomes)
        if tx[t]["tipo"] == "atual" and tx[t]["clado"] not in ("outros",)
    ]
    res: dict = {
        "matriz": id_,
        "fonte": m.fonte,
        "n_taxa_pcoa": len(nomes),
        "removidos_por_sobreposicao": fora,
        "pcoa_var_eixo1_2": [round(float(expl[0]), 3), round(float(expl[1]), 3)],
        "autovalores_negativos": int((val < -1e-10).sum()),
    }
    if len(tell) >= 8:
        dt = d[np.ix_(tell, tell)]
        clado = [tx[nomes[i]]["clado"] for i in tell]
        guil = [tx[nomes[i]]["guilda"] for i in tell]
        guil2 = [tx[nomes[i]]["guilda_bico_gancho"] for i in tell]
        rng = np.random.default_rng(20261017)
        res["n_atuais_Telluraves"] = len(tell)
        res["n_rapinantes"] = guil.count("rapinante")
        res["PERMANOVA_clado"] = permanova(dt, clado, None, rng)
        res["PERMANOVA_guilda"] = permanova(dt, guil, None, rng)
        res["PERMANOVA_guilda_dado_clado"] = permanova(dt, guil, clado, rng)
        res["PERMANOVA_clado_dado_guilda"] = permanova(dt, clado, guil, rng)
        res["PERMANOVA_bico_gancho"] = permanova(dt, guil2, None, rng)
        same_c = np.array([[float(a == b) for b in clado] for a in clado])
        same_g = np.array([[float(a == b) for b in guil] for a in guil])
        # distancia ~ -similaridade: usamos 1 - mesma para correlacao positiva com distancia
        res["Mantel_guilda_dado_clado"] = mantel_parcial(dt, 1 - same_g, 1 - same_c, rng)
        res["Mantel_clado_dado_guilda"] = mantel_parcial(dt, 1 - same_c, 1 - same_g, rng)
        # distancia media entre rapinantes de clados diferentes x entre nao-rapinantes de clados diferentes
        rap = [k for k, g in enumerate(guil) if g == "rapinante"]
        par_rr = [
            dt[a, b] for ai, a in enumerate(rap) for b in rap[ai + 1 :] if clado[a] != clado[b]
        ]
        par_cc = [
            dt[a, b]
            for a in range(len(tell))
            for b in range(a + 1, len(tell))
            if clado[a] != clado[b] and not (guil[a] == guil[b] == "rapinante")
        ]
        if par_rr:
            res["MORD_medio_rapinantes_entre_clados"] = round(float(np.mean(par_rr)), 3)
            res["MORD_medio_outros_pares_entre_clados"] = round(float(np.mean(par_cc)), 3)
            # nulo: o mesmo numero de "rapinantes" sorteado entre os atuais de
            # Telluraves; estatistica = distancia media entre sorteados de clados
            # diferentes. p unilateral (rapinantes mais proximos que o acaso).
            obs = float(np.mean(par_rr))
            nulos = []
            for _ in range(NPERM):
                sor = rng.choice(len(tell), size=len(rap), replace=False)
                pp = [
                    dt[a, b]
                    for ai, a in enumerate(sor)
                    for b in sor[ai + 1 :]
                    if clado[a] != clado[b]
                ]
                if pp:
                    nulos.append(np.mean(pp))
            res["p_rapinantes_mais_proximos_que_acaso"] = (1 + sum(x <= obs for x in nulos)) / (
                1 + len(nulos)
            )
    # --- tabelas e figura ---
    with open(SAIDA / f"{id_}_pcoa.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["terminal", "tipo", "clado", "guilda", "PCo1", "PCo2", "PCo3"])
        for i, t in enumerate(nomes):
            w.writerow(
                [
                    t,
                    tx[t]["tipo"],
                    tx[t]["clado"],
                    tx[t]["guilda"],
                    *[f"{coords[i, k]:.4f}" for k in range(min(3, coords.shape[1]))],
                ]
            )
    np.savetxt(
        SAIDA / f"{id_}_MORD.tsv",
        d,
        delimiter="\t",
        fmt="%.4f",
        header="\t".join(nomes),
        comments="",
    )
    np.savetxt(
        SAIDA / f"{id_}_Gower.tsv",
        dg,
        delimiter="\t",
        fmt="%.4f",
        header="\t".join(nomes),
        comments="",
    )
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for i, t in enumerate(nomes):
        r = tx[t]
        cor = CORES.get(r["clado"], "#999999")
        mk = "*" if r["tipo"] == "fossil" else ("^" if r["guilda"] == "rapinante" else "o")
        ax.scatter(
            coords[i, 0],
            coords[i, 1],
            c=cor,
            marker=mk,
            s=90 if mk != "o" else 40,
            edgecolors="k" if r["tipo"] == "fossil" else "none",
            linewidths=0.5,
            zorder=3,
        )
        if r["guilda"] == "rapinante" or r["tipo"] == "fossil" or r["clado"] == "Cariamiformes":
            ax.annotate(
                t.replace("_", " ")[:22],
                (coords[i, 0], coords[i, 1]),
                fontsize=6,
                xytext=(3, 3),
                textcoords="offset points",
            )
    for k, v in CORES.items():
        ax.scatter([], [], c=v, marker="o", label=k)
    ax.scatter([], [], c="w", edgecolors="k", marker="^", label="rapinante (atual)")
    ax.scatter([], [], c="k", marker="*", label="fossil")
    ax.legend(fontsize=6, loc="best", frameon=False)
    ax.set_xlabel(f"PCo1 ({100 * expl[0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({100 * expl[1]:.1f}%)")
    ax.set_title(
        f"{id_}: PCoA de MORD, n = {len(nomes)} taxa ({len(fora)} fora por sobreposicao)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(cl.BASE / "figuras" / f"morfoespaco_{id_}.png", dpi=200)
    plt.close(fig)
    return res


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    (cl.BASE / "figuras").mkdir(exist_ok=True)
    ms = mi.todas()
    todos = []
    for id_ in ["K19", "M11", "MC03", "C25", "D21", "VZ26", "M05", "LZ07"]:
        r = analisar(id_, ms[id_])
        todos.append(r)
        print(json.dumps(r, ensure_ascii=False))
    (SAIDA / "testes_morfoespaco.json").write_text(
        json.dumps(todos, indent=1, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
