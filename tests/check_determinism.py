# /// script
# requires-python = ">=3.11"
# dependencies = ["scanpy==1.11.5", "numpy==2.4.6", "igraph>=0.11", "leidenalg>=0.10"]
# ///
"""
Verifica que semente fixada realmente produz o mesmo resultado nesta instalacao.

Isto nao e paranoia. UMAP e Leiden sao estocasticos, e "fixei a semente" e uma
afirmacao que quase ninguem testa — a pessoa roda uma vez, ve um numero, e
assume. Quando o resultado muda entre execucoes, muda silenciosamente: nenhum
erro, so uma figura diferente e uma contagem de cluster diferente.

Roda sobre matriz sintetica de proposito. Determinismo e propriedade do codigo,
nao do dado; usar dado real tornaria o teste lento e dependente de download.

    uv run --script tests/test_determinism.py
"""

from __future__ import annotations

import sys
import warnings

import numpy as np
import scanpy as sc

SEMENTE = 0
N_CELULAS, N_GENES = 300, 200

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0


def matriz():
    """
    Contagens sinteticas com estrutura de grupo, para que o Leiden tenha o que
    particionar. Sem estrutura, todas as celulas caem num cluster so e o teste
    passaria sem exercitar nada.

    Gerador proprio e explicito: `pytest-randomly` reseeda apenas o RandomState
    legado do numpy, nunca um Generator. Confiar nele aqui seria testar o
    plugin, nao o pipeline.
    """
    rng = np.random.default_rng(12345)
    base = rng.negative_binomial(5, 0.3, size=(N_CELULAS, N_GENES)).astype("float32")
    base[: N_CELULAS // 2, : N_GENES // 4] *= 4  # um subgrupo com perfil distinto
    ad = sc.AnnData(base)
    ad.var_names = [f"gene{i}" for i in range(N_GENES)]
    ad.obs_names = [f"cel{i}" for i in range(N_CELULAS)]
    return ad


def pipeline():
    """Normaliza, PCA, vizinhos, UMAP e Leiden — tudo com semente declarada."""
    ad = matriz()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.pca(ad, n_comps=20, random_state=SEMENTE)
    sc.pp.neighbors(ad, n_neighbors=15, random_state=SEMENTE)
    sc.tl.umap(ad, random_state=SEMENTE)
    # flavor, n_iterations e directed fixados junto com a semente: o padrao do
    # scanpy para o backend do Leiden ja mudou entre versoes, e particoes
    # diferem entre backends. Semente sozinha nao protege disso.
    sc.tl.leiden(
        ad, resolution=1.0, random_state=SEMENTE,
        flavor="igraph", n_iterations=2, directed=False,
    )
    return ad


def main() -> int:
    falhas: list[str] = []

    def check(descricao: str, condicao: bool, detalhe: str = "") -> None:
        print(f"{'ok   ' if condicao else 'FALHA'} {descricao}{'  ' + detalhe if detalhe else ''}")
        if not condicao:
            falhas.append(descricao)

    a, b = pipeline(), pipeline()

    check("PCA reproduz bit a bit", np.array_equal(a.obsm["X_pca"], b.obsm["X_pca"]))
    check("UMAP reproduz bit a bit", np.array_equal(a.obsm["X_umap"], b.obsm["X_umap"]))
    check(
        "Leiden reproduz os mesmos rotulos",
        (a.obs["leiden"].to_numpy() == b.obs["leiden"].to_numpy()).all(),
        f"({a.obs['leiden'].nunique()} clusters)",
    )

    # Semente diferente PRECISA mudar o resultado. Sem esta checagem, um
    # pipeline que ignorasse random_state por completo passaria nas tres acima.
    global SEMENTE
    SEMENTE = 7
    c = pipeline()
    SEMENTE = 0
    check(
        "semente diferente muda o embedding (prova que random_state e lido)",
        not np.array_equal(a.obsm["X_umap"], c.obsm["X_umap"]),
    )

    # O scanpy tipa random_state como _LegacyRandom. Passar um Generator moderno
    # levanta ValueError — e por isso que a regra da stack e: `random_state=<int>`
    # em chamada de biblioteca, `default_rng` so em codigo proprio.
    ad = matriz()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    try:
        # O type: ignore aqui e o proprio ponto do teste. O basedpyright ja
        # recusa esta linha estaticamente ("Generator nao e _LegacyRandom") —
        # o que confirma a regra sem nem rodar. A checagem em runtime existe
        # para o caso de o scanpy passar a aceitar e silenciosamente ignorar.
        sc.pp.pca(ad, n_comps=20, random_state=np.random.default_rng(0))  # type: ignore[arg-type]
        recusou = False
    except (ValueError, TypeError):
        recusou = True
    check("scanpy recusa np.random.default_rng em random_state", recusou)

    total = 5
    print(f"\n{total - len(falhas)}/{total} passaram")
    return 1 if falhas else 0


if __name__ == "__main__":
    sys.exit(main())
