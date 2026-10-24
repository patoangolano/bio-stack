"""A auditoria de `nsig`: amostrar alvos e exato, e as constantes batem o preset.

Roda no python do `envs/cell-eval2`:

    ./envs/cell-eval2/.venv/Scripts/python.exe -I -m pytest tests/vcc/test_auditoria.py

Tres afirmacoes, e as tres sao o que faz a auditoria significar algo:

1. **Amostrar e exato.** `fdr_scope=per_pert` diz que o BH e refeito dentro de
   cada alvo, entao o `nsig` de um alvo nao pode depender de quais outros alvos
   entraram. Se dependesse, uma amostra de 40 alvos nao diria nada sobre os 300
   e o modulo estaria prometendo exatidao que nao tem.
2. **As constantes duplicadas batem o preset instalado.** `P_ADJ`, `GATE_CPM` e
   companhia estao copiadas do `configs/vcc2026.yaml` em vez de lidas, porque a
   chave pode mudar de lugar entre versoes. Copia sem conferencia e divergencia
   silenciosa esperando a hora.
3. **A leitura por fatia acerta as celulas.** Ler 400 de 120.000 linhas de um
   CSR por fatia minima e onde um erro de indice sairia como um `nsig`
   plausivel e errado.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip(
    "cell_eval2",
    reason=(
        "estes testes CHAMAM o cell_eval2 instalado, que e o ponto deles. Rodam com: "
        "./envs/cell-eval2/.venv/Scripts/python.exe -I -m pytest tests/vcc/test_auditoria.py"
    ),
)

import anndata as ad  # noqa: E402
import auditoria as aud  # noqa: E402
import cell_eval2  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import yaml  # noqa: E402

N_GENES = 120
N_CELULAS = 40
N_CTRL = 400
ALVOS = [f"G{i:03d}" for i in range(6)]


def escrever(caminho: Path, rotulos, contextos, matriz) -> Path:
    obs = pd.DataFrame(
        {
            "target_gene": pd.Categorical(rotulos),
            "context": pd.Categorical(contextos),
        },
        index=pd.Index([f"c{i}" for i in range(len(rotulos))], dtype=object),
    )
    var = pd.DataFrame(index=pd.Index([f"G{i:03d}" for i in range(N_GENES)], dtype=object))
    ad.AnnData(X=sp.csr_matrix(matriz.astype(np.float32)), obs=obs, var=var).write_h5ad(caminho)
    return caminho


def perfil(rng, deslocar: slice | None, forca: float = 4.0) -> np.ndarray:
    """Perfil quase plano, com um bloco de genes deslocado — o bastante para dar DE."""
    peso = np.full(N_GENES, 1.0 / N_GENES)
    if deslocar is not None:
        peso = peso.copy()
        peso[deslocar] *= forca
        peso /= peso.sum()
    return rng.multinomial(int(rng.integers(20_000, 30_000)), peso)


@pytest.fixture
def mundo(tmp_path):
    """Submissao de 6 alvos e um controle real separado, no formato do controls-final."""
    rng = np.random.default_rng(3)
    linhas, rotulos, ctxs = [], [], []
    for i, alvo in enumerate(ALVOS):
        fatia = slice(i * 15, i * 15 + 12)
        for _ in range(N_CELULAS):
            linhas.append(perfil(rng, fatia))
            rotulos.append(alvo)
            ctxs.append("D")
    submissao = escrever(tmp_path / "sub.h5ad", rotulos, ctxs, np.asarray(linhas))

    ctrl = np.asarray([perfil(rng, None) for _ in range(N_CTRL)])
    controles = escrever(
        tmp_path / "context_D.h5ad",
        ["non-targeting"] * N_CTRL,
        ["D"] * N_CTRL,
        ctrl,
    )
    return submissao, controles


def test_as_constantes_batem_o_preset_instalado():
    """Constante copiada sem conferencia e divergencia silenciosa esperando a hora."""
    raiz = cell_eval2.__file__
    assert raiz is not None
    caminho = Path(raiz).parent / "configs" / "vcc2026.yaml"
    preset = yaml.safe_load(caminho.read_text(encoding="utf-8"))
    assert preset["de"]["p_adj_threshold"] == aud.P_ADJ
    assert preset["de"]["min_abs_log2fc"] == aud.MIN_ABS_LOG2FC
    assert preset["de"]["sort_by"] == aud.SORT_BY
    assert preset["de"]["method"] == aud.METODO
    assert preset["de"]["nan_lfc_policy"] == aud.NAN_LFC
    assert preset["de"]["mean_calc"] == aud.MEAN_CALC
    assert preset["de"]["fdr_scope"] == aud.FDR_SCOPE
    assert preset["filter"]["filter_gene_min_cpm_cell"] == aud.GATE_CPM
    assert preset["target_sum"] == aud.TARGET_SUM
    assert preset["input_type"] == aud.TIPO_DE_ENTRADA
    assert preset["control"] == aud.CONTROLE


def test_amostrar_alvos_nao_muda_o_nsig_de_nenhum(mundo):
    """A afirmacao que torna a amostra exata, e nao uma estimativa.

    Com `fdr_scope=per_pert` o BH e refeito dentro de cada alvo sobre os
    p-valores dele. Logo auditar 2 alvos tem de dar, NESSES DOIS, o mesmo numero
    que auditar os 6. Se o escopo virasse global, este teste falharia — e seria a
    hora de parar de prometer exatidao.
    """
    submissao, controles = mundo
    _, todos = aud.auditar(submissao, controles, "D", n_alvos=0, semente=0, threads=1)
    _, parcial = aud.auditar(submissao, controles, "D", n_alvos=2, semente=0, threads=1)
    assert len(parcial) == 2
    assert set(parcial) <= set(todos)
    for alvo, valor in parcial.items():
        assert valor == todos[alvo], f"{alvo}: {valor} na amostra, {todos[alvo]} no conjunto"


def test_as_duas_contagens_de_nsig_concordam(mundo):
    """`p_adj < 0,05` agrupado e a largura da coluna na matriz de postos.

    A segunda e a que os membros de DE consomem. Se divergissem, a definicao de
    "significativo" da auditoria nao seria a do placar, e o relatorio diria.
    """
    submissao, controles = mundo
    laudo, _ = aud.auditar(submissao, controles, "D", n_alvos=0, semente=0, threads=1)
    assert laudo.duas_contagens_concordam
    assert laudo.maior_divergencia == 0


def test_o_gate_de_cpm_corta_genes_e_e_o_mesmo_para_todos(mundo):
    """O gate vem SO da referencia (#351), entao e um conjunto unico por contexto."""
    submissao, controles = mundo
    laudo, _ = aud.auditar(submissao, controles, "D", n_alvos=0, semente=0, threads=1)
    assert laudo.genes_no_eixo == N_GENES
    # O perfil do controle e plano em 120 genes com ~25.000 UMI: cada gene fica
    # em ~8.300 CPM, muito acima de 5, entao o gate nao deve cortar nada aqui.
    # O que se afirma e que o numero EXISTE e e consistente, nao que corte.
    assert 0 < laudo.genes_apos_o_gate <= N_GENES


def test_leitura_por_fatia_acerta_as_celulas_pedidas(mundo):
    """Ler 2 alvos de 6 tem de dar exatamente as linhas deles, na ordem do arquivo."""
    submissao, _ = mundo
    pedidos = {ALVOS[1], ALVOS[4]}
    recorte = aud.ler_celulas(submissao, "D", pedidos)
    assert recorte.n_obs == 2 * N_CELULAS
    assert set(recorte.obs["target_gene"].astype(str)) == pedidos
    inteiro = ad.read_h5ad(submissao)
    esperado = inteiro[inteiro.obs["target_gene"].astype(str).isin(sorted(pedidos))]
    assert np.array_equal(
        sp.csr_matrix(recorte.X).toarray(), sp.csr_matrix(esperado.X).toarray()
    )


def test_leitura_por_fatia_acerta_quando_o_grupo_nao_e_contiguo(mundo, tmp_path):
    """Um arquivo de terceiro pode intercalar os grupos; o recorte nao pode errar."""
    submissao, _ = mundo
    inteiro = ad.read_h5ad(submissao)
    rng = np.random.default_rng(9)
    embaralhada = tmp_path / "embaralhada.h5ad"
    inteiro[rng.permutation(inteiro.n_obs)].copy().write_h5ad(embaralhada)
    pedidos = {ALVOS[0], ALVOS[3]}
    recorte = aud.ler_celulas(embaralhada, "D", pedidos)
    assert recorte.n_obs == 2 * N_CELULAS
    assert set(recorte.obs["target_gene"].astype(str)) == pedidos
    # As somas por gene nao dependem da ordem das celulas: compara o conjunto.
    de_la = ad.read_h5ad(embaralhada)
    esperado = de_la[de_la.obs["target_gene"].astype(str).isin(sorted(pedidos))]
    assert np.array_equal(
        np.asarray(sp.csr_matrix(recorte.X).sum(axis=0)).ravel(),
        np.asarray(sp.csr_matrix(esperado.X).sum(axis=0)).ravel(),
    )


def test_semente_fixa_da_a_mesma_amostra(mundo):
    """`default_rng(semente)` e nao `np.random.seed`: a amostra e reproduzivel."""
    submissao, controles = mundo
    _, a = aud.auditar(submissao, controles, "D", n_alvos=3, semente=11, threads=1)
    _, b = aud.auditar(submissao, controles, "D", n_alvos=3, semente=11, threads=1)
    _, c = aud.auditar(submissao, controles, "D", n_alvos=3, semente=12, threads=1)
    assert a == b
    assert set(a) != set(c) or len(ALVOS) <= 3
