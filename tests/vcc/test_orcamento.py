"""A auditoria do orcamento de `expr_mse`: a estrategia de memoria e exata, e a razao e um limite.

Roda no python do `envs/cell-eval2`, porque o modulo CHAMA o `cell_eval2` em vez
de reimplementa-lo e e justamente isso que vale testar:

    ./envs/cell-eval2/.venv/Scripts/python.exe -I -m pytest tests/vcc/test_orcamento.py

Duas afirmacoes load-bearing, e as duas sao sobre o pacote e nao sobre o nosso
codigo:

1. `jackknife_correction` e decomponivel por grupo. A submissao tem 682 milhoes
   de nao nulos por contexto e o modulo varre um grupo por vez; se o jackknife
   de um grupo visse os outros, a varredura daria numero errado sem falhar.
2. Limitar `C_p` por `min(C_p, K * C_p_real)` so pode REDUZIR o pedido, e por
   isso a razao sem limite e um limite inferior da do avaliador. E o unico
   motivo pelo qual uma auditoria sem o lado real de D/E/F significa algo.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip(
    "cell_eval2",
    reason=(
        "estes testes CHAMAM o cell_eval2 instalado, que e o ponto deles. Rodam com: "
        "./envs/cell-eval2/.venv/Scripts/python.exe -I -m pytest tests/vcc/test_orcamento.py"
    ),
)

import anndata as ad  # noqa: E402
import orcamento as orc  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from cell_eval2 import moments  # noqa: E402
from cell_eval2.metrics import delta  # noqa: E402

N_GENES = 180
N_CELULAS = 30
ALVOS = [f"G{i:03d}" for i in range(6)]


def mundo(caminho, *, n_celulas: int = N_CELULAS, contextos=("D",)):
    """Uma submissao pequena com a FORMA da real: contagens inteiras, CSR, obs com contexto."""
    rng = np.random.default_rng(13)
    peso = 1.0 / (1.0 + np.arange(N_GENES)) ** 1.1
    peso /= peso.sum()
    blocos, rotulos, ctxs = [], [], []
    for contexto in contextos:
        for alvo in ALVOS:
            # Cada alvo desloca um bloco de genes proprio: sem dispersao ENTRE
            # alvos o orcamento seria zero e a razao degeneraria.
            proprio = peso.copy()
            inicio = (ALVOS.index(alvo) * 17) % (N_GENES - 20)
            proprio[inicio : inicio + 20] *= 3.0
            proprio /= proprio.sum()
            for _ in range(n_celulas):
                blocos.append(rng.multinomial(int(rng.integers(8_000, 25_000)), proprio))
                rotulos.append(alvo)
                ctxs.append(contexto)
    x = sp.csr_matrix(np.asarray(blocos, dtype=np.float32))
    obs = pd.DataFrame(
        {
            "target_gene": pd.Categorical(rotulos),
            "context": pd.Categorical(ctxs),
        },
        index=pd.Index([f"c{i}" for i in range(len(rotulos))], dtype=object),
    )
    var = pd.DataFrame(index=pd.Index([f"G{i:03d}" for i in range(N_GENES)], dtype=object))
    ad.AnnData(X=x, obs=obs, var=var).write_h5ad(caminho)
    return caminho


@pytest.fixture
def submissao(tmp_path):
    return mundo(tmp_path / "sub.h5ad")


def test_jackknife_por_grupo_e_identico_ao_conjunto():
    """A propriedade de que a varredura por grupo depende, afirmada sobre o pacote.

    Se um dia o `cell_eval2` fizesse o jackknife olhar o conjunto todo, este
    teste falharia e a estrategia de memoria do `orcamento.py` precisaria mudar.
    Sem ele a divergencia sairia como um `C_p` plausivel e errado.
    """
    rng = np.random.default_rng(0)
    n_grupos, n_celulas, n_genes = 5, 40, 250
    peso = 1.0 / (1.0 + np.arange(n_genes)) ** 1.1
    peso /= peso.sum()
    blocos, codigos = [], []
    for p in range(n_grupos):
        for _ in range(n_celulas):
            blocos.append(rng.multinomial(int(rng.integers(8_000, 25_000)), peso))
        codigos += [p] * n_celulas
    x = sp.csr_matrix(np.asarray(blocos, dtype=np.float32))
    codigos = np.asarray(codigos)

    juntos = np.asarray(moments.jackknife_correction(x, codigos, n_grupos, orc.ALVO_SOMA))
    um_a_um = np.array(
        [
            float(
                np.asarray(
                    moments.jackknife_correction(
                        sp.csr_matrix(x[codigos == p]),
                        np.zeros(int((codigos == p).sum()), dtype=np.int64),
                        1,
                        orc.ALVO_SOMA,
                    )
                )[0]
            )
            for p in range(n_grupos)
        ]
    )
    assert np.array_equal(juntos, um_a_um)


def test_limitar_reduz_o_pedido_e_aumenta_a_razao():
    """Por que a razao sem limite e um LIMITE INFERIOR da do avaliador.

    O avaliador faz `pred_term = np.minimum(pred_term, K * real_tn)`. `minimum`
    nunca sobe, `row_w` e positivo, logo o pedido limitado <= pedido sem limite,
    logo razao_limitada >= razao_sem_limite. Afirmado contra `PRED_TRACE_CAP_K`
    do pacote, que e 1,0 — se virasse maior que 1 o limite morderia menos e a
    desigualdade continuaria no mesmo sentido.
    """
    rng = np.random.default_rng(7)
    jk = rng.uniform(0.05, 0.5, size=40)
    real_tn = rng.uniform(0.02, 0.4, size=40)
    row_w = np.full(40, 1.0 / 180.0)
    pedido_cru = float((jk * row_w).sum())
    limitado = np.minimum(jk, delta.PRED_TRACE_CAP_K * real_tn)
    pedido_limitado = float((limitado * row_w).sum())
    assert pedido_limitado <= pedido_cru
    orcamento = 0.3
    assert orcamento / pedido_limitado >= orcamento / pedido_cru
    assert delta.PRED_TRACE_CAP_K > 0.0


def test_laudo_reproduz_o_orcamento_calculado_a_mao_pelo_pacote(submissao):
    """Fim a fim: o laudo bate a composicao direta das mesmas chamadas do pacote.

    Nao e tautologia — o caminho do modulo le o `.h5ad` por fatias, factoriza o
    grupo a partir de `obs` e monta o pseudobulk pelo nosso `bulk_lognorm`,
    enquanto este teste carrega tudo de uma vez pelo `anndata`. Se a leitura por
    fatia desalinhasse grupo e linha, os dois divergiriam.
    """
    laudo = orc.auditar(submissao, "D", chunk=64)
    adata = ad.read_h5ad(submissao)
    adata = adata[adata.obs["context"].astype(str) == "D"]
    rotulos_celula = adata.obs["target_gene"].astype(str).to_numpy()
    codigos, unicos = pd.factorize(rotulos_celula)
    genes = adata.var_names.astype(str).tolist()
    somas = np.zeros((len(unicos), len(genes)), dtype=np.float64)
    matriz = sp.csr_matrix(adata.X)
    for p in range(len(unicos)):
        somas[p] = np.asarray(matriz[codigos == p].sum(axis=0)).ravel()
    medias = orc.bulk_lognorm(somas)
    jk = np.asarray(
        moments.jackknife_correction(matriz, np.asarray(codigos), len(unicos), orc.ALVO_SOMA, chunk=64)
    )
    rotulos = [str(u) for u in unicos]
    excl = delta._exclusion_cols(  # noqa: SLF001
        rotulos, genes, target_gene_map=None, gate_labels=rotulos, who="teste", n_genes=len(genes)
    )
    row_w = delta._row_weights(len(rotulos), len(genes), excl)  # noqa: SLF001
    esperado = float(delta._across_pert_budget(medias, excl, row_w))  # noqa: SLF001

    assert laudo.orcamento == pytest.approx(esperado, rel=1e-12)
    assert laudo.pedido_sem_limite == pytest.approx(float((jk * row_w).sum()), rel=1e-12)
    assert laudo.alvos_pontuados == len(ALVOS)


def test_os_seis_alvos_resolvem_o_proprio_gene(submissao):
    """Os rotulos SAO genes do eixo, como no painel VCC, entao `excl` cobre todos.

    Isso importa porque `row_w` vale `1/(G-1)` na linha que resolve e `1/G` na
    que nao: um painel com resolucao parcial tem duas taxas de cambio, e o
    docstring do `_row_weights` mostra como isso e farmavel. Conferir a cobertura
    e conferir que estamos no caso de taxa unica.
    """
    laudo = orc.auditar(submissao, "D", chunk=64)
    assert laudo.alvos_com_alvo_resolvido == laudo.alvos_pontuados


def test_veredito_diz_paga_quando_a_razao_passa_de_um(submissao, monkeypatch):
    """O veredito e sobre a razao, e as tres faixas tem de ser distinguiveis."""
    laudo = orc.auditar(submissao, "D", chunk=64)
    if laudo.razao_sem_limite >= 1.0:
        assert "paga por inteiro" in laudo.veredito
        assert laudo.r_sem_limite == 1.0
    elif laudo.razao_sem_limite >= laudo.flag:
        assert "descontada" in laudo.veredito
    else:
        assert "EM RISCO" in laudo.veredito
    assert laudo.flag == pytest.approx(0.7)


def test_contexto_inexistente_falha_em_vez_de_pontuar_vazio(submissao):
    """Pedir um contexto ausente tem de parar: pontuar zero celulas daria razao NaN."""
    with pytest.raises(ValueError, match="nenhuma celula no contexto"):
        orc.auditar(submissao, "Z", chunk=64)


def test_grupo_nao_contiguo_da_o_mesmo_laudo(tmp_path):
    """A submissao nossa vem ordenada por grupo; um arquivo de terceiro pode intercalar.

    O modulo le a fatia minima que cobre as linhas do grupo e descarta o resto.
    Se esse recorte errasse, o pseudobulk misturaria alvos e o orcamento sairia
    maior (mais dispersao aparente) sem nada falhar.
    """
    caminho = mundo(tmp_path / "ordenada.h5ad")
    laudo_ordenado = orc.auditar(caminho, "D", chunk=64)

    adata = ad.read_h5ad(caminho)
    rng = np.random.default_rng(5)
    ordem = rng.permutation(adata.n_obs)
    embaralhada = tmp_path / "embaralhada.h5ad"
    adata[ordem].copy().write_h5ad(embaralhada)
    laudo_embaralhado = orc.auditar(embaralhada, "D", chunk=64)

    assert laudo_embaralhado.orcamento == pytest.approx(laudo_ordenado.orcamento, rel=1e-10)
    assert laudo_embaralhado.pedido_sem_limite == pytest.approx(
        laudo_ordenado.pedido_sem_limite, rel=1e-10
    )
