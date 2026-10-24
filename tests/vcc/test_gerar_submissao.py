"""O gerador de ponta a ponta: do arquivo de delta ao `.h5ad` que passa os dois validadores.

Os testes rodam num painel PEQUENO -- 2 contextos, 40 alvos, 300 genes, 20
celulas por grupo -- por uma razao de desenho, nao de conveniencia: o painel
real tem 360.000 celulas e 15 GB, e um teste que precisa dele nao roda. O que o
painel pequeno nao testa (o tempo, o pico de memoria, o teto de 2^31-1 no nnz)
esta nas guardas numericas do proprio `gerar_submissao`, medidas em producao.

O portao que importa e o T16: `validar_numerico` contra o delta EFETIVO tem de
aprovar. O delta efetivo e a predicao que a submissao realmente faz, e e ele que
o gerador grava ao lado da saida -- nao o delta pedido, que pode conter o
impossivel.
"""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from conftest import escrever_controle_pequeno
from formato_delta import Delta, ler_delta, salvar_delta
from gerar_submissao import gerar_submissao, main

N_ALVOS = 40
N_GENES = 300
N_CEL = 20
CONTEXTOS = ("D", "E")


@pytest.fixture(scope="module")
def painel(tmp_path_factory):
    """Painel pequeno completo: controles por contexto + arquivo de delta."""
    pasta = tmp_path_factory.mktemp("painel")
    controles = pasta / "controles"
    controles.mkdir()
    for c in CONTEXTOS:
        escrever_controle_pequeno(controles / f"context_{c}.h5ad", n_celulas=120, n_genes=N_GENES)
        # `escrever_controle_pequeno` rotula o contexto como D; reescreve o rotulo.
        if c != "D":
            a = ad.read_h5ad(controles / f"context_{c}.h5ad")
            a.obs["context"] = c
            a.write_h5ad(controles / f"context_{c}.h5ad")

    genes = [f"G{i:04d}" for i in range(N_GENES)]
    alvos = genes[:N_ALVOS]  # alvos SAO genes do eixo, como na competicao
    # Referencia do painel pequeno, no mesmo formato do `controls-final`: e o que
    # permite rodar o validador de CONTRATO aqui. Ele confere o painel contra
    # esses tres arquivos, nao contra constantes no codigo.
    pd.DataFrame({"gene_name": genes}).to_csv(controles / "gene_names.csv", index=False)
    pd.DataFrame({"target_gene": alvos}).to_csv(controles / "pert_counts.csv", index=False)
    (controles / "manifest.json").write_text(
        json.dumps(
            {
                "season": "2026",
                "partition": "test",
                "contexts": list(CONTEXTOS),
                "pert_col": "target_gene",
                "context_col": "context",
                "control_label": "non-targeting",
                "n_genes": N_GENES,
                "n_constructs": N_ALVOS,
                "cells_per_pert": N_CEL,
            }
        ),
        encoding="utf-8",
    )
    rng = np.random.default_rng(21)
    por_contexto = {}
    for c in CONTEXTOS:
        d = np.zeros((N_ALVOS, N_GENES), dtype=np.float32)
        # Esparso como o real: 2% dos pares com delta, e o alvo sempre derrubado.
        for i in range(N_ALVOS):
            onde = rng.choice(N_GENES, size=6, replace=False)
            d[i, onde] = rng.normal(0.0, 0.05, size=6)
            d[i, i] = -0.5
        por_contexto[c] = d
    caminho = pasta / "delta.npz"
    salvar_delta(caminho, Delta(genes, alvos, list(CONTEXTOS), por_contexto))
    return pasta, controles, caminho, genes, alvos


def ler(caminho) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]:
    a = ad.read_h5ad(caminho)
    x, obs, var = a.X, a.obs, a.var
    assert isinstance(x, sp.csr_matrix)
    assert isinstance(obs, pd.DataFrame) and isinstance(var, pd.DataFrame)
    return x, obs, var


def gerar(painel, **kw):
    pasta, controles, delta, _, _ = painel
    saida = pasta / kw.pop("nome", "sub.h5ad")
    return saida, gerar_submissao(
        delta=delta, controles=controles, saida=saida, n_celulas=N_CEL, **kw
    )


# --------------------------------------------------------------------------


def test_gera_o_painel_inteiro_com_a_forma_certa(painel):
    saida, r = gerar(painel)
    _, _, _, genes, alvos = painel
    assert saida.exists()
    x, obs, var = ler(saida)
    assert x.shape == (len(CONTEXTOS) * N_ALVOS * N_CEL, N_GENES)
    assert sorted(set(obs["context"].astype(str))) == sorted(CONTEXTOS)
    assert sorted(set(obs["target_gene"].astype(str))) == sorted(alvos)
    assert var.index.astype(str).tolist() == genes
    # Exatamente `n_celulas` por (contexto, alvo), sem linha de controle.
    contagem = obs.groupby(["context", "target_gene"], observed=True).size()
    assert set(contagem.unique().tolist()) == {N_CEL}
    assert "non-targeting" not in set(obs["target_gene"].astype(str))
    assert r.n_celulas == len(CONTEXTOS) * N_ALVOS * N_CEL
    assert r.nnz_total == x.nnz
    assert len(r.grupos) == len(CONTEXTOS) * N_ALVOS


def test_grava_o_delta_efetivo_ao_lado(painel):
    """O delta efetivo e a predicao que a submissao FAZ, e nao o que foi pedido."""
    saida, r = gerar(painel, nome="sub_ef.h5ad")
    efetivo = saida.with_suffix(".efetivo.npz")
    assert efetivo.exists()
    assert r.delta_efetivo == str(efetivo)
    lido = ler_delta(efetivo)  # carrega sem pickle, por construcao do formato
    _, _, _, genes, alvos = painel
    assert lido.genes == genes
    assert lido.alvos == alvos
    assert sorted(lido.contextos) == sorted(CONTEXTOS)
    for c in CONTEXTOS:
        assert lido.por_contexto[c].shape == (N_ALVOS, N_GENES)
        assert np.all(np.isfinite(lido.por_contexto[c]))


def test_contrato_28_de_28(painel):
    """T15 -- o validador de contrato aprova os 28 itens, sem nenhum numerico."""
    import validar_submissao

    saida, _ = gerar(painel, nome="sub_contrato.h5ad")
    _, controles, _, _, _ = painel
    laudo = validar_submissao.validar(saida, referencia=controles)
    reprovados = [i.nome for i in laudo.itens if not i.passou]
    assert reprovados == []
    assert laudo.veredito == "aprovado"
    assert laudo.itens_testados == 28


def test_validador_numerico_contra_o_delta_efetivo(painel):
    """T16 -- PORTAO VINCULANTE: a inclinacao contra o delta efetivo tem de aprovar."""
    import validar_numerico

    saida, _ = gerar(painel, nome="sub_num.h5ad")
    pasta, controles, _, _, _ = painel
    medida = validar_numerico.medir(saida, saida.with_suffix(".efetivo.npz"), controles)
    assert medida.veredito == "aprovado"
    for a in medida.ajustes:
        assert 0.90 <= a.inclinacao <= 1.10, a
        assert a.deslocamento_por_gene <= validar_numerico.LIMITE_DESLOCAMENTO
        # `razao_deslocamento` ~1 e a assinatura de ruido INDEPENDENTE entre
        # alvos. Um molde de celulas compartilhado levaria a razao a ~raiz(P).
        assert 0.3 <= a.razao_deslocamento <= 3.0, a


def test_delta_infactivel_e_recusado_sem_o_flag(painel):
    """T17 -- o gerador nao "conserta" a inclinacao; ele recusa e diz de quem e o conserto."""
    pasta, controles, _, genes, alvos = painel
    rng = np.random.default_rng(31)
    # Infactivel mas NAO degenerado: delta negativo grande em metade dos genes,
    # ruido pequeno no resto, para que a variancia nao seja zero.
    ruim = {}
    for c in CONTEXTOS:
        d = rng.normal(0.0, 0.02, size=(N_ALVOS, N_GENES)).astype(np.float32)
        d[:, : N_GENES // 2] = -9.0
        ruim[c] = d
    caminho = pasta / "delta_ruim.npz"
    salvar_delta(caminho, Delta(genes, alvos, list(CONTEXTOS), ruim))

    with pytest.raises(ValueError, match="infactivel|inclinacao"):
        gerar_submissao(
            delta=caminho, controles=controles, saida=pasta / "x.h5ad", n_celulas=N_CEL
        )
    # Com o flag, emite e REPORTA -- a decisao e do operador, registrada no laudo.
    relatorio = gerar_submissao(
        delta=caminho,
        controles=controles,
        saida=pasta / "sub_ruim.h5ad",
        n_celulas=N_CEL,
        aceitar_delta_infactivel=True,
    )
    assert (pasta / "sub_ruim.h5ad").exists()
    assert all(v.inclinacao_contra_delta_cru < 0.85 for v in relatorio.viabilidade)
    assert all(v.fracao_de_pares_infactiveis > 0.4 for v in relatorio.viabilidade)


def test_delta_degenerado_tambem_e_recusado(painel):
    """T17, controle -- `NaN < 0,85` e `False`, logo o portao TEM de testar finitude.

    Um delta constante em todos os pares tem variancia zero e inclinacao NaN. E o
    pior arquivo possivel -- nao predis nada -- e um portao escrito como
    `inclinacao < piso` o deixaria passar em silencio.
    """
    pasta, controles, _, genes, alvos = painel
    constante = {c: np.full((N_ALVOS, N_GENES), -9.0, dtype=np.float32) for c in CONTEXTOS}
    caminho = pasta / "delta_constante.npz"
    salvar_delta(caminho, Delta(genes, alvos, list(CONTEXTOS), constante))
    from conversao import viabilidade_do_delta
    from momentos import momentos_do_controle

    v = viabilidade_do_delta(
        momentos_do_controle(controles / "context_D.h5ad"), constante["D"]
    )
    assert not np.isfinite(v.inclinacao_contra_delta_cru)
    with pytest.raises(ValueError, match="infactivel|inclinacao"):
        gerar_submissao(
            delta=caminho, controles=controles, saida=pasta / "y.h5ad", n_celulas=N_CEL
        )


def test_relatorio_tem_a_viabilidade_de_cada_contexto(painel):
    """T18 -- a viabilidade e medida antes de emitir e sai no laudo."""
    pasta, _, _, _, _ = painel
    caminho = pasta / "laudo.json"
    saida, r = gerar(painel, nome="sub_laudo.h5ad", relatorio=caminho)
    assert len(r.viabilidade) == len(CONTEXTOS)
    assert sorted(v.contexto for v in r.viabilidade) == sorted(CONTEXTOS)
    for v in r.viabilidade:
        assert v.inclinacao_contra_delta_cru >= 0.85
    assert caminho.exists()
    laudo = json.loads(caminho.read_text(encoding="utf-8"))
    assert laudo["semente_base"] == 0
    assert len(laudo["viabilidade"]) == len(CONTEXTOS)
    assert len(laudo["grupos"]) == len(CONTEXTOS) * N_ALVOS
    assert laudo["maior_soma_por_celula"] > 0
    # O laudo nao carrega matriz nenhuma: tem de ser pequeno.
    assert caminho.stat().st_size < 400_000


def test_determinismo_byte_a_byte(painel):
    """T21 -- duas execucoes com a mesma semente dao `data` e `indices` identicos."""
    a, _ = gerar(painel, nome="sub_a.h5ad")
    b, _ = gerar(painel, nome="sub_b.h5ad")
    xa, _, _ = ler(a)
    xb, _, _ = ler(b)
    assert np.array_equal(xa.data, xb.data)
    assert np.array_equal(xa.indices, xb.indices)
    c, _ = gerar(painel, nome="sub_c.h5ad", semente_base=7)
    xc, _, _ = ler(c)
    assert not np.array_equal(xa.data, xc.data)


def test_indptr_em_int64_e_nnz_sob_o_teto(painel):
    """T20 -- `indptr` int64 e o nnz conferido contra o teto declarado."""
    import h5py

    saida, r = gerar(painel, nome="sub_dtype.h5ad")
    with h5py.File(saida, "r") as f:
        x = f["X"]
        assert isinstance(x, h5py.Group)
        indptr = x["indptr"]
        assert isinstance(indptr, h5py.Dataset)
        assert indptr.dtype == np.int64
        assert indptr.size == r.n_celulas + 1
    assert r.nnz_total < 2_147_483_647


def test_teto_de_nnz_recusa_antes_de_escrever(painel):
    pasta, controles, delta, _, _ = painel
    with pytest.raises(ValueError, match="nnz"):
        gerar_submissao(
            delta=delta,
            controles=controles,
            saida=pasta / "estouro.h5ad",
            n_celulas=N_CEL,
            teto_de_nnz=1000,
        )


def test_um_contexto_so(painel):
    """`--contextos D` emite um painel de um contexto, que e como o teste cabe no tempo."""
    saida, r = gerar(painel, nome="sub_um.h5ad", contextos=["D"])
    x, obs, _ = ler(saida)
    assert set(obs["context"].astype(str)) == {"D"}
    assert x.shape == (N_ALVOS * N_CEL, N_GENES)
    assert len(r.viabilidade) == 1


def test_cli_devolve_zero_e_grava_o_laudo(painel, tmp_path):
    pasta, controles, delta, _, _ = painel
    saida = tmp_path / "cli.h5ad"
    laudo = tmp_path / "cli.json"
    codigo = main(
        [
            str(delta),
            "--controles", str(controles),
            "--saida", str(saida),
            "--n-celulas", str(N_CEL),
            "--relatorio", str(laudo),
        ]
    )
    assert codigo == 0
    assert saida.exists() and laudo.exists()


def test_gene_alvo_sai_derrubado_no_emitido(painel):
    """O alvo nao e variavel de folga, mas o knockdown pedido tem de aparecer nele."""
    saida, r = gerar(painel, nome="sub_alvo.h5ad")
    por_grupo = {(g.contexto, g.alvo): g for g in r.grupos}
    negativos = [g.log2fc_do_alvo for g in por_grupo.values() if np.isfinite(g.log2fc_do_alvo)]
    assert len(negativos) >= N_ALVOS
    # Pedimos -0,5 em bulk no alvo de todo grupo; o emitido tem de ser negativo.
    assert float(np.median(negativos)) < 0.0


def test_o_delta_gravado_e_o_EFETIVO_e_nao_o_pedido(painel):
    """T16/T17 -- o teste que distingue os dois, que num delta factivel coincidem.

    Com um delta parcialmente INFACTIVEL os dois arquivos divergem, e a divergencia
    tem consequencia medivel: `validar_numerico` aprova contra o efetivo e reprova
    contra o pedido. Num painel em que todo delta e factivel os dois sao iguais e
    gravar o errado nao aparece em teste nenhum -- foi assim que o teste de
    mutacao pegou esta lacuna.
    """
    import validar_numerico

    pasta, controles, _, genes, alvos = painel
    rng = np.random.default_rng(41)
    parcial = {}
    for c in CONTEXTOS:
        d = rng.normal(0.0, 0.03, size=(N_ALVOS, N_GENES)).astype(np.float32)
        # Um quarto dos genes pedindo o impossivel: delta muito abaixo de -b_ctrl.
        d[:, : N_GENES // 4] = -6.0
        parcial[c] = d
    pedido = pasta / "delta_parcial.npz"
    salvar_delta(pedido, Delta(genes, alvos, list(CONTEXTOS), parcial))

    saida = pasta / "sub_parcial.h5ad"
    gerar_submissao(
        delta=pedido,
        controles=controles,
        saida=saida,
        n_celulas=N_CEL,
        aceitar_delta_infactivel=True,
    )
    efetivo = ler_delta(saida.with_suffix(".efetivo.npz"))
    # Os dois arquivos TEM de diferir, e diferir no lugar certo.
    for c in CONTEXTOS:
        assert not np.allclose(efetivo.por_contexto[c], parcial[c])
        atingido = parcial[c][:, : N_GENES // 4]
        assert float(np.abs(efetivo.por_contexto[c][:, : N_GENES // 4] - atingido).max()) > 1.0

    aprovado = validar_numerico.medir(saida, saida.with_suffix(".efetivo.npz"), controles)
    reprovado = validar_numerico.medir(saida, pedido, controles)
    assert aprovado.veredito == "aprovado"
    assert reprovado.veredito == "reprovado"
    assert all(a.inclinacao < validar_numerico.FAIXA[0] for a in reprovado.ajustes)


def test_maior_soma_por_celula_e_lida_do_arquivo(painel):
    """A guarda do teto e medida sobre o que foi ESCRITO, nunca reconstruida.

    No codigo de referencia do organizador a reconstrucao `max(kernel, controle)`
    divergiu do real em 2 de 4 sementes. Aqui o oraculo e o proprio arquivo, lido
    por outro caminho.
    """
    saida, r = gerar(painel, nome="sub_teto.h5ad")
    x, _, _ = ler(saida)
    esperado = float(np.asarray(x.sum(axis=1)).ravel().max())
    assert r.maior_soma_por_celula == pytest.approx(esperado, rel=0, abs=1e-6)
    assert r.celulas_acima_do_teto == 0
    # E coerente com o maximo que os grupos reportaram, sem ser derivada dele.
    dos_grupos = max(g.soma_maxima_por_celula for g in r.grupos)
    assert r.maior_soma_por_celula == pytest.approx(float(dos_grupos), rel=0, abs=1e-6)
