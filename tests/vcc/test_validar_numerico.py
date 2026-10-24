"""O validador numerico pega submissao que passa no contrato com os numeros errados.

Mundo sintetico pequeno (60 genes, 12 alvos, 2 contextos) e composicionalmente
coerente: o delta pretendido ja sai renormalizado, entao uma submissao amostrada
dele, sem defeito, tem de dar inclinacao ~1. Os defeitos sao os do registro:
alvos no controle (a v1 tinha 143 de 300 assim e passou 28/28 no contrato), o
stub de delta zero, e conferir contra o `--delta` da versao errada.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from conftest import RAIZ
from formato_delta import Delta, ler_delta, salvar_delta
from pseudobulk import bulk_lognorm
from validar_numerico import FAIXA, medir

SCRIPT = RAIZ / "scripts" / "vcc" / "validar_numerico.py"
PY_CELL_EVAL2 = RAIZ / "envs" / "cell-eval2" / ".venv" / "Scripts" / "python.exe"

NG = 60
GENES = [f"G{i}" for i in range(NG)]
ALVOS = GENES[:12]  # alvos estao no eixo de genes, como no desafio
CONTEXTOS = ["D", "E"]
CELULAS = 300
PROFUNDIDADE = 4000


def _escrever(caminho: Path, contagens: np.ndarray, alvos: list[str], contextos: list[str]) -> Path:
    obs = pd.DataFrame(
        {"target_gene": pd.Categorical(alvos), "context": pd.Categorical(contextos)},
        index=pd.Index([f"c{i}" for i in range(len(alvos))], dtype=object),
    )
    x = sp.csr_matrix(contagens.astype(np.float32))
    ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=pd.Index(GENES, dtype=object))).write_h5ad(caminho)
    return caminho


def _composicional(b_ctrl: np.ndarray, log_razao: np.ndarray) -> np.ndarray:
    """Delta em bulk_lognorm a partir de uma razao por gene, renormalizado para somar 5e4."""
    m = np.expm1(b_ctrl) * np.exp(log_razao)
    m = m * (5e4 / m.sum(axis=-1, keepdims=True))
    return np.log1p(m) - b_ctrl


class Mundo:
    """Controles, delta pretendido e um amostrador de submissao, por contexto."""

    def __init__(self, pasta: Path, massa_dos_alvos: float = 1.0):
        self.rng = np.random.default_rng(11)
        self.pasta_ctrl = pasta / "controles"
        self.pasta_ctrl.mkdir()
        self.b_ctrl: dict[str, np.ndarray] = {}
        self.delta: dict[str, np.ndarray] = {}
        for c in CONTEXTOS:
            p0 = self.rng.gamma(0.6, size=NG)
            p0[: len(ALVOS)] *= massa_dos_alvos
            p0 /= p0.sum()
            ctrl = self.rng.multinomial(PROFUNDIDADE, p0, size=3000)
            _escrever(self.pasta_ctrl / f"context_{c}.h5ad", ctrl, ["non-targeting"] * 3000, [c] * 3000)
            self.b_ctrl[c] = bulk_lognorm(ctrl.sum(axis=0, keepdims=True))[0]
            log_razao = self.rng.normal(0, 0.4, size=(len(ALVOS), NG))
            log_razao[np.arange(len(ALVOS)), np.arange(len(ALVOS))] = np.log(0.1)  # knockdown do alvo
            self.delta[c] = _composicional(self.b_ctrl[c], log_razao)

    def gravar_delta(self, caminho: Path, escala: float = 1.0) -> Path:
        salvar_delta(
            caminho,
            Delta(GENES, ALVOS, CONTEXTOS, {c: (escala * d).astype(np.float32) for c, d in self.delta.items()}),
        )
        return caminho

    def submissao(
        self,
        caminho: Path,
        emitido: dict[str, np.ndarray],
        deslocamento_de_base: dict[str, np.ndarray] | None = None,
    ) -> Path:
        """Amostra CELULAS por (contexto, alvo) da composicao `b_ctrl + emitido`.

        `deslocamento_de_base` simula um gerador que partiu da linha de base errada.
        """
        blocos, alvos, contextos = [], [], []
        for c in CONTEXTOS:
            base = self.b_ctrl[c] + (deslocamento_de_base[c] if deslocamento_de_base else 0.0)
            for i, a in enumerate(ALVOS):
                m = np.expm1(base + emitido[c][i])
                blocos.append(self.rng.multinomial(PROFUNDIDADE, m / m.sum(), size=CELULAS))
                alvos += [a] * CELULAS
                contextos += [c] * CELULAS
        return _escrever(caminho, np.vstack(blocos), alvos, contextos)


@pytest.fixture
def mundo(tmp_path) -> Mundo:
    return Mundo(tmp_path)


def test_submissao_fiel_ao_delta_aprova_com_inclinacao_perto_de_1(mundo, tmp_path):
    sub = mundo.submissao(tmp_path / "fiel.h5ad", mundo.delta)
    medida = medir(sub, mundo.gravar_delta(tmp_path / "d.npz"), mundo.pasta_ctrl)
    assert medida.veredito == "aprovado"
    assert all(abs(a.inclinacao - 1) < 0.05 for a in medida.ajustes)


def test_stub_de_delta_zero_reprova(mundo, tmp_path):
    zero = {c: np.zeros_like(d) for c, d in mundo.delta.items()}
    sub = mundo.submissao(tmp_path / "stub.h5ad", zero)
    medida = medir(sub, mundo.gravar_delta(tmp_path / "d.npz"), mundo.pasta_ctrl)
    assert medida.veredito == "reprovado"
    assert all(abs(a.inclinacao) < 0.1 for a in medida.ajustes)


def test_metade_dos_alvos_no_controle_reprova(mundo, tmp_path):
    # O defeito da v1: parte dos alvos emitida como controle passa no contrato.
    meio = {c: np.where(np.arange(len(ALVOS))[:, None] < 6, d, 0.0) for c, d in mundo.delta.items()}
    sub = mundo.submissao(tmp_path / "meio.h5ad", meio)
    medida = medir(sub, mundo.gravar_delta(tmp_path / "d.npz"), mundo.pasta_ctrl)
    assert medida.veredito == "reprovado"
    assert all(a.inclinacao < FAIXA[0] for a in medida.ajustes)


def test_conferir_contra_o_delta_da_versao_errada_reprova(mundo, tmp_path):
    # Arquivo bom, `--delta` de outra versao (k diferente): tem de reprovar, nao passar calado.
    sub = mundo.submissao(tmp_path / "fiel.h5ad", mundo.delta)
    medida = medir(sub, mundo.gravar_delta(tmp_path / "metade.npz", escala=0.5), mundo.pasta_ctrl)
    assert medida.veredito == "reprovado"
    assert all(a.inclinacao > FAIXA[1] for a in medida.ajustes)


def test_gene_alvo_fica_fora_da_regressao(tmp_path):
    # O gene-alvo e excluido de todas as metricas do placar; aqui ele vem SUPERexpresso
    # em vez de derrubado, e isso nao pode mover a inclinacao.
    #
    # Os alvos tem massa pequena nos controles DE PROPOSITO. Com massa normal, o
    # alvo superexpresso toma em mediana 11% da soma (ate 66%), e a normalizacao
    # para 5e4 empurra todos os outros genes para baixo: intercepto -0,25 e
    # inclinacao 0,94 no contexto D, medido. Isso e vazamento composicional real
    # -- o gene-alvo sai das metricas mas NAO sai da normalizacao -- e nao e o que
    # este teste mede. Com massa pequena, sobra so a exclusao: se o alvo entrasse
    # na regressao, 12 pontos com x ~ -1,4 e y = +3 derrubariam a inclinacao.
    mundo = Mundo(tmp_path, massa_dos_alvos=0.02)
    errado = {c: d.copy() for c, d in mundo.delta.items()}
    for d in errado.values():
        d[np.arange(len(ALVOS)), np.arange(len(ALVOS))] = 3.0
    sub = mundo.submissao(tmp_path / "alvo.h5ad", errado)
    medida = medir(sub, mundo.gravar_delta(tmp_path / "d.npz"), mundo.pasta_ctrl)
    assert medida.veredito == "aprovado"
    assert all(abs(a.inclinacao - 1) < 0.05 for a in medida.ajustes)


def test_linha_de_base_errada_no_gerador_reprova(mundo, tmp_path):
    # A armadilha registrada: o `ctrl_means` do PIE e media de log1p POR CELULA, nao
    # log1p do pseudobulk. Um gerador que parte dessa linha de base desloca cada gene
    # por uma constante sem relacao com o delta. A inclinacao nao ve isso -- o
    # deslocamento e descorrelacionado do x -- e o arquivo sai errado em todo gene.
    desvio = {c: np.random.default_rng(5).normal(0, 0.3, NG) for c in CONTEXTOS}
    sub = mundo.submissao(tmp_path / "base.h5ad", mundo.delta, deslocamento_de_base=desvio)
    medida = medir(sub, mundo.gravar_delta(tmp_path / "d.npz"), mundo.pasta_ctrl)
    assert medida.veredito == "reprovado"


def test_cada_contexto_e_julgado_separado(mundo, tmp_path):
    misto = {"D": mundo.delta["D"], "E": np.zeros_like(mundo.delta["E"])}
    sub = mundo.submissao(tmp_path / "misto.h5ad", misto)
    medida = medir(sub, mundo.gravar_delta(tmp_path / "d.npz"), mundo.pasta_ctrl)
    assert medida.veredito == "reprovado"
    assert {a.contexto: a.passou for a in medida.ajustes} == {"D": True, "E": False}


def test_blocos_pequenos_dao_a_mesma_medida(mundo, tmp_path):
    sub = mundo.submissao(tmp_path / "fiel.h5ad", mundo.delta)
    d = mundo.gravar_delta(tmp_path / "d.npz")
    pequeno = medir(sub, d, mundo.pasta_ctrl, linhas_por_bloco=7)
    grande = medir(sub, d, mundo.pasta_ctrl, linhas_por_bloco=1_000_000)
    for a, b in zip(pequeno.ajustes, grande.ajustes, strict=True):
        assert a.inclinacao == pytest.approx(b.inclinacao, rel=1e-9)


def test_ler_delta_recusa_arquivo_que_exige_pickle(tmp_path):
    caminho = tmp_path / "objeto.npz"
    np.savez(
        caminho,
        genes=np.array(GENES, dtype=object),
        alvos=np.array(ALVOS, dtype=object),
        contextos=np.array(["D"], dtype=object),
        delta_D=np.zeros((len(ALVOS), NG), np.float32),
    )
    with pytest.raises(ValueError, match="pickle"):
        ler_delta(caminho)


def test_salvar_e_ler_delta_preserva_tudo(mundo, tmp_path):
    lido = ler_delta(mundo.gravar_delta(tmp_path / "d.npz"))
    assert (lido.genes, lido.alvos, lido.contextos) == (GENES, ALVOS, CONTEXTOS)
    for c in CONTEXTOS:
        np.testing.assert_array_equal(lido.por_contexto[c], mundo.delta[c].astype(np.float32))


@pytest.mark.skipif(not PY_CELL_EVAL2.exists(), reason="envs/cell-eval2 ausente")
def test_bulk_lognorm_identico_ao_do_cell_eval2_instalado(tmp_path):
    rng = np.random.default_rng(3)
    somas = rng.poisson(rng.gamma(0.5, 30, size=(40, 500))).astype(np.float64)
    somas[5] = 0  # grupo todo nulo: o cell-eval2 devolve zeros
    np.save(tmp_path / "somas.npy", somas)
    codigo = (
        "import sys, numpy as np; from cell_eval2.prep import bulk_lognorm_means; "
        "np.save(sys.argv[2], bulk_lognorm_means(np.load(sys.argv[1]), 5e4))"
    )
    subprocess.run(
        [str(PY_CELL_EVAL2), "-I", "-c", codigo, str(tmp_path / "somas.npy"), str(tmp_path / "ref.npy")],
        check=True,
    )
    np.testing.assert_array_equal(bulk_lognorm(somas), np.load(tmp_path / "ref.npy"))


@pytest.mark.parametrize(("defeituosa", "codigo"), [(False, 0), (True, 1)])
def test_cli_sai_com_codigo_que_o_portao_de_envio_le(mundo, tmp_path, defeituosa, codigo):
    emitido = {c: np.zeros_like(d) for c, d in mundo.delta.items()} if defeituosa else mundo.delta
    sub = mundo.submissao(tmp_path / "s.h5ad", emitido)
    saida = tmp_path / "medida.json"
    r = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(sub),
            "--delta",
            str(mundo.gravar_delta(tmp_path / "d.npz")),
            "--controles",
            str(mundo.pasta_ctrl),
            "--json",
            str(saida),
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    assert r.returncode == codigo, r.stderr
    assert json.loads(saida.read_text(encoding="utf-8"))["veredito"] == (
        "reprovado" if defeituosa else "aprovado"
    )
