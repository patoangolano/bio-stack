"""A ponte entre o preditor e o gerador: `varredura_k.py` grava no `formato_delta`.

O preditor recuperado gravava `delta_kd<k>.npz` com rotulos em array de objeto,
que so carregam com pickle. O `gerar_submissao` le pelo `formato_delta`, que
recusa pickle de proposito. Este teste fixa a ponte e, de quebra, o invariante
que o docstring do `varredura_k` afirma e nada conferia: k=1 reproduz o delta de
entrada, porque a decomposicao media + desvio e exata em log2 da razao.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
import varredura_k as vk
from formato_delta import ler_delta

GENES = [f"G{i}" for i in range(12)]
ALVOS = ["G0", "G1", "G2", "G3"]  # G2 e G3 sem evidencia: carregam so a media
CONTEXTOS = ("D", "E", "F")


@pytest.fixture
def entradas(tmp_path, monkeypatch):
    rng = np.random.default_rng(17)
    transfer, eixo, ctrl = tmp_path / "transfer", tmp_path / "eixo", tmp_path / "ctrl"
    for p in (transfer, eixo, ctrl):
        p.mkdir()
    por_gene = {"gene": GENES}
    deltas = {}
    for c in CONTEXTOS:
        w = rng.uniform(0.5, 400.0, size=len(GENES))
        w[-1] = 0.0  # um gene inalcancavel, como os 2.975 de D
        w = w / w.sum() * 50_000.0
        b = np.log1p(w)
        por_gene[f"b_ctrl_{c}"] = b
        por_gene[f"no_universo_{c}"] = w > 100.0
        media = rng.normal(0.0, 0.3, size=len(GENES))
        log2r = np.tile(media, (len(ALVOS), 1))
        log2r[0] += rng.normal(0.0, 0.5, size=len(GENES))  # desvio especifico
        log2r[1] += rng.normal(0.0, 0.5, size=len(GENES))
        deltas[f"delta_{c}"] = vk.para_delta(log2r, np.expm1(b))
    np.savez_compressed(
        transfer / "delta_previsto.npz",
        genes=np.array(GENES, dtype=object),
        alvos=np.array(ALVOS, dtype=object),
        allow_pickle=True,
        **deltas,
    )
    pd.DataFrame(por_gene).to_csv(ctrl / "por_gene.tsv", sep="\t", index=False)
    (eixo / "resumo.json").write_text(
        json.dumps({"evidencia_doadora": {"alvos_sem_evidencia": ["G2", "G3"]}}), encoding="utf-8"
    )
    monkeypatch.setattr(vk, "TRANSFER", transfer)
    monkeypatch.setattr(vk, "EIXO", eixo)
    monkeypatch.setattr(vk, "CTRL", ctrl)
    monkeypatch.setattr(vk, "SAIDA", transfer / "varredura_k")
    monkeypatch.setattr(vk, "ESCALAS_DESVIO", (1.0, 4.0))
    return transfer, deltas


def test_saida_carrega_sem_pickle_no_formato_delta(entradas):
    transfer, _ = entradas
    vk.main()
    d = ler_delta(transfer / "varredura_k" / "delta_kd4.npz")
    assert d.genes == GENES
    assert d.alvos == ALVOS
    assert d.contextos == list(CONTEXTOS)
    for c in CONTEXTOS:
        assert d.por_contexto[c].shape == (len(ALVOS), len(GENES))
        assert d.por_contexto[c].dtype == np.float32


def test_k1_reproduz_o_delta_de_entrada(entradas):
    """A decomposicao e exata: k=1 devolve o `delta_previsto`, a menos do float32."""
    transfer, deltas = entradas
    vk.main()
    d = ler_delta(transfer / "varredura_k" / "delta_kd1.npz")
    for c in CONTEXTOS:
        assert np.allclose(d.por_contexto[c], deltas[f"delta_{c}"], rtol=0, atol=2e-6)


def test_k4_e_factivel_e_mexe_so_no_desvio(entradas):
    """`log1p(.) - b >= -b` por construcao, e os alvos sem evidencia nao mudam com k."""
    transfer, _ = entradas
    vk.main()
    k1 = ler_delta(transfer / "varredura_k" / "delta_kd1.npz")
    k4 = ler_delta(transfer / "varredura_k" / "delta_kd4.npz")
    b = pd.read_csv(vk.CTRL / "por_gene.tsv", sep="\t")
    for c in CONTEXTOS:
        bc = np.asarray(b[f"b_ctrl_{c}"], dtype=np.float64)
        assert np.all(bc[None, :] + k4.por_contexto[c].astype(np.float64) >= -1e-6)
        # G2 e G3 carregam so a media: identicos em k=1 e k=4.
        assert np.allclose(k1.por_contexto[c][2:], k4.por_contexto[c][2:], atol=2e-6)
        # G0 e G1 tem desvio: k=4 os afasta mais da media que k=1.
        dist1 = np.abs(k1.por_contexto[c][:2] - k1.por_contexto[c][2]).sum()
        dist4 = np.abs(k4.por_contexto[c][:2] - k4.por_contexto[c][2]).sum()
        assert dist4 > 2.0 * dist1
