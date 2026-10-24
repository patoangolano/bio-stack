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


def montar_mundo(tmp_path, monkeypatch, escala_do_desvio: float, genes: list[str] = GENES):
    """Um mundo pequeno com a FORMA do real, parametrizado na magnitude do desvio.

    `escala_do_desvio` existe porque duas propriedades diferentes vivem em
    regimes diferentes. A ponte do formato e a reproducao de k=1 valem em
    qualquer magnitude. A proporcionalidade do fator POR GENE, medida sobre o
    delta emitido, só e limpa no regime quase linear: o fechamento composicional
    soma uma constante POR LINHA, e essa constante cresce com o quanto a
    perturbacao desloca a massa. Com desvio de 0,5 em log2 e k=4, sao 4x em
    razao e a constante e grande o bastante para embaralhar a ordem dos genes.
    """
    rng = np.random.default_rng(17)
    transfer, eixo, ctrl = tmp_path / "transfer", tmp_path / "eixo", tmp_path / "ctrl"
    for p in (transfer, eixo, ctrl):
        p.mkdir()
    por_gene = {"gene": genes}
    deltas = {}
    for c in CONTEXTOS:
        w = rng.uniform(0.5, 400.0, size=len(genes))
        w[-1] = 0.0  # um gene inalcancavel, como os 2.975 de D
        w = w / w.sum() * 50_000.0
        b = np.log1p(w)
        por_gene[f"b_ctrl_{c}"] = b
        por_gene[f"no_universo_{c}"] = w > 100.0
        media = rng.normal(0.0, 0.3 * min(escala_do_desvio / 0.5, 1.0), size=len(genes))
        log2r = np.tile(media, (len(ALVOS), 1))
        log2r[0] += rng.normal(0.0, escala_do_desvio, size=len(genes))  # desvio especifico
        log2r[1] += rng.normal(0.0, escala_do_desvio, size=len(genes))
        deltas[f"delta_{c}"] = vk.para_delta(log2r, np.expm1(b))
    np.savez_compressed(
        transfer / "delta_previsto.npz",
        genes=np.array(genes, dtype=object),
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
    # Peso por gene como o `encolhimento.py` grava: media geometrica 1 sobre os
    # medidos, e os quatro ultimos genes SEM medida (entradas = 0), que e o caso
    # de 37% do painel real.
    rng_peso = np.random.default_rng(5)
    medidos = len(genes) - 4
    log_peso = rng_peso.normal(0.0, 0.45, size=medidos)
    log_peso -= log_peso.mean()
    peso = np.concatenate([np.exp(log_peso), np.ones(4)])
    entradas = np.concatenate([np.full(medidos, 800.0), np.zeros(4)])
    np.savez(
        transfer / "pesos_por_gene.npz",
        genes=np.array(genes, dtype=str),
        peso=peso,
        entradas=entradas,
        k_global_medido=np.float64(5.615),
    )
    return transfer, deltas


@pytest.fixture
def entradas(tmp_path, monkeypatch):
    """O mundo de sempre: desvio de 0,5 em log2, que e a ordem do real sob k=4."""
    return montar_mundo(tmp_path, monkeypatch, escala_do_desvio=0.5)


GENES_LARGO = [f"G{i}" for i in range(400)]


@pytest.fixture
def entradas_lineares(tmp_path, monkeypatch):
    """O mesmo mundo em regime quase linear e com 400 genes.

    Duas contaminacoes somam na medida fim a fim, e as duas caem com o numero de
    genes. A constante de fechamento vale `log2(massa_alvo/massa)` e e dominada
    pela fracao de massa do gene mais abundante. O residuo da centragem vale
    `k*mean_g(fator*desvio)`, que e uma media sobre genes de valores de media
    zero e cai com `1/sqrt(n_genes)`. Com os 12 genes do mundo padrao as duas
    ficam da ORDEM do efeito que se quer medir; com 400 caem abaixo de 2%.
    """
    return montar_mundo(tmp_path, monkeypatch, escala_do_desvio=0.01, genes=GENES_LARGO)


def test_saida_carrega_sem_pickle_no_formato_delta(entradas):
    transfer, _ = entradas
    vk.main([])
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
    vk.main([])
    d = ler_delta(transfer / "varredura_k" / "delta_kd1.npz")
    for c in CONTEXTOS:
        assert np.allclose(d.por_contexto[c], deltas[f"delta_{c}"], rtol=0, atol=2e-6)


def test_k4_e_factivel_e_mexe_so_no_desvio(entradas):
    """`log1p(.) - b >= -b` por construcao, e os alvos sem evidencia nao mudam com k."""
    transfer, _ = entradas
    vk.main([])
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


# --------------------------------------------------------------------------
# Encolhimento POR GENE: `alpha` redistribui a mesma massa, e so isso.
# --------------------------------------------------------------------------


def test_alpha_zero_reproduz_o_braco_escalar_bit_a_bit(entradas):
    """O invariante que torna a varredura de alpha interpretavel.

    Em `alpha = 0` o fator por gene e 1 e a constante de variancia e 1, entao o
    braco TEM de sair identico ao escalar — nao proximo, identico. Se saisse
    diferente, a varredura mediria a diferenca entre dois caminhos de codigo em
    vez de medir alpha. E o mesmo papel de `k=1` reproduzir o `delta_previsto`.
    """
    transfer, _ = entradas
    vk.main([])
    escalar = ler_delta(transfer / "varredura_k" / "delta_kd4.npz")
    vk.main(["--k", "4", "--alphas", "0"])
    por_gene = ler_delta(transfer / "varredura_k" / "delta_kd4_a0.npz")
    for c in CONTEXTOS:
        assert np.array_equal(escalar.por_contexto[c], por_gene.por_contexto[c])


def test_fator_por_gene_preserva_a_soma_de_quadrados_e_segue_o_peso():
    """`c(alpha)` existe para que alpha nao mexa no nivel. Confere o que promete.

    E a ordem: o fator e monotono no peso, que e a unica coisa que alpha faz.
    Subir `k` move tudo; alpha troca QUEM move, mantendo o total.
    """
    rng = np.random.default_rng(23)
    desvio = rng.normal(0.0, 0.4, size=(5, len(GENES)))
    peso = np.exp(rng.normal(0.0, 0.5, size=len(GENES)))
    referencia = float(np.square(desvio).sum())
    for alpha in (0.0, 0.25, 0.5, 1.0, 2.0):
        fator, c = vk.fator_por_gene(peso, alpha, desvio)
        escalado = float(np.square(fator[None, :] * desvio).sum())
        assert escalado == pytest.approx(referencia, rel=1e-12)
        if alpha == 0.0:
            assert c == 1.0
            assert np.array_equal(fator, np.ones(len(GENES)))
        else:
            # Ordem estrita: argsort do fator e argsort do peso.
            assert np.array_equal(np.argsort(fator), np.argsort(peso))
            # E alpha maior afasta mais do neutro, em log.
            assert np.abs(np.log(fator)).sum() > 0.0


def test_alpha_maior_espalha_mais_sem_mover_o_total():
    """Alpha e um eixo monotono de dispersao, nao um liga-desliga."""
    rng = np.random.default_rng(29)
    desvio = rng.normal(0.0, 0.4, size=(7, len(GENES)))
    peso = np.exp(rng.normal(0.0, 0.5, size=len(GENES)))
    espalhamento = []
    for alpha in (0.0, 0.25, 0.5, 0.75, 1.0):
        fator, _ = vk.fator_por_gene(peso, alpha, desvio)
        espalhamento.append(float(np.std(np.log(fator))))
    assert espalhamento == sorted(espalhamento)
    assert espalhamento[0] == 0.0
    assert espalhamento[-1] > 0.3


def test_alpha_um_redistribui_na_ordem_do_peso(entradas_lineares):
    """Fim a fim: o delta emitido carrega o fator por gene que foi pedido.

    Medido em LOG2 DA RAZAO, que e onde alpha se aplica. A primeira versao deste
    teste mediu em espaco de delta e a correlacao saiu -0,34: `log1p` e quase
    linear no gene abundante e quase identidade na proporcao do gene raro, e os
    dois regimes convivem no mesmo eixo.

    `para_log2r(para_delta(L)) = L + log2(massa_alvo/massa)` — uma constante POR
    LINHA, do fechamento composicional. Ela e removida centrando a linha, e e
    por isso que o mundo aqui e o LINEAR: a centragem subtrai uma constante
    diferente em cada braco, e o residuo disso cresce com a magnitude do desvio.
    No mundo de desvio 0,5 sob k=4 — 4x em razao — esse residuo leva a
    correlacao de Pearson a 0,89 e o Spearman a 0,79; em regime linear o fator
    reaparece com dois digitos de precisao, que e o que se afirma.
    """
    transfer, _ = entradas_lineares
    vk.main(["--k", "4", "--alphas", "0,1"])
    b = pd.read_csv(vk.CTRL / "por_gene.tsv", sep="\t")
    a0 = ler_delta(transfer / "varredura_k" / "delta_kd4_a0.npz")
    a1 = ler_delta(transfer / "varredura_k" / "delta_kd4_a1.npz")
    with np.load(transfer / "pesos_por_gene.npz", allow_pickle=False) as z:
        peso = np.asarray(z["peso"])
        tem_medida = np.asarray(z["entradas"]) > 0

    for c in CONTEXTOS:
        proporcao = np.expm1(np.asarray(b[f"b_ctrl_{c}"], dtype=np.float64))
        alcancavel = proporcao > 0
        magnitudes = []
        for arm in (a0, a1):
            log2r = vk.para_log2r(arm.por_contexto[c].astype(np.float64), proporcao)
            # Linha 2 (G2) nao tem desvio: carrega so a parcela media.
            desvio = log2r[:2] - log2r[2][None, :]
            desvio -= desvio[:, alcancavel].mean(axis=1, keepdims=True)
            magnitudes.append(np.sqrt(np.square(desvio).sum(axis=0)))
        d0, d1 = magnitudes
        # O fator PEDIDO vem do desvio ORIGINAL, que e o do braco alpha=0 a menos
        # do `k` — e `c` e invariante a escala, entao o `k` cancela. Tirar `c` do
        # desvio ja escalado do braco alpha=1 daria um fator errado por uma
        # constante global: foi o que a primeira versao deste teste fez, e a
        # razao medida/pedido saiu em 1,45 uniforme.
        log2r_a0 = vk.para_log2r(a0.por_contexto[c].astype(np.float64), proporcao)
        desvio_a0 = log2r_a0[:2] - log2r_a0[2][None, :]
        f0, _ = vk.fator_por_gene(peso, 0.0, desvio_a0)
        f1, _ = vk.fator_por_gene(peso, 1.0, desvio_a0)
        move = alcancavel & tem_medida & (d0 > 1e-9)
        medido = d1[move] / d0[move]
        pedido = f1[move] / f0[move]
        # A ordem e exata: correlacao em log acima de 0,99 sobre todos os genes
        # que se movem.
        assert np.corrcoef(np.log(medido), np.log(pedido))[0, 1] > 0.99
        # O VALOR casa com mediana de erro de 1 a 2%.
        erro = np.abs(medido / pedido - 1.0)
        assert np.median(erro) < 0.03
        # E o que prova que o residuo e a CENTRAGEM e nao um fator errado: o erro
        # cai com a magnitude do gene. Um termo aditivo por linha pesa em
        # proporcao inversa ao desvio; um fator errado pesaria igual em todos.
        # Medido: mediana de 1,1-1,8% em todos, 0,5-0,9% no quartil de cima.
        por_faixa = [
            float(np.median(erro[d0[move] >= np.quantile(d0[move], q)]))
            for q in (0.0, 0.5, 0.75)
        ]
        assert por_faixa == sorted(por_faixa, reverse=True)
        # Gene sem medida tem peso 1 e sente apenas a constante, nao a forma:
        # esses genes andam JUNTOS. A afirmacao e relativa porque o mesmo residuo
        # aditivo por linha os separa em alguns por cento — o que se nega e que
        # eles se espalhem como os medidos, e o contraste e de uma ordem de
        # grandeza.
        neutro = alcancavel & ~tem_medida & (d0 > 1e-9)
        if neutro.sum() >= 2:
            espalhamento_neutro = float(np.ptp(d1[neutro] / d0[neutro]))
            assert espalhamento_neutro < 0.1
            assert espalhamento_neutro < 0.1 * float(np.ptp(medido))


def test_alpha_exige_o_arquivo_de_pesos(entradas, monkeypatch):
    """Sem `encolhimento.py` rodado, a varredura para em vez de inventar peso 1."""
    transfer, _ = entradas
    (transfer / "pesos_por_gene.npz").unlink()
    with pytest.raises(FileNotFoundError, match="encolhimento"):
        vk.main(["--alphas", "1"])
