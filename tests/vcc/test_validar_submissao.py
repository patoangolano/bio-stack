"""O validador de contrato reproduz os 15 laudos do validador antigo.

O oraculo sao os `laudo_*.json` resgatados: veredito, quantidade de itens
testados e, para cada item reprovado, o detalhe exato -- que carrega as
contagens (8000 zeros, 400 celulas nulas, a linha da primeira celula afetada)
e por isso prova que a verificacao mediu a coisa certa, nao so que disparou.
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
from conftest import RAIZ, VARIANTES, escrever_forja, laudo_antigo
from validar_submissao import validar

SCRIPT = RAIZ / "scripts" / "vcc" / "validar_submissao.py"


def _sem_caminho(laudo: dict) -> dict:
    return {k: v for k, v in laudo.items() if k != "arquivo"}


def test_nomes_e_ordem_dos_28_itens(forja, ref):
    laudo = validar(forja("limpo"), ref, contextos=["D"])
    esperado = [i["nome"] for i in laudo_antigo("limpo")["itens"]]
    assert [i.nome for i in laudo.itens] == esperado


@pytest.mark.parametrize("variante", VARIANTES)
def test_veredito_reproduz_o_laudo_antigo(forja, ref, variante):
    laudo = validar(forja(variante), ref, contextos=["D"])
    assert laudo.veredito == laudo_antigo(variante)["veredito"]


@pytest.mark.parametrize("variante", VARIANTES)
def test_itens_reprovados_e_detalhes_reproduzem_o_laudo_antigo(forja, ref, variante):
    laudo = validar(forja(variante), ref, contextos=["D"])
    obtido = {i.nome: i.detalhe for i in laudo.itens if not i.passou}
    esperado = {i["nome"]: i["detalhe"] for i in laudo_antigo(variante)["itens"] if not i["passou"]}
    assert obtido == esperado


@pytest.mark.parametrize("variante", VARIANTES)
def test_quantidade_de_itens_testados_reproduz_o_laudo_antigo(forja, ref, variante):
    # 18, nao 28, quando o indptr esta quebrado: o mapeamento linha -> valor e
    # indefinido e ler os valores daria numeros sem sentido.
    laudo = validar(forja(variante), ref, contextos=["D"])
    assert laudo.itens_testados == laudo_antigo(variante)["itens_testados"]


@pytest.mark.parametrize(
    "variante",
    [
        "limpo",
        "desordem",
        "repetido",
        "indice_fora",
        "celula_nula",
        "pert_nula",
        "fracionario",  # primeira linha afetada 5: bloco 0
        "negativo",  # linha 10: bloco 1
        "nan",  # linha 15: bloco 2
    ],
)
def test_blocos_pequenos_que_cortam_linhas_dao_o_mesmo_laudo(forja, ref, variante):
    # A submissao real tem ~2,1e9 nao nulos e e lida em blocos de linhas. As
    # verificacoes por linha (ordem, soma, primeira linha afetada) sao as que
    # quebram na fronteira entre blocos; 7 linhas por bloco forca muitas.
    # `negativo` e `nan` sao obrigatorias: o defeito delas cai FORA do primeiro
    # bloco, unico jeito de pegar um deslocamento de bloco esquecido -- no bloco
    # padrao de 4.096 linhas todo defeito das falsificacoes cai no bloco 0, e o
    # teste de mutacao mostrou que, sem elas, esse erro passava nos 58 testes.
    pequeno = validar(forja(variante), ref, contextos=["D"], linhas_por_bloco=7)
    grande = validar(forja(variante), ref, contextos=["D"], linhas_por_bloco=1_000_000)
    assert _sem_caminho(pequeno.para_dict()) == _sem_caminho(grande.para_dict())


def test_contextos_padrao_vem_do_manifest_e_exigem_o_painel_completo(forja, ref):
    # Sem `contextos=`, um arquivo so com D nao pode passar: o painel final e D/E/F.
    laudo = validar(forja("limpo"), ref)
    assert laudo.contextos_exigidos == ["D", "E", "F"]
    assert laudo.painel_completo is True
    assert {"contextos", "total_de_celulas", "alvos_por_contexto"} <= set(laudo.falhas)


def test_sem_coluna_context_reprova_sem_quebrar(tmp_path, ref):
    laudo = validar(escrever_forja("sem_contexto", tmp_path / "s.h5ad"), ref, contextos=["D"])
    assert laudo.veredito == "reprovado"
    assert {"colunas_obs", "contextos", "alvos_por_contexto", "celulas_por_alvo"} <= set(
        laudo.falhas
    )


def test_x_denso_reprova_e_nao_le_estrutura_csr(tmp_path, ref):
    genes = pd.read_csv(ref / "gene_names.csv")["gene_name"].astype(str).tolist()
    obs = pd.DataFrame(
        {"target_gene": ["ABCC1"] * 4, "context": ["D"] * 4},
        index=pd.Index([f"c{i}" for i in range(4)], dtype=object),
    )
    x = np.ones((4, len(genes)), dtype=np.float32)
    caminho = tmp_path / "denso.h5ad"
    ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=pd.Index(genes, dtype=object))).write_h5ad(
        caminho
    )
    laudo = validar(caminho, ref, contextos=["D"])
    assert laudo.veredito == "reprovado"
    assert "matriz_csr" in laudo.falhas
    assert laudo.itens_testados == 12


@pytest.mark.parametrize(("variante", "codigo"), [("limpo", 0), ("nan", 1)])
def test_cli_sai_com_codigo_que_o_portao_de_envio_le(forja, ref, tmp_path, variante, codigo):
    saida = tmp_path / "laudo.json"
    r = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(forja(variante)),
            "--ref",
            str(ref),
            "--contextos",
            "D",
            "--json",
            str(saida),
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    assert r.returncode == codigo, r.stderr
    gravado = json.loads(Path(saida).read_text(encoding="utf-8"))
    assert gravado["veredito"] == laudo_antigo(variante)["veredito"]
