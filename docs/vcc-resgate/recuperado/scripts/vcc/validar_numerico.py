"""
Confere que as celulas emitidas carregam o delta que pretendemos submeter.

## O buraco que este script fecha

`validar_submissao.py` testa 28 itens de contrato — estrutura do CSR, rotulos,
eixo de genes, tetos. Nenhum deles e numerico. A prova esta no disco:
`runs/vcc/submissao/validacao_stub.json` diz **aprovado, 28 de 28** para o
arquivo de delta zero nos 300 alvos, que e o pior arquivo que este pipeline
consegue produzir.

Com quatro versoes de submissao em disco e duas submissoes por dia, "o arquivo
certo esta formalmente valido" nao e a pergunta. A pergunta e "este arquivo e o
do modelo atual".

## O teste

Para cada (contexto, alvo), recomputar o pseudobulk das 400 celulas emitidas,
normalizar para CP50k, aplicar log1p, e subtrair `b_ctrl`. Isso recupera o
delta que as celulas de fato carregam. Comparar com o delta do `.npz` que se
pretende submeter.

Tres estatisticas, e a primeira e a que pega troca de arquivo:

`inclinacao`  regressao de recuperado sobre pedido, sem intercepto. Deve ser
              ~1. Um arquivo emitido de outro delta da inclinacao proxima de 0.
`correlacao`  entre recuperado e pedido, no universo de DE.
`vies`        media de (recuperado - pedido). Deve ser ~0; desvio sistematico
              aponta erro de normalizacao ou de arredondamento.

A comparacao e restrita ao **universo de DE** de cada contexto: fora dele o
delta e minusculo e o ruido de amostragem de 400 celulas domina, o que
afrouxaria o teste exatamente onde ele precisa ser apertado.

## Limite honesto

Isto NAO verifica que o modelo esta certo — verifica que o arquivo corresponde
ao modelo. Um delta errado emitido fielmente passa aqui, e e para isso que
existe o `avaliador_proxy`.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/validar_numerico.py \
        runs/vcc/submissao/submissao_v4.h5ad \
        --delta runs/vcc/transfer/delta_previsto.npz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

CTRL = Path("runs/vcc/ctrl_bulk/por_gene.tsv")
ESCALA_PSEUDOBULK = 5e4
BLOCO_LINHAS = 20_000

# Limites. A inclinacao e o teste de identidade do arquivo: 0,85 a 1,15 aceita
# ruido de amostragem de 400 celulas e recusa um arquivo emitido de outro delta,
# que da inclinacao perto de zero.
INCLINACAO_MINIMA = 0.85
INCLINACAO_MAXIMA = 1.15
# Vies sistematico tolerado em log1p de CP50k. Acima disso ha erro de
# normalizacao ou de arredondamento, nao ruido.
VIES_MAXIMO = 0.005


def _grupo(no: h5py.Group | h5py.File, nome: str) -> h5py.Group:
    """Subgrupo, com o tipo estreitado: `h5py` devolve uniao com Datatype."""
    alvo = no[nome]
    if not isinstance(alvo, h5py.Group):
        raise TypeError(f"{nome} nao e um grupo HDF5")
    return alvo


def _dados(no: h5py.Group, nome: str) -> h5py.Dataset:
    alvo = no[nome]
    if not isinstance(alvo, h5py.Dataset):
        raise TypeError(f"{nome} nao e um dataset HDF5")
    return alvo


def _ler_rotulos(arquivo: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    """(contexto, alvo) por celula, resolvendo categoria do anndata."""
    obs = _grupo(arquivo, "obs")
    saidas: list[np.ndarray] = []
    for coluna in ("context", "target_gene"):
        no = obs[coluna]
        if isinstance(no, h5py.Group):
            categorias = np.asarray(_dados(no, "categories")).astype(str)
            codigos = np.asarray(_dados(no, "codes"), dtype=np.int64)
            saidas.append(categorias[codigos])
        elif isinstance(no, h5py.Dataset):
            saidas.append(np.asarray(no).astype(str))
        else:
            raise TypeError(f"obs/{coluna} em formato inesperado")
    return saidas[0], saidas[1]


def _sha256(caminho: Path) -> str:
    digestor = hashlib.sha256()
    with caminho.open("rb") as fonte:
        for bloco in iter(lambda: fonte.read(1 << 20), b""):
            digestor.update(bloco)
    return digestor.hexdigest()


def conferir_procedencia(h5ad: Path, delta: Path) -> tuple[bool, str]:
    """O `--delta` passado e o mesmo que o emissor usou? Checa IDENTIDADE.

    O teste de regressao abaixo NAO distingue dois deltas parecidos. Medido: v5
    e v6 diferem so na parcela da media, tem correlacao +0,996, e a inclinacao
    cruzada entre elas da 0,996 — dentro de [0,85; 1,15]. Nenhuma janela estreita
    o suficiente para separa-las sobreviveria ao ruido de 400 celulas por grupo,
    entao apertar o limite e o conserto errado.

    O emissor grava `delta` e `delta_sha256` em `medicao_<nome>.json`. Comparar
    com isso pega a troca por identidade, que e o que a pergunta "este arquivo e
    do modelo atual" realmente exige. Arquivo de medicao ausente ou sem o campo
    nao reprova — e a situacao das emissoes anteriores a este campo — mas vira
    aviso, porque a cobertura entao e so a da regressao.
    """
    medicao = h5ad.parent / f"medicao_{h5ad.stem}.json"
    if not medicao.exists():
        return True, f"sem {medicao.name}: procedencia NAO conferida, so a regressao vale"
    registro = json.loads(medicao.read_text(encoding="utf-8"))
    gravado_sha = registro.get("delta_sha256")
    gravado_caminho = registro.get("delta")
    if gravado_sha:
        atual = _sha256(delta)
        if atual == gravado_sha:
            return True, f"procedencia OK por sha256 ({atual[:12]}...)"
        return False, (
            f"o emissor usou um delta de sha256 {str(gravado_sha)[:12]}... e o "
            f"--delta passado tem {atual[:12]}...; sao arquivos DIFERENTES"
        )
    if gravado_caminho:
        mesmo = Path(str(gravado_caminho)).name == delta.name
        if mesmo:
            return True, f"procedencia OK por nome ({delta.name}); sem sha256 gravado"
        return False, (
            f"o emissor usou {Path(str(gravado_caminho)).name} e o --delta passado "
            f"e {delta.name}; sao arquivos DIFERENTES"
        )
    return True, f"{medicao.name} sem campo de delta: procedencia NAO conferida"


def main() -> int:
    analisador = argparse.ArgumentParser(description=__doc__ or "")
    analisador.add_argument("arquivo", help="h5ad da submissao")
    analisador.add_argument("--delta", required=True, help=".npz do delta pretendido")
    analisador.add_argument("--laudo", default=None)
    argumentos = analisador.parse_args()

    caminho = Path(argumentos.arquivo)
    pacote = np.load(Path(argumentos.delta), allow_pickle=True)
    genes_delta = [str(g) for g in pacote["genes"]]
    alvos_delta = [str(a) for a in pacote["alvos"]]
    indice_alvo = {a: i for i, a in enumerate(alvos_delta)}

    por_gene = pd.read_csv(CTRL, sep="\t")
    if list(por_gene["gene"].astype(str)) != genes_delta:
        print("!! eixo de genes do ctrl_bulk difere do eixo do delta", file=sys.stderr)
        return 1

    print(f"conferindo {caminho.name} contra {Path(argumentos.delta).name}")
    procedencia_ok, procedencia_nota = conferir_procedencia(caminho, Path(argumentos.delta))
    print(f"  [{'ok  ' if procedencia_ok else 'FALHA'}] {procedencia_nota}")
    with h5py.File(caminho, "r") as arquivo:
        contexto_por_celula, alvo_por_celula = _ler_rotulos(arquivo)
        x = _grupo(arquivo, "X")
        dados = _dados(x, "data")
        indices = _dados(x, "indices")
        indptr = np.asarray(_dados(x, "indptr"), dtype=np.int64)
        n_celulas = int(indptr.size) - 1
        n_genes = len(genes_delta)

        chaves = np.char.add(np.char.add(contexto_por_celula, "|"), alvo_por_celula)
        unicas, grupo_por_celula = np.unique(chaves, return_inverse=True)
        soma = np.zeros((unicas.size, n_genes), dtype=np.float64)
        celulas_por_grupo = np.bincount(grupo_por_celula, minlength=unicas.size)
        print(f"  {n_celulas:,} celulas em {unicas.size} grupos (contexto, alvo)")

        for inicio in range(0, n_celulas, BLOCO_LINHAS):
            fim = min(inicio + BLOCO_LINHAS, n_celulas)
            a, b = int(indptr[inicio]), int(indptr[fim])
            if b <= a:
                continue
            valores = np.asarray(dados[a:b], dtype=np.float64)
            colunas = np.asarray(indices[a:b], dtype=np.int64)
            comprimentos = np.diff(indptr[inicio : fim + 1])
            linha_de_entrada = np.repeat(
                grupo_por_celula[inicio:fim], comprimentos.astype(np.int64)
            )
            np.add.at(soma, (linha_de_entrada, colunas), valores)

    if int(celulas_por_grupo.min()) != 400 or int(celulas_por_grupo.max()) != 400:
        print(
            f"!! celulas por grupo entre {celulas_por_grupo.min()} e "
            f"{celulas_por_grupo.max()}; esperado exatamente 400",
            file=sys.stderr,
        )
        return 1

    resultados: dict[str, dict[str, float]] = {}
    reprovou = False
    for contexto in ("D", "E", "F"):
        chave_delta = f"delta_{contexto}"
        if chave_delta not in pacote.files:
            continue
        b_ctrl = np.asarray(por_gene[f"b_ctrl_{contexto}"], dtype=np.float64)
        universo = np.asarray(por_gene[f"no_universo_{contexto}"], dtype=bool)
        pedido_todos = np.asarray(pacote[chave_delta], dtype=np.float64)

        recuperados: list[np.ndarray] = []
        pedidos: list[np.ndarray] = []
        for posicao, chave in enumerate(unicas):
            ctx, alvo = str(chave).split("|", 1)
            if ctx != contexto:
                continue
            total = soma[posicao].sum()
            if total <= 0:
                continue
            b_emitido = np.log1p(ESCALA_PSEUDOBULK * soma[posicao] / total)
            recuperados.append((b_emitido - b_ctrl)[universo])
            pedidos.append(pedido_todos[indice_alvo[alvo]][universo])

        if not recuperados:
            continue
        recuperado = np.concatenate(recuperados)
        pedido = np.concatenate(pedidos)
        denominador = float(np.square(pedido).sum())
        inclinacao = float(pedido @ recuperado) / denominador if denominador > 0 else 0.0
        correlacao = float(np.corrcoef(pedido, recuperado)[0, 1])
        vies = float(np.mean(recuperado - pedido))
        rms = float(np.sqrt(np.mean(np.square(recuperado - pedido))))

        ok_inclinacao = INCLINACAO_MINIMA <= inclinacao <= INCLINACAO_MAXIMA
        ok_vies = abs(vies) <= VIES_MAXIMO
        marca = "ok  " if (ok_inclinacao and ok_vies) else "FALHA"
        reprovou = reprovou or not (ok_inclinacao and ok_vies)
        print(
            f"  [{marca}] {contexto}: inclinacao={inclinacao:+.4f} "
            f"correlacao={correlacao:+.4f} vies={vies:+.6f} rms={rms:.4f} "
            f"({len(recuperados)} alvos x {int(universo.sum())} genes)"
        )
        resultados[contexto] = {
            "inclinacao": inclinacao,
            "correlacao": correlacao,
            "vies": vies,
            "rms_do_residuo": rms,
            "alvos": float(len(recuperados)),
            "genes_do_universo": float(universo.sum()),
            "inclinacao_ok": float(ok_inclinacao),
            "vies_ok": float(ok_vies),
        }

    reprovou = reprovou or not procedencia_ok
    laudo = {
        "arquivo": str(caminho),
        "procedencia": {"aprovada": procedencia_ok, "nota": procedencia_nota},
        "delta": str(argumentos.delta),
        "limites": {
            "inclinacao": [INCLINACAO_MINIMA, INCLINACAO_MAXIMA],
            "vies_maximo": VIES_MAXIMO,
        },
        "por_contexto": resultados,
        "aprovado": not reprovou,
        "o_que_isto_nao_testa": (
            "Que o modelo esta certo. Um delta errado emitido fielmente passa "
            "aqui. Para isso existe o avaliador_proxy. E a REGRESSAO sozinha nao "
            "distingue dois deltas parecidos: entre v5 e v6, que diferem so na "
            "parcela da media, a inclinacao cruzada da 0,996. Quem pega isso e a "
            "checagem de procedencia por sha256, nao o limite da inclinacao."
        ),
    }
    destino = Path(argumentos.laudo) if argumentos.laudo else None
    if destino is not None:
        destino.parent.mkdir(parents=True, exist_ok=True)
        destino.write_text(json.dumps(laudo, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"laudo em {destino}")

    print("\n" + ("REPROVADO" if reprovou else "APROVADO: o arquivo carrega o delta pretendido"))
    return 1 if reprovou else 0


if __name__ == "__main__":
    raise SystemExit(main())
