#!/usr/bin/env python
"""Gera a submissao do VCC 2026 a partir de um arquivo de delta.

    uv run python scripts/vcc/gerar_submissao.py delta.npz --saida submissao.h5ad

Um contexto por vez, um grupo por vez, nada materializado. Para cada contexto:
le os momentos do pool de controle, mede a VIABILIDADE do delta antes de emitir
uma celula, e entao emite os 300 grupos de 400 celulas direto para o `.h5ad`.

Grava, ao lado da saida, `<saida>.efetivo.npz`: o delta EFETIVO, que e a predicao
que a submissao realmente faz. Essa distincao e o ponto. O arquivo de delta pode
pedir o impossivel -- contagem negativa -- e o clamp na fronteira factivel muda a
predicao. Validar a submissao contra o delta PEDIDO mede a distancia entre o que
o preditor quis e o que e fisicamente emissivel; contra o delta EFETIVO mede se o
gerador emitiu o que prometeu. As duas medidas sao uteis e sao diferentes, e e o
delta efetivo que vai ao `validar_numerico.py` como portao.

Tres portoes, todos antes de escrever:

1. **Disco.** 2,0e9 nao nulos sao 15 GB sem compressao; exige 25 GB livres.
2. **Viabilidade.** Se a inclinacao analitica contra o delta cru cai abaixo de
   `validar_numerico.FAIXA[0]`, o arquivo de delta esta pedindo o impossivel em
   volume bastante para reprovar a submissao -- e o conserto NAO e do gerador.
   E uma linha no preditor: gravar `delta >= -b_ctrl[g]`. Recusa, a menos de
   `--aceitar-delta-infactivel`, e o laudo registra a decisao.
3. **Teto de nnz.** 2,1e9 e 97,8% de 2^31-1. Se cruzar, a acao documentada e
   restringir o ramo `pi == 0` aos genes do gate.

A execucao de producao vai por `bootstrap/tarefa_longa.ps1`, nao por job do Bash
tool: sao cerca de 11 min de calculo e o job morre de SIGHUP ao fechar a sessao.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

# `-I` tira o diretorio do script do sys.path, e os cinco imports abaixo sao
# todos de modulos irmaos. O mesmo idioma do `varredura_k.py` e do
# `pool_media.py`, e aqui ele faltava: sem isso o script morre em
# `ModuleNotFoundError: No module named 'conversao'` justamente sob `-I`, que e
# o modo que o projeto exige para script que le dado baixado.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from conversao import PISO_LOG2, TETO_LOG2, Viabilidade, viabilidade_do_delta  # noqa: E402
from escritor import EscritorCSR, conferir_teto_por_celula  # noqa: E402
from formato_delta import Delta, ler_delta, salvar_delta  # noqa: E402
from kernel import N_CELULAS, GrupoEmitido, emitir_grupo  # noqa: E402
from momentos import carregar_pool, momentos_do_controle  # noqa: E402

RAIZ = Path(__file__).resolve().parents[2]
CONTROLES_PADRAO = RAIZ / "data" / "raw" / "vcc" / "controls-final"
TETO_DE_NNZ: int = 2_100_000_000
DISCO_MINIMO_GB: float = 25.0
PISO_DE_INCLINACAO: float = 0.85  # == validar_numerico.FAIXA[0]


@dataclass
class Relatorio:
    saida: str
    delta: str
    delta_efetivo: str
    semente_base: int
    piso_log2: float
    teto_log2: float
    n_celulas_por_grupo: int
    viabilidade: list[Viabilidade] = field(default_factory=list)
    grupos: list[GrupoEmitido] = field(default_factory=list)
    nnz_total: int = 0
    maior_soma_por_celula: float = 0.0
    celulas_acima_do_teto: int = 0
    n_celulas: int = 0
    delta_infactivel_aceito: bool = False

    def para_dict(self) -> dict:
        return {
            "saida": self.saida,
            "delta": self.delta,
            "delta_efetivo": self.delta_efetivo,
            "semente_base": self.semente_base,
            "piso_log2": self.piso_log2,
            "teto_log2": self.teto_log2,
            "n_celulas_por_grupo": self.n_celulas_por_grupo,
            "delta_infactivel_aceito": self.delta_infactivel_aceito,
            "nnz_total": self.nnz_total,
            "maior_soma_por_celula": self.maior_soma_por_celula,
            "celulas_acima_do_teto": self.celulas_acima_do_teto,
            "n_celulas": self.n_celulas,
            "viabilidade": [asdict(v) for v in self.viabilidade],
            "grupos": [
                {
                    "contexto": g.contexto,
                    "alvo": g.alvo,
                    "semente": g.semente,
                    "nnz": g.nnz,
                    "soma_maxima_por_celula": g.soma_maxima_por_celula,
                    "soma_minima_por_celula": g.soma_minima_por_celula,
                    "linhas_nulas": g.linhas_nulas,
                    "n_genes_criados": g.n_genes_criados,
                    "zerofrac_no_gate": g.zerofrac_no_gate,
                    "residuo_bulk_rms": g.residuo_bulk_rms,
                    "residuo_bulk_rms_no_gate": g.residuo_bulk_rms_no_gate,
                    "log2fc_do_alvo": g.log2fc_do_alvo,
                    "lambda_de_fechamento": g.alvo_de_composicao.lambda_de_fechamento,
                    "n_no_piso_factivel": g.alvo_de_composicao.n_no_piso_factivel,
                    "n_no_piso_do_gate": g.alvo_de_composicao.n_no_piso_do_gate,
                    "n_no_teto_do_gate": g.alvo_de_composicao.n_no_teto_do_gate,
                    "log2fc_min_no_gate": g.alvo_de_composicao.log2fc_min_no_gate,
                    "log2fc_max_no_gate": g.alvo_de_composicao.log2fc_max_no_gate,
                }
                for g in self.grupos
            ],
        }


def _conferir_disco(saida: Path, minimo_gb: float) -> None:
    pasta = saida.parent if saida.parent.exists() else Path.cwd()
    livre_gb = shutil.disk_usage(pasta).free / 1024**3
    if livre_gb < minimo_gb:
        raise ValueError(
            f"{livre_gb:.1f} GB livres em {pasta}, abaixo do minimo de {minimo_gb} GB "
            f"-- a submissao completa sao ~15 GB sem compressao"
        )


def gerar_submissao(
    delta: Path,
    controles: Path = CONTROLES_PADRAO,
    saida: Path = Path("submissao.h5ad"),
    semente_base: int = 0,
    piso_log2: float = PISO_LOG2,
    teto_log2: float = TETO_LOG2,
    contextos: Sequence[str] | None = None,
    n_celulas: int = N_CELULAS,
    teto_de_nnz: int = TETO_DE_NNZ,
    aceitar_delta_infactivel: bool = False,
    disco_minimo_gb: float = DISCO_MINIMO_GB,
    relatorio: Path | None = None,
) -> Relatorio:
    delta = Path(delta)
    controles = Path(controles)
    saida = Path(saida)
    pretendido = ler_delta(delta)
    quais = list(contextos) if contextos is not None else list(pretendido.contextos)
    faltam = [c for c in quais if c not in pretendido.por_contexto]
    if faltam:
        raise ValueError(f"contexto(s) {faltam} ausentes em {delta}")
    _conferir_disco(saida, disco_minimo_gb)

    caminho_efetivo = saida.with_suffix(".efetivo.npz")
    laudo = Relatorio(
        saida=str(saida),
        delta=str(delta),
        delta_efetivo=str(caminho_efetivo),
        semente_base=semente_base,
        piso_log2=piso_log2,
        teto_log2=teto_log2,
        n_celulas_por_grupo=n_celulas,
        delta_infactivel_aceito=aceitar_delta_infactivel,
    )

    # Portao 1: viabilidade de TODOS os contextos antes de emitir uma celula.
    # Custa segundos e e a diferenca entre descobrir agora e descobrir depois de
    # 11 minutos de calculo e 15 GB de escrita.
    momentos_por_contexto = {}
    for c in quais:
        momentos_por_contexto[c] = momentos_do_controle(controles / f"context_{c}.h5ad")
        v = viabilidade_do_delta(momentos_por_contexto[c], pretendido.por_contexto[c])
        laudo.viabilidade.append(v)
    # `not (x >= piso)` e nao `x < piso`: um delta degenerado -- constante em
    # todos os pares -- da variancia zero e inclinacao NaN, e `NaN < 0,85` e
    # `False`, logo o portao passaria o pior arquivo possivel em silencio.
    ruins = [
        v
        for v in laudo.viabilidade
        if not np.isfinite(v.inclinacao_contra_delta_cru)
        or not v.inclinacao_contra_delta_cru >= PISO_DE_INCLINACAO
    ]
    if ruins and not aceitar_delta_infactivel:
        linhas = "; ".join(
            f"{v.contexto}: inclinacao {v.inclinacao_contra_delta_cru:+.4f}, "
            f"{v.fracao_dos_nao_nulos_infactiveis:.1%} das predicoes infactiveis"
            for v in ruins
        )
        raise ValueError(
            f"arquivo de delta infactivel -- {linhas}. A inclinacao contra o delta cru fica "
            f"abaixo de {PISO_DE_INCLINACAO}, e `validar_numerico.py` reprovaria. O conserto "
            f"e no PREDITOR, nao aqui: gravar `delta >= -b_ctrl[g]` leva a inclinacao a 1,0. "
            f"Para emitir assim mesmo, passe --aceitar-delta-infactivel."
        )

    efetivo: dict[str, np.ndarray] = {}
    escritor = EscritorCSR(saida, pretendido.genes)
    try:
        for c in quais:
            momentos = momentos_por_contexto[c]
            if momentos.genes != pretendido.genes:
                raise ValueError(f"{c}: eixo de genes do controle difere do arquivo de delta")
            pool = carregar_pool(controles / f"context_{c}.h5ad")
            matriz_efetiva = np.zeros(
                (len(pretendido.alvos), len(pretendido.genes)), dtype=np.float32
            )
            linha_de_delta = pretendido.por_contexto[c]
            for i, alvo in enumerate(pretendido.alvos):
                g = emitir_grupo(
                    momentos,
                    pool,
                    np.asarray(linha_de_delta[i], dtype=np.float64),
                    alvo,
                    semente_base,
                    piso_log2,
                    teto_log2,
                    n_celulas,
                )
                laudo.nnz_total += g.nnz
                if laudo.nnz_total > teto_de_nnz:
                    raise ValueError(
                        f"nnz acumulado {laudo.nnz_total:,} acima do teto de {teto_de_nnz:,} "
                        f"em {c}/{alvo}. Acao documentada: restringir o ramo `pi == 0` aos "
                        f"genes do gate."
                    )
                escritor.acrescentar(g.matriz, c, alvo)
                matriz_efetiva[i] = g.alvo_de_composicao.delta_efetivo
                laudo.grupos.append(g)
            efetivo[c] = matriz_efetiva
            del pool
        laudo.n_celulas = escritor.fechar()  # devolve nnz; o n_obs vem logo abaixo
    except BaseException:
        escritor.fechar()
        raise

    laudo.n_celulas = len(quais) * len(pretendido.alvos) * n_celulas
    laudo.maior_soma_por_celula, laudo.celulas_acima_do_teto = conferir_teto_por_celula(saida)
    salvar_delta(
        caminho_efetivo,
        Delta(pretendido.genes, pretendido.alvos, quais, efetivo),
    )
    if relatorio is not None:
        Path(relatorio).write_text(
            json.dumps(laudo.para_dict(), ensure_ascii=False, indent=2), encoding="utf-8"
        )
    return laudo


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("delta", type=Path, help="o .npz do formato_delta")
    ap.add_argument("--controles", type=Path, default=CONTROLES_PADRAO)
    ap.add_argument("--saida", type=Path, default=Path("submissao.h5ad"))
    ap.add_argument("--semente-base", type=int, default=0)
    ap.add_argument("--piso-log2", type=float, default=PISO_LOG2)
    ap.add_argument("--teto-log2", type=float, default=TETO_LOG2)
    ap.add_argument("--contextos", nargs="*", default=None)
    ap.add_argument("--n-celulas", type=int, default=N_CELULAS)
    ap.add_argument("--teto-de-nnz", type=int, default=TETO_DE_NNZ)
    ap.add_argument("--aceitar-delta-infactivel", action="store_true")
    ap.add_argument("--disco-minimo-gb", type=float, default=DISCO_MINIMO_GB)
    ap.add_argument("--relatorio", type=Path, default=None)
    a = ap.parse_args(argv)

    laudo = gerar_submissao(
        delta=a.delta,
        controles=a.controles,
        saida=a.saida,
        semente_base=a.semente_base,
        piso_log2=a.piso_log2,
        teto_log2=a.teto_log2,
        contextos=a.contextos,
        n_celulas=a.n_celulas,
        teto_de_nnz=a.teto_de_nnz,
        aceitar_delta_infactivel=a.aceitar_delta_infactivel,
        disco_minimo_gb=a.disco_minimo_gb,
        relatorio=a.relatorio,
    )
    for v in laudo.viabilidade:
        print(
            f"{v.contexto}: inclinacao contra o delta cru {v.inclinacao_contra_delta_cru:+.4f}, "
            f"{v.fracao_dos_nao_nulos_infactiveis:.2%} das predicoes infactiveis, "
            f"{v.genes_atingidos:,} genes atingidos"
        )
    residuos = [g.residuo_bulk_rms for g in laudo.grupos]
    print(
        f"\n{len(laudo.grupos):,} grupos, {laudo.n_celulas:,} celulas, "
        f"{laudo.nnz_total:,} nao nulos ({laudo.nnz_total / 4.75e9:.1%} do teto do prep)\n"
        f"maior soma por celula {laudo.maior_soma_por_celula:,.0f} "
        f"({laudo.celulas_acima_do_teto} acima do teto)\n"
        f"residuo de pseudobulk rms mediano {float(np.median(residuos)):.4f}\n"
        f"delta efetivo em {laudo.delta_efetivo}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
