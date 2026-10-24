#!/usr/bin/env python
"""Validador NUMERICO da submissao do VCC 2026: os numeros emitidos sao os pretendidos?

O validador de contrato (`validar_submissao.py`) confere estrutura, rotulo e eixo
e aprova 28/28 um arquivo de delta zero. Este fecha o buraco: recomputa o
pseudobulk de cada (contexto, alvo) no espaco do placar (`bulk_lognorm`), subtrai
o pseudobulk dos controles do mesmo contexto, e regride esse delta EMITIDO contra
o delta PRETENDIDO do arquivo de delta. Inclinacao em [0,85; 1,15] aprova.

E a inclinacao que prova fidelidade, nao a correlacao. Com 400 celulas por
grupo, o ruido de amostragem por gene (rms ~0,038 na v4) supera o |delta| tipico
(~0,0018), entao a correlacao de um arquivo correto fica em ~0,4 -- mas ruido
independente do x nao enviesa a inclinacao. Medido no registro: v4 correta
+0,998; v1, com 143 dos 300 alvos emitidos como controle, +0,809 e reprovada;
stub de delta zero -0,001. Pega erro PARCIAL, nao so total.

A inclinacao e CEGA a deslocamento por gene descorrelacionado do delta -- e e
exatamente o que um gerador produz se partir da linha de base errada (a armadilha
registrada: o `ctrl_means` do PIE e media de log1p por celula, nao log1p do
pseudobulk). Medido no teste: linha de base errada por N(0; 0,3) em cada gene da
inclinacao 1,005 e passaria. Por isso o segundo criterio: o residuo medio POR
GENE ao longo das perturbacoes, em rms sobre os genes, tem de ficar abaixo de
`LIMITE_DESLOCAMENTO`. Com ruido de amostragem independente ele encolhe com
raiz(P): na escala da v4, ~0,038/raiz(300) = 0,002. Um gerador que reusa as
mesmas celulas-molde em toda perturbacao compartilha o ruido e da ~0,04. Linha
de base errada da >> 0,1. O limite de 0,1 e PROVISORIO: calibrar na primeira
submissao real e registrar o valor medido. A `razao_deslocamento` (= rms por
gene x raiz(P) / rms do residuo) e diagnostico: ~1 com ruido independente, ~raiz(P)
com molde compartilhado.

O gene-alvo de cada perturbacao fica fora da regressao, como fica fora de todas
as metricas do placar.

O `--delta` MUDA com a versao: conferir um arquivo bom contra o delta de outro k
reprova, e esse e o comportamento certo.

Uso:
    uv run python scripts/vcc/validar_numerico.py submissao.h5ad --delta delta.npz
Sai com 0 se aprovado e 1 se reprovado.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from formato_delta import ler_delta
from pseudobulk import LINHAS_POR_BLOCO, bulk_lognorm, somas_por_grupo

RAIZ = Path(__file__).resolve().parents[2]
CONTROLES_PADRAO = RAIZ / "data" / "raw" / "vcc" / "controls-final"
FAIXA = (0.85, 1.15)
LIMITE_DESLOCAMENTO = 0.1  # provisorio, ver o docstring


@dataclass(frozen=True)
class Ajuste:
    contexto: str
    n_pontos: int
    inclinacao: float
    intercepto: float
    correlacao: float
    rms_ruido: float
    rms_delta: float
    deslocamento_por_gene: float
    razao_deslocamento: float
    passou: bool


@dataclass
class MedidaNumerica:
    arquivo: str
    delta: str
    ajustes: list[Ajuste] = field(default_factory=list)

    @property
    def veredito(self) -> str:
        return "aprovado" if self.ajustes and all(a.passou for a in self.ajustes) else "reprovado"

    def para_dict(self) -> dict:
        return {
            "arquivo": self.arquivo,
            "delta": self.delta,
            "faixa_de_inclinacao": list(FAIXA),
            "limite_deslocamento_por_gene": LIMITE_DESLOCAMENTO,
            "ajustes": [asdict(a) for a in self.ajustes],
            "veredito": self.veredito,
        }


def _ajustar(contexto: str, pretendido: np.ndarray, emitido: np.ndarray, usar: np.ndarray) -> Ajuste:
    """Regressao do emitido no pretendido sobre `usar` ([alvos, genes], gene-alvo fora)."""
    x, y = pretendido[usar], emitido[usar]
    xm, ym = float(x.mean()), float(y.mean())
    vx = float(((x - xm) ** 2).mean())
    vy = float(((y - ym) ** 2).mean())
    cxy = float(((x - xm) * (y - ym)).mean())
    inclinacao = cxy / vx if vx > 0 else math.nan
    intercepto = ym - inclinacao * xm
    correlacao = cxy / math.sqrt(vx * vy) if vx > 0 and vy > 0 else math.nan
    residuo = emitido - (intercepto + inclinacao * pretendido)
    rms_ruido = float(np.sqrt((residuo[usar] ** 2).mean()))

    # Residuo medio por gene ao longo das perturbacoes: o que a inclinacao nao ve.
    contagem = usar.sum(axis=0)
    tem = contagem > 0
    por_gene = np.where(usar, residuo, 0.0).sum(axis=0)[tem] / contagem[tem]
    deslocamento = float(np.sqrt((por_gene**2).mean())) if por_gene.size else math.nan
    razao = deslocamento * math.sqrt(pretendido.shape[0]) / rms_ruido if rms_ruido > 0 else math.nan
    return Ajuste(
        contexto=contexto,
        n_pontos=int(x.size),
        inclinacao=inclinacao,
        intercepto=intercepto,
        correlacao=correlacao,
        rms_ruido=rms_ruido,
        rms_delta=float(np.sqrt((x**2).mean())),
        deslocamento_por_gene=deslocamento,
        razao_deslocamento=razao,
        # NaN em qualquer um dos dois reprova.
        passou=FAIXA[0] <= inclinacao <= FAIXA[1] and deslocamento <= LIMITE_DESLOCAMENTO,
    )


def medir(
    submissao: Path,
    delta: Path,
    controles: Path = CONTROLES_PADRAO,
    linhas_por_bloco: int = LINHAS_POR_BLOCO,
    col_alvo: str = "target_gene",
    col_contexto: str = "context",
) -> MedidaNumerica:
    pretendido = ler_delta(delta)
    chaves, somas, genes = somas_por_grupo(submissao, [col_contexto, col_alvo], linhas_por_bloco)
    if genes != pretendido.genes:
        raise ValueError("eixo de genes da submissao difere do arquivo de delta")
    linha_de = {k: i for i, k in enumerate(chaves)}
    coluna_de = {g: j for j, g in enumerate(genes)}

    medida = MedidaNumerica(str(submissao), str(delta))
    for c in pretendido.contextos:
        chaves_c, somas_c, genes_c = somas_por_grupo(
            Path(controles) / f"context_{c}.h5ad", [col_contexto], linhas_por_bloco
        )
        if genes_c != genes:
            raise ValueError(f"eixo de genes dos controles de {c} difere da submissao")
        b_ctrl = bulk_lognorm(somas_c[[chaves_c.index((c,))]])[0]

        faltam = [a for a in pretendido.alvos if (c, a) not in linha_de]
        if faltam:
            raise ValueError(f"{c}: {len(faltam)} alvo(s) ausentes na submissao, ex. {faltam[:5]}")
        emitido = bulk_lognorm(somas[[linha_de[(c, a)] for a in pretendido.alvos]]) - b_ctrl

        fora_do_alvo = np.ones(emitido.shape, dtype=bool)
        for i, a in enumerate(pretendido.alvos):
            if a in coluna_de:
                fora_do_alvo[i, coluna_de[a]] = False
        x = pretendido.por_contexto[c].astype(np.float64)
        medida.ajustes.append(_ajustar(c, x, emitido, fora_do_alvo))
    return medida


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("arquivo", type=Path, help="o .h5ad da submissao")
    ap.add_argument("--delta", type=Path, required=True, help="o arquivo de delta DESTA versao")
    ap.add_argument("--controles", type=Path, default=CONTROLES_PADRAO, help="pasta com context_<c>.h5ad")
    ap.add_argument("--json", type=Path, default=None, help="grava a medida completa neste arquivo")
    ap.add_argument("--linhas-por-bloco", type=int, default=LINHAS_POR_BLOCO)
    args = ap.parse_args(argv)

    medida = medir(args.arquivo, args.delta, args.controles, args.linhas_por_bloco)
    for a in medida.ajustes:
        print(
            f"{'ok   ' if a.passou else 'FALHA'}  {a.contexto}: inclinacao {a.inclinacao:+.3f} "
            f"(faixa {FAIXA[0]}-{FAIXA[1]}), intercepto {a.intercepto:+.4f}, corr {a.correlacao:.3f}, "
            f"ruido rms {a.rms_ruido:.4f}, |delta| rms {a.rms_delta:.4f}, n {a.n_pontos:,}\n"
            f"       deslocamento por gene {a.deslocamento_por_gene:.4f} (limite {LIMITE_DESLOCAMENTO}), "
            f"razao {a.razao_deslocamento:.2f} (~1 ruido independente, ~raiz(P) molde compartilhado)"
        )
    print(f"\n{medida.veredito.upper()} contra {args.delta.name}")
    if args.json is not None:
        args.json.write_text(json.dumps(medida.para_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if medida.veredito == "aprovado" else 1


if __name__ == "__main__":
    sys.exit(main())
