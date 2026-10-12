"""Pontua uma edicao planejada contra a referencia com Evo 2, COM CONTROLES.

O delta de log-verossimilhanca sozinho nao e interpretavel: nao existe escala
absoluta, e o valor depende do comprimento, do conteudo GC e de quanto o
modelo conhece aquele clado. Este script so reporta um delta ao lado de dois
controles medidos na MESMA sequencia:

  - controle negativo: mutacoes sinonimas/neutras em numero igual, sorteadas
    fora do sitio de interesse. Da o ruido de fundo.
  - controle positivo: lesao que sabidamente quebra a sequencia (destruir o
    ATG inicial, ou introduzir parada prematura). Da a escala do "quebrado".

A edicao so e informativa se o delta dela sair FORA da faixa do controle
negativo. Se cair dentro, o resultado e "o modelo nao distingue", nao
"a edicao e segura".

LIMITE HONESTO: a validacao publicada de modelos de linguagem genomica
concentra-se em humano e procarioto. Em sauropsida nao-modelo isto e sinal
exploratorio. Reporte a comparacao, nunca o valor absoluto.

Uso (ASPAS OBRIGATORIAS: sem elas o shell le o '>' como redirecionamento e
cria um arquivo chamado 'T' em vez de passar a edicao):

    python 02_evo2_edicao.py alvo.fasta --edicao "150:A>T" --edicao "151:G>C"
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parents[2]
SAIDA = RAIZ / "runs" / "atavismo"
COMPLEMENTO = {"A": "T", "T": "A", "G": "C", "C": "G"}


@dataclass
class Edicao:
    """Uma substituicao pontual 1-based: posicao, base de referencia, base nova."""

    pos: int
    de: str
    para: str

    @classmethod
    def parse(cls, s: str) -> Edicao:
        # formato "150:A>T"
        try:
            pos_s, troca = s.split(":")
            de, para = troca.split(">")
            return cls(int(pos_s), de.strip().upper(), para.strip().upper())
        except ValueError as e:
            raise ValueError(f"edicao mal formada: {s!r} (esperado '150:A>T')") from e

    def aplica(self, seq: str) -> str:
        i = self.pos - 1
        if not 0 <= i < len(seq):
            raise ValueError(f"posicao {self.pos} fora da sequencia ({len(seq)} nt)")
        if seq[i] != self.de:
            raise ValueError(
                f"posicao {self.pos}: referencia tem {seq[i]}, edicao diz {self.de}. "
                "Coordenada errada ou fita errada."
            )
        return seq[:i] + self.para + seq[i + 1 :]


@dataclass
class Resultado:
    rotulo: str
    seq: str
    mean_ll: float = 0.0
    perplexidade: float = 0.0
    n_sitios: int = 0
    detalhe: dict = field(default_factory=dict)


def le_fasta(caminho: Path) -> tuple[str, str]:
    linhas = caminho.read_text(encoding="utf-8").strip().splitlines()
    cab = linhas[0].lstrip(">") if linhas and linhas[0].startswith(">") else caminho.name
    seq = "".join(ln.strip() for ln in linhas if not ln.startswith(">")).upper()
    seq = seq.replace("U", "T")
    invalidas = set(seq) - set("ACGTN")
    if invalidas:
        raise ValueError(f"bases invalidas na entrada: {sorted(invalidas)}")
    return cab, seq


def controle_negativo(seq: str, n: int, evitar: set[int], rng: np.random.Generator) -> str:
    """n transversoes sorteadas fora das posicoes de interesse."""
    livres = [i for i in range(len(seq)) if i not in evitar and seq[i] in COMPLEMENTO]
    if len(livres) < n:
        raise ValueError("sequencia curta demais para o controle negativo")
    escolhidas = rng.choice(livres, size=n, replace=False)
    lista = list(seq)
    for i in escolhidas:
        atual = lista[int(i)]
        opcoes = [b for b in "ACGT" if b != atual]
        lista[int(i)] = str(rng.choice(opcoes))
    return "".join(lista)


def controle_positivo(seq: str) -> tuple[str, str]:
    """Lesao que sabidamente quebra: destroi o primeiro ATG, ou insere parada."""
    i = seq.find("ATG")
    if i >= 0:
        return seq[:i] + "ACG" + seq[i + 3 :], f"ATG->ACG na posicao {i + 1}"
    # sem ATG: introduz parada no meio, em quadro a partir do inicio
    meio = (len(seq) // 6) * 3
    return seq[:meio] + "TAA" + seq[meio + 3 :], f"parada TAA inserida na posicao {meio + 1}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("fasta", type=Path, help="sequencia de referencia")
    ap.add_argument(
        "--edicao",
        action="append",
        default=[],
        metavar="POS:DE>PARA",
        help="substituicao 1-based, repetivel (ex: 150:A>T)",
    )
    ap.add_argument(
        "--replicas-controle",
        type=int,
        default=5,
        help="quantos sorteios do controle negativo (padrao 5)",
    )
    ap.add_argument(
        "--seed", type=int, default=0, help="semente do gerador; NUNCA usar np.random.seed global"
    )
    args = ap.parse_args(argv)

    if not args.edicao:
        ap.error("informe ao menos uma --edicao")

    cab, ref = le_fasta(args.fasta)
    edicoes = [Edicao.parse(s) for s in args.edicao]

    variante = ref
    for e in edicoes:
        variante = e.aplica(variante)
    alvos = {e.pos - 1 for e in edicoes}

    rng = np.random.default_rng(args.seed)
    dir_ent = SAIDA / "evo2_entradas"
    dir_ent.mkdir(parents=True, exist_ok=True)

    conjunto: list[tuple[str, str, str]] = [
        ("referencia", ref, "sequencia original"),
        ("variante", variante, "; ".join(f"{e.pos}:{e.de}>{e.para}" for e in edicoes)),
    ]
    for r in range(args.replicas_controle):
        conjunto.append(
            (
                f"ctrl_neg_{r}",
                controle_negativo(ref, len(edicoes), alvos, rng),
                f"{len(edicoes)} substituicoes sorteadas fora do alvo, seed {args.seed}+{r}",
            )
        )
    seq_pos, desc_pos = controle_positivo(ref)
    conjunto.append(("ctrl_pos", seq_pos, desc_pos))

    manifesto = []
    for rotulo, seq, desc in conjunto:
        arq = dir_ent / f"{rotulo}.fasta"
        arq.write_text(f">{rotulo} {desc}\n{seq}\n", encoding="utf-8")
        manifesto.append(
            {
                "rotulo": rotulo,
                "arquivo": str(arq),
                "descricao": desc,
                "comprimento": len(seq),
                "n_diferencas_vs_ref": sum(1 for a, b in zip(ref, seq, strict=True) if a != b),
            }
        )

    (SAIDA / "evo2_manifesto.json").write_text(
        json.dumps(
            {
                "referencia_cabecalho": cab,
                "fonte": str(args.fasta),
                "seed": args.seed,
                "entradas": manifesto,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print(f"Referencia: {cab}  ({len(ref)} nt)")
    print(f"Edicoes: {len(edicoes)}")
    print(f"\n{len(conjunto)} sequencias escritas em {dir_ent}\n")
    for m in manifesto:
        print(f"  {m['rotulo']:<14s} {m['n_diferencas_vs_ref']:>3d} dif  {m['descricao']}")
    print(f"\nManifesto: {SAIDA / 'evo2_manifesto.json'}")
    print("\nProximo passo: pontuar cada FASTA com mcp__nvidia-bio__evo2_likelihood")
    print("e rodar 03_evo2_consolida.py com os mean_log_likelihood.")
    print("\nLEITURA: a edicao so e informativa se o delta dela sair FORA da faixa")
    print("dos controles negativos. Dentro da faixa significa 'o modelo nao")
    print("distingue' -- nao significa 'a edicao e segura'.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
