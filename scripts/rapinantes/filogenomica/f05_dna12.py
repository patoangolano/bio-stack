"""
f05 — matriz de nucleotideos retrotraduzida do alinhamento de AA, posicoes 1+2.

Para cada loco da execucao AA do pipeline:
  1. le o alinhamento de AA do MAFFT (02_alinhamento/<loco>.fas) e o log do
     ClipKIT (03_limpeza/<loco>.limpo.fa.log: uma linha por coluna, keep/trim);
  2. le o CDS de cada taxon (~/rapinantes/locos/cds_final/<loco>.fna) e confere
     que a traducao dele e identica a sequencia de AA sem gaps (se nao for, o
     loco para com erro: o codon estaria no lugar errado);
  3. troca cada residuo pelo codon correspondente e cada gap por '---';
  4. aplica a mesma mascara de colunas do ClipKIT (coluna de AA mantida -> os
     tres nucleotideos mantidos), entao AA e DNA12 cobrem os mesmos sitios;
  5. tira a 3a posicao de cada codon e grava <saida>/<loco>.fa (cabecalho = taxon).

A matriz sai alinhada e limpa: roda no pipeline com --alinhador nenhum
--limpeza nenhuma (limpar de novo removeria colunas diferentes das do AA).

    uv run python scripts/rapinantes/filogenomica/f05_dna12.py <run_id_AA> <dir_saida>
"""

from __future__ import annotations

import sys
from pathlib import Path

from Bio.Seq import Seq

RAIZ = Path(__file__).resolve().parents[3]
CDS = Path(r"\\wsl.localhost\Ubuntu-24.04\home\manec\rapinantes\locos\cds_final")


def ler_fasta(caminho: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    nome, partes = None, []
    with caminho.open(encoding="utf-8") as fh:
        for linha in fh:
            if linha.startswith(">"):
                if nome is not None:
                    seqs[nome] = "".join(partes)
                nome, partes = linha[1:].split()[0], []
            else:
                partes.append(linha.strip())
    if nome is not None:
        seqs[nome] = "".join(partes)
    return seqs


def mascara_clipkit(log: Path, n_col: int) -> list[bool]:
    manter = []
    for linha in log.open(encoding="utf-8"):
        c = linha.split()
        if len(c) >= 2:
            manter.append(c[1] == "keep")
    if len(manter) != n_col:
        raise SystemExit(f"{log}: {len(manter)} colunas no log, {n_col} no alinhamento")
    return manter


def main() -> None:
    run = RAIZ / "runs" / "filogenia" / sys.argv[1]
    saida = Path(sys.argv[2])
    saida.mkdir(parents=True, exist_ok=True)
    n_loc = n_sitios = 0
    for aln_f in sorted((run / "02_alinhamento").glob("*.fas")):
        loco = aln_f.stem
        aln = {k: v.upper() for k, v in ler_fasta(aln_f).items()}
        cds = ler_fasta(CDS / f"{loco}.fna")
        n_col = len(next(iter(aln.values())))
        manter = mascara_clipkit(run / "03_limpeza" / f"{loco}.limpo.fa.log", n_col)
        # O PREPARA_LOCO transforma '*' em nada; os CDS ja vem sem stop, entao o
        # residuo i do AA sem gap e o codon i do CDS.
        linhas = []
        for taxon, seq_aa in aln.items():
            nt = cds[taxon]
            if str(Seq(nt).translate(table="Standard")) != seq_aa.replace("-", ""):
                raise SystemExit(f"{loco}/{taxon}: traducao do CDS difere do AA alinhado")
            k, cods = 0, []
            for aa, m in zip(seq_aa, manter, strict=True):
                if aa == "-":
                    cod = "---"
                else:
                    cod = nt[3 * k: 3 * k + 3]
                    k += 1
                if m:
                    cods.append(cod[:2])
            linhas.append(f">{taxon}\n{''.join(cods)}\n")
        (saida / f"{loco}.fa").write_text("".join(linhas), encoding="utf-8", newline="\n")
        n_loc += 1
        n_sitios += 2 * sum(manter)
    print(f"{n_loc} locos DNA12 em {saida}; {n_sitios} sitios (posicoes 1+2)")


if __name__ == "__main__":
    main()
