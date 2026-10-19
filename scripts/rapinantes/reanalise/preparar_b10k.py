"""
Prepara as arvores de genes do B10K (Stiller et al. 2024, Nature 629:851) para o ASTER.

    python3 -I preparar_b10k.py <dir brutos> <dir saida> <folhas_63K.txt>

Por tipo de dado (intergenico 63K, exons C12, introns, UCEs):
  - tira a primeira coluna (nome do loco) dos arquivos *.named.*;
  - troca o rotulo de suporte ")/0.999" por ")0.999": o wASTRAL le so o
    primeiro numero do rotulo e aborta se ele comecar com "/" (docs/FILOGENIA.md);
  - traduz os nomes de folha para o nome de especie usado na arvore de especies
    63K (exons: codigo de 6 letras GENESP; UCEs: uce-NNNN_<amostra> -> mapa
    uce_name_map.txt); folha sem traducao unica e contada e removida da arvore
    nao, o script para com erro.
Saida: <tipo>.alrt.tre (uma arvore por linha, suporte aLRT 0-1) e
<tipo>.mapa.json (traducoes e contagens).
Somente biblioteca padrao: roda no python3 do WSL sem venv.
"""

from __future__ import annotations

import collections
import gzip
import json
import re
import sys
from pathlib import Path

FOLHA = re.compile(r"(?<=[(,])([^(),:;]+)(?=[:),])")
SUPORTE_BARRA = re.compile(r"\)/")


def ler_folhas(caminho: Path) -> list[str]:
    return [x.strip() for x in caminho.read_text().splitlines() if x.strip()]


def mapa_exon(especies: list[str]) -> dict[str, str]:
    cod: dict[str, list[str]] = collections.defaultdict(list)
    for e in especies:
        g, s = e.split("_", 1)
        cod[(g[:3] + s[:3]).upper()].append(e)
    m = {k: v[0] for k, v in cod.items() if len(v) == 1}
    # Codigos que nao seguem GENESP do nome usado na arvore 63K: os cormoroes
    # estao com o genero atual (Urile, Nannopterum) e os dois Pomato- (colisao
    # POMRUF) ganharam a 3a letra distinta. A troca entre os dois Pomato-, se
    # estiver invertida, nao altera nenhum ramo fora de Passeriformes: os dois
    # caem sempre no mesmo grupo dos quatro em volta de um ramo nao passeriforme.
    m.update({
        "URIPEL": "Phalacrocorax_pelagicus",
        "NANBRA": "Phalacrocorax_brasilianus",
        "NANHAR": "Phalacrocorax_harrisi",
        "NANAUR": "Phalacrocorax_auritus",
        "POSRUF": "Pomatostomus_ruficeps",
        "PORRUF": "Pomatorhinus_ruficollis",
    })
    return {k: v for k, v in m.items() if v in especies}


def mapa_uce(caminho: Path) -> dict[str, str]:
    m = {}
    for linha in caminho.read_text().splitlines():
        p = linha.split()
        if len(p) == 2:
            m[p[0]] = p[1]
    return m


def processar(entrada: Path, saida: Path, traduz, nome: str) -> dict:
    cont = collections.Counter()
    sem_traducao = collections.Counter()
    n = 0
    with gzip.open(entrada, "rt") as fin, open(saida, "w") as fout:
        for linha in fin:
            linha = linha.strip()
            if not linha:
                continue
            arv = linha.split("\t")[-1].strip()
            arv = SUPORTE_BARRA.sub(")", arv)

            def sub(m: re.Match) -> str:
                rot = m.group(1)
                novo = traduz(rot)
                if novo is None:
                    sem_traducao[rot] += 1
                    return rot
                cont[novo] += 1
                return novo

            arv = FOLHA.sub(sub, arv)
            fout.write(arv + "\n")
            n += 1
    if sem_traducao:
        raise SystemExit(f"[{nome}] {len(sem_traducao)} rotulos sem traducao, ex.: {list(sem_traducao)[:10]}")
    return {"arvores": n, "ocorrencias_por_especie": dict(sorted(cont.items()))}


def main() -> None:
    brutos, saida, folhas = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    saida.mkdir(parents=True, exist_ok=True)
    especies = ler_folhas(folhas)
    conj = set(especies)
    mex = mapa_exon(especies)
    muce = mapa_uce(brutos / "01_alignments_and_gene_trees_uces_uce_name_map.txt")
    uce_pref = re.compile(r"^uce-\d+_")

    def t_ident(r: str):
        return r if r in conj else None

    def t_exon(r: str):
        return mex.get(r)

    def t_uce(r: str):
        e = muce.get(uce_pref.sub("", r))
        return e if e in conj else None

    tarefas = {
        "intergenico63K": ("01_alignments_and_gene_trees_intergenic_regions_63430.named.gene.trees.gz", t_ident),
        "intergenico63K_c95": ("01_alignments_and_gene_trees_intergenic_regions_63430.aLRT-0.95-collapsed.gene.trees.gz", t_ident),
        "exon_c12": ("01_alignments_and_gene_trees_exons_named.trees.gz", t_exon),
        "intron": ("01_alignments_and_gene_trees_introns_named.trees.gz", t_ident),
        "uce": ("01_alignments_and_gene_trees_uces_named.trees.gz", t_uce),
    }
    alvo = sys.argv[4:] or list(tarefas)
    for nome in alvo:
        arq, f = tarefas[nome]
        info = processar(brutos / arq, saida / f"{nome}.alrt.tre", f, nome)
        info["fonte"] = arq
        (saida / f"{nome}.mapa.json").write_text(json.dumps(info, indent=1))
        print(nome, info["arvores"], "arvores;", len(info["ocorrencias_por_especie"]), "especies")


if __name__ == "__main__":
    main()
