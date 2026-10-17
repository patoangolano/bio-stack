"""Nivel 4 da escada (skill atavismo): lesao do ENAM COMPARTILHADA entre aves.

Le os dois alinhamentos multiplos do s31_aves_nt.sh (MAFFT L-INS-i e MUSCLE 5)
de 3 crocodilianos com o gene funcional (Alligator, Gavialis, Crocodylus) e 6
aves com reliquia (3 neognatas: Strigops, Anas, Taeniopygia; 3 paleognatas:
Dromaius, Struthio, Apteryx), todos na fita codificante.

LESAO de ave, polarizada pelos crocodilianos:
  delecao   corrida de gap na ave onde a maioria dos crocodilianos tem base
  insercao  corrida de base na ave onde a maioria dos crocodilianos tem gap
  frameshift  delecao ou insercao de tamanho nao multiplo de 3
  stop      codon de parada no quadro do jacare (bloco = aa 1-849 contiguos)

COMPARTILHADA: mesmo tipo, mesma posicao no jacare e mesmo tamanho em >= 2 aves.
ROBUSTA: aparece nos DOIS alinhadores. Indel em homopolimero (>= 4 bases iguais
encostadas) e marcado e nao conta: e onde montagem e alinhador erram.

O que o compartilhamento data (Palaeognathae x Neognathae e a divisao basal de
Neornithes): lesao em paleognata E neognata -> anterior ao ancestral comum das
aves vivas; so em neognatas -> anterior a Neognathae.

    uv run python scripts/atavismo/s31_nivel4.py
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

N4 = Path(__file__).resolve().parents[2] / "runs" / "atavismo" / "nivel4"
GENE = "ENAM"
CROCS = ["Amis", "Ggan", "Cpor"]
PALEO = {"Dnov", "Scam", "Arow"}
NEO = {"Shab", "Apla", "Tgut"}
STOPS = {"TAA", "TAG", "TGA"}
# Clados das 6 aves. Paleognatas: Struthio irmao das demais; Apteryx + Dromaius
# em Novaeratitae. Neognatas: Anas (Galloanserae) irmao de Strigops + Taeniopygia
# (Psittacopasserae).
CLADOS = [
    ("Neornithes", PALEO | NEO),
    ("Palaeognathae", PALEO),
    ("Novaeratitae (Apteryx+Dromaius)", {"Arow", "Dnov"}),
    ("Neognathae", NEO),
    ("Psittacopasserae (Strigops+Taeniopygia)", {"Shab", "Tgut"}),
]


def le_fasta(caminho: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    nome = ""
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        if linha.startswith(">"):
            nome = linha[1:].split()[0].replace(f"{GENE}_", "")
            seqs[nome] = []
        elif nome:
            seqs[nome].append(linha.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def homopolimero(seq: str, i: int, j: int) -> bool:
    """Ha corrida >= 4 da mesma base encostada no trecho [i, j) da sequencia sem gap?"""
    trecho = seq[max(0, i - 4): j + 4].replace("-", "")
    return any(trecho[k:k + 4] == trecho[k] * 4 for k in range(len(trecho) - 3))


def lesoes(msa: dict[str, str]) -> dict[str, set[tuple[str, int, int, bool]]]:
    """ave -> {(tipo, posicao no jacare, tamanho, em_homopolimero)}."""
    ref = msa["Amis"]
    ncol = len(ref)
    # posicao no jacare (0-based) de cada coluna; colunas de gap no jacare herdam a anterior
    pos_ref, p = [], -1
    for c in ref:
        if c != "-":
            p += 1
        pos_ref.append(p)
    croc_base = [sum(msa[cr][c] != "-" for cr in CROCS) >= 2 for c in range(ncol)]
    saida: dict[str, set[tuple[str, int, int, bool]]] = {}
    for ave, s in msa.items():
        if ave in CROCS:
            continue
        cob = [i for i, b in enumerate(s) if b != "-"]
        if not cob:
            continue
        ini, fim = cob[0], cob[-1]
        ev: set[tuple[str, int, int, bool]] = set()
        c = ini
        while c <= fim:
            if s[c] == "-" and croc_base[c]:  # delecao na ave
                d = c
                while d <= fim and not (s[d] != "-" and croc_base[d]) and (s[d] == "-" or not croc_base[d]):
                    d += 1
                n = sum(1 for k in range(c, d) if croc_base[k] and s[k] == "-")
                ev.add(("del", pos_ref[c], n, homopolimero(msa["Amis"], c, d)))
                c = d
            elif s[c] != "-" and not croc_base[c]:  # insercao na ave
                d = c
                while d <= fim and not croc_base[d]:
                    d += 1
                n = sum(1 for k in range(c, d) if s[k] != "-")
                if n:
                    ev.add(("ins", pos_ref[c], n, homopolimero(s, c, d)))
                c = d
            else:
                c += 1
        # stops no quadro do jacare: codon k = posicoes 3k..3k+2 do jacare
        col_de: dict[int, int] = {}
        for col, b in enumerate(ref):
            if b != "-":
                col_de[pos_ref[col]] = col
        for k in range(len(ref.replace("-", "")) // 3):
            cols = [col_de.get(3 * k + x) for x in range(3)]
            if None in cols or not (ini <= cols[0] and cols[2] <= fim):  # type: ignore[operator]
                continue
            codon = "".join(s[x] for x in cols if x is not None)
            if len(codon) == 3 and "-" not in codon and codon in STOPS:
                ev.add(("stop", 3 * k, 3, False))
        saida[ave] = ev
    return saida


def main() -> int:
    alinhadores = {"mafft": le_fasta(N4 / f"{GENE}_mafft.fna"), "muscle": le_fasta(N4 / f"{GENE}_muscle.fna")}
    por_alinhador = {}
    for nome, msa in alinhadores.items():
        ev = lesoes(msa)
        por_alinhador[nome] = ev
        tot = {a: sum(1 for e in v if e[0] == "stop" or e[2] % 3) for a, v in ev.items()}
        print(f"{nome:6s}: lesoes (frameshift + stop) por ave: "
              + ", ".join(f"{a} {n}" for a, n in sorted(tot.items())))

    # compartilhamento: chave = (tipo, posicao, tamanho); robusto = nos dois alinhadores
    def compartilhadas(ev: dict[str, set]) -> dict[tuple, set[str]]:
        quem: dict[tuple, set[str]] = defaultdict(set)
        for ave, eventos in ev.items():
            for tipo, pos, n, homo in eventos:
                if homo:
                    continue
                if tipo == "stop" or n % 3:
                    quem[(tipo, pos, n)].add(ave)
        return {k: v for k, v in quem.items() if len(v) >= 2}

    c1, c2 = compartilhadas(por_alinhador["mafft"]), compartilhadas(por_alinhador["muscle"])
    robustas = {k: c1[k] & c2[k] for k in c1.keys() & c2.keys() if len(c1[k] & c2[k]) >= 2}

    print(f"\nlesoes compartilhadas por >= 2 aves: mafft {len(c1)}, muscle {len(c2)}, "
          f"nos DOIS (robustas) {len(robustas)}")

    # Quem COBRE a posicao e nao tem a lesao. Datacao so vale se a lesao e
    # exatamente a intersecao das aves que cobrem com um clado da arvore.
    msa = alinhadores["mafft"]
    ref = msa["Amis"]
    col_de, p = {}, -1
    for col, b in enumerate(ref):
        if b != "-":
            p += 1
            col_de[p] = col

    def cobre(ave: str, pos: int, n: int) -> bool:
        s = msa[ave]
        cols = [col_de[x] for x in range(pos - 3, pos + n + 3) if x in col_de]
        return bool(cols) and all(s[c] != "-" for c in (cols[0], cols[-1]))

    def codon(nome: str, k: int) -> str:
        return "".join(msa[nome][col_de[3 * k + x]] for x in range(3))

    print(f"\n{'tipo':5s} {'codon':>6s} {'tam':>4s}  {'tem':20s} {'cobre e NAO tem':16s} "
          f"{'jacare':6s} {'stop nas aves':16s} {'CpG':4s} leitura")
    linhas = []
    for (tipo, pos, n), aves in sorted(robustas.items(), key=lambda x: x[0][1]):
        cobrem = {a for a in PALEO | NEO if a in msa and cobre(a, pos, n)}
        sem = cobrem - aves
        # Menor clado que contem as aves com a lesao (o MRCA delas); a lesao data
        # esse clado so se TODA ave que cobre a posicao dentro dele a tiver.
        menor = min((x for x in CLADOS if aves <= x[1]), key=lambda x: len(x[1]))
        clado = menor[0] if cobrem & menor[1] == aves else None
        k = pos // 3
        anc = codon("Amis", k) if tipo == "stop" else "-"
        stops_aves = ",".join(sorted({codon(a, k) for a in aves})) if tipo == "stop" else "-"
        cpg = "sim" if tipo == "stop" and anc == "CGA" else ""
        if clado is None:
            leitura = "INCONGRUENTE com a arvore (homoplasia ou artefato)"
        elif tipo == "stop" and (cpg or len(stops_aves.split(",")) > 1):
            leitura = f"{clado} — fraco ({'CpG' if cpg else 'stops diferentes'})"
        else:
            leitura = f"data {clado}"
        print(f"{tipo:5s} {k + 1:6d} {n:4d}  {','.join(sorted(aves)):20s} {','.join(sorted(sem)) or '-':16s} "
              f"{anc:6s} {stops_aves:16s} {cpg:4s} {leitura}")
        linhas.append({"gene": GENE, "tipo": tipo, "pos_jacare_nt": pos + 1, "codon": k + 1,
                       "tamanho": n, "aves": ",".join(sorted(aves)), "cobrem_sem": ",".join(sorted(sem)),
                       "codon_jacare": anc, "stops_aves": stops_aves, "cpg": cpg, "leitura": leitura})

    # controle: os crocodilianos tambem nao deveriam ter frameshift entre si
    for nome, msa in alinhadores.items():
        croc_fs = 0
        for cr in ("Ggan", "Cpor"):
            sub = {"Amis": msa["Amis"], cr: msa[cr], **{x: msa[x] for x in CROCS if x not in ("Amis", cr)}}
            ev = lesoes({k: v for k, v in sub.items()})
            croc_fs += sum(1 for e in ev.get(cr, set()) if (e[0] == "stop" or e[2] % 3))
        print(f"controle {nome}: lesoes nos crocodilianos com dente = {croc_fs} (esperado 0)")

    with (N4 / f"s31_{GENE}_compartilhadas.tsv").open("w", encoding="utf-8", newline="") as f:
        campos = list(linhas[0]) if linhas else ["gene"]
        w = csv.DictWriter(f, fieldnames=campos, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
