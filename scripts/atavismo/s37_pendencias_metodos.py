# /// script
# requires-python = ">=3.11"
# dependencies = ["requests>=2.31", "pandas>=2.0", "xlrd>=2.0.1", "biopython>=1.80"]
# ///
"""Fecha as pendencias 2-6 de docs/metodos.md com dado, nao com memoria.

2  acessos das proteinas do Anolis usadas como 2a consulta no s33 — o
   cabecalho do consultas.faa so guardou o simbolo. Refaz a mesma busca do
   s33 (s27.gene_pos + s27.tabela_features) e exige que a sequencia baixada
   seja IDENTICA a que entrou no consultas.faa.
3  ENAM (QOX06214.1) e AMTN (KYO35092.1) do jacare sao GenBank: compara com
   o XP da anotacao RefSeq da GCF_030867095.1 (alinhamento global).
4  data de download de cada sequencia: mtime UTC dos arquivos que a
   analise leu.
5  versao do multiz100way e do RefSeq curado do hg38: Last-Modified do
   arquivo no UCSC e dataTime da trilha na API.
6  a reliquia esta em GRCg7b e o mapa de Elferink 2010 em galGal3: ancora
   a reliquia no galGal3 por 32-meros exatos (mesma especie) e le o bin.

uv run --script scripts/atavismo/s37_pendencias_metodos.py
"""

from __future__ import annotations

import csv
import importlib.util
import statistics
import sys
from datetime import UTC, datetime
from pathlib import Path

import requests
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

RAIZ = Path(__file__).resolve().parents[2]
AT = RAIZ / "runs" / "atavismo"
OUT = AT / "pendencias"
WSL = Path(r"\\wsl.localhost\Ubuntu-24.04\home\manec\genomas")
UCSC = "https://api.genome.ucsc.edu"
K = 32


def importa(nome: str, arquivo: str):  # noqa: ANN201 — modulo carregado em tempo de execucao
    spec = importlib.util.spec_from_file_location(nome, Path(__file__).parent / arquivo)
    if spec is None or spec.loader is None:
        raise SystemExit(f"nao carreguei {arquivo}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[nome] = mod  # as dataclasses do s27 procuram o proprio modulo aqui
    spec.loader.exec_module(mod)
    return mod


s27 = importa("s27", "s27_painel_intervalos.py")
s28 = importa("s28", "s28_recombinacao_z.py")


def le_faa(caminho: Path) -> dict[str, tuple[str, str]]:
    """nome -> (cabecalho completo, sequencia)."""
    seqs: dict[str, tuple[str, list[str]]] = {}
    nome = ""
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        if linha.startswith(">"):
            nome = linha[1:].split()[0]
            seqs[nome] = (linha[1:].strip(), [])
        elif nome:
            seqs[nome][1].append(linha.strip())
    return {k: (c, "".join(s)) for k, (c, s) in seqs.items()}


def grava(nome: str, linhas: list[dict]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / nome).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)


def p2_anolis() -> None:
    print("\n=== 2. proteinas do Anolis da segunda consulta (s33) ===")
    usadas = le_faa(AT / "indice_nulo" / "consultas.faa")
    linhas = []
    for g in ["ENAM", "AMBN", "AMTN", "LPL"]:
        p = s27.gene_pos(g, 28377)
        gene = (
            next(
                (
                    x
                    for x in s27.tabela_features(p["acc"], p["ini"], p["fim"])
                    if x.sym.upper() == g and x.proteina
                ),
                None,
            )
            if p
            else None
        )
        if not gene or not gene.proteina:
            linhas.append({"gene": g, "acesso": "", "aa": 0, "identica_ao_s33": False})
            continue
        baixada = next(iter(le_faa_texto(s27.proteinas({g: gene.proteina})).values()))
        usada = usadas[f"{g}_Acar"][1]
        linhas.append(
            {
                "gene": g,
                "acesso": gene.proteina,
                "aa": len(baixada),
                "identica_ao_s33": baixada == usada,
            }
        )
        print(
            f"  {g:5s} {gene.proteina:16s} {len(baixada):5d} aa  identica a usada no s33: {baixada == usada}"
        )
    grava("s37_anolis_consultas.tsv", linhas)


def le_faa_texto(texto: str) -> dict[str, str]:
    tmp = OUT / "_tmp.faa"
    OUT.mkdir(parents=True, exist_ok=True)
    tmp.write_text(texto, encoding="utf-8")
    seqs = {k: s for k, (_, s) in le_faa(tmp).items()}
    tmp.unlink()
    return seqs


def blocos(com_gap: str, outra: str) -> str:
    """Corridas de gap em com_gap, na coordenada (1-based) dos residuos de outra."""
    saida, pos, ini = [], 0, None
    for a, b in zip(com_gap, outra, strict=True):
        if b != "-":
            pos += 1
        if a == "-" and ini is None:
            ini = pos
        elif a != "-" and ini is not None:
            saida.append(f"{ini}-{pos - 1}")
            ini = None
    if ini is not None:
        saida.append(f"{ini}-{pos}")
    return ";".join(saida)


def p3_genbank_refseq() -> None:
    print("\n=== 3. consultas GenBank do jacare x RefSeq da GCF_030867095.1 ===")
    usadas = le_faa(AT / "lncrna" / "jacare_proteinas.faa")
    proteoma = le_faa(AT / "painel" / "proteoma_Amis" / "protein.faa")
    alinhador = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    alinhador.substitution_matrix = substitution_matrices.load("BLOSUM62")
    linhas = []
    for nome in ["ENAM_Amis", "AMTN_Amis"]:
        cab, seq = usadas[nome]
        g = nome.split("_")[0]
        p = s27.gene_pos(g, 8496)
        gene = (
            next(
                (
                    x
                    for x in s27.tabela_features(p["acc"], p["ini"], p["fim"])
                    if x.sym.upper() == g and x.proteina
                ),
                None,
            )
            if p
            else None
        )
        if not gene or not gene.proteina or gene.proteina not in proteoma:
            print(f"  {g}: XP da anotacao nao achado no proteoma local")
            continue
        ref = proteoma[gene.proteina][1]
        aln = alinhador.align(seq, ref)[0]
        gb, rs = str(aln[0]), str(aln[1])
        pares = [(a, b) for a, b in zip(gb, rs, strict=True) if a != "-" and b != "-"]
        iguais = sum(a == b for a, b in pares)
        # identidade sobre as colunas pareadas; o que sobra e residuo que so uma tem
        linha = {
            "gene": g,
            "genbank": cab.split()[1],
            "aa_genbank": len(seq),
            "refseq": gene.proteina,
            "aa_refseq": len(ref),
            "pareados": len(pares),
            "identidade_pareados": round(iguais / len(pares), 4),
            "so_no_genbank": sum(b == "-" for b in rs),
            "so_no_refseq": sum(a == "-" for a in gb),
            "blocos_so_refseq": blocos(gb, rs),
            "blocos_so_genbank": blocos(rs, gb),
        }
        linhas.append(linha)
        print(
            f"  {g}: {linha['genbank']} ({len(seq)} aa) x {gene.proteina} ({len(ref)} aa): "
            f"{iguais}/{len(pares)} identicos nos pareados ({iguais / len(pares):.2%}); "
            f"so no GenBank {linha['so_no_genbank']} aa, so no RefSeq {linha['so_no_refseq']} aa"
        )
    if linhas:
        grava("s37_genbank_refseq.tsv", linhas)


def mtime(caminho: Path) -> str:
    try:
        return datetime.fromtimestamp(caminho.stat().st_mtime, UTC).strftime("%Y-%m-%d %H:%M")
    except OSError:
        return "inacessivel"


def p4_datas() -> None:
    print("\n=== 4. data de download (mtime UTC do arquivo lido pela analise) ===")
    linhas = []
    for f in sorted((AT / "lncrna").glob("intervalo_*.fna")):
        acesso = f.read_text(encoding="utf-8").split("\n", 1)[0].split()[1]
        linhas.append({"arquivo": f"lncrna/{f.name}", "sequencia": acesso, "baixado_utc": mtime(f)})
    painel = sorted((AT / "painel").glob("intervalo_*.fna"))
    datas = sorted(mtime(f) for f in painel)
    linhas.append(
        {
            "arquivo": f"painel/intervalo_*.fna ({len(painel)})",
            "sequencia": "intervalos do s27",
            "baixado_utc": f"{datas[0]} a {datas[-1]}",
        }
    )
    extra = [
        (AT / "painel" / "proteoma_Amis" / "protein.faa", "proteoma GCF_030867095.1"),
        (WSL / "Amis.fna", "genoma GCF_030867095.1"),
        (WSL / "Ggal.fna", "genoma GCF_016699485.2"),
        (WSL / "Shab.fna", "genoma GCF_004027225.2"),
        (WSL / "multiz" / "chr4.maf.gz", "multiz100way hg38 chr4"),
        (AT / "multiz" / "cds_hg38.bed", "ncbiRefSeqCurated hg38 (API UCSC)"),
        (AT / "recombinacao" / s28.ARQ, "Elferink 2010 S5 (Europe PMC)"),
    ]
    for f, rotulo in extra:
        linhas.append(
            {
                "arquivo": str(f.relative_to(AT))
                if f.is_relative_to(AT)
                else f"~/genomas/{f.relative_to(WSL)}",
                "sequencia": rotulo,
                "baixado_utc": mtime(f),
            }
        )
    for r in linhas:
        print(f"  {r['baixado_utc']:33s} {r['sequencia']:36s} {r['arquivo']}")
    grava("s37_datas_acesso.tsv", linhas)


def p5_versoes() -> None:
    print("\n=== 5. versoes do multiz100way e do RefSeq curado (hg38) ===")
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/chr4.maf.gz"
    h = requests.head(url, timeout=120, allow_redirects=True)
    h.raise_for_status()
    local = (WSL / "multiz" / "chr4.maf.gz").stat().st_size
    r = requests.get(
        f"{UCSC}/getData/track",
        params={
            "genome": "hg38",
            "track": "ncbiRefSeqCurated",
            "chrom": "chr4",
            "start": 70_500_000,
            "end": 70_510_000,
        },
        timeout=120,
    )
    r.raise_for_status()
    meta = {k: v for k, v in r.json().items() if not isinstance(v, (list, dict))}
    linhas = [
        {
            "item": "multiz100way chr4.maf.gz",
            "last_modified": h.headers.get("Last-Modified", ""),
            "bytes_servidor": h.headers.get("Content-Length", ""),
            "bytes_local": local,
        },
        {
            "item": "ncbiRefSeqCurated hg38",
            "last_modified": meta.get("dataTime", ""),
            "bytes_servidor": "",
            "bytes_local": "",
        },
    ]
    for x in linhas:
        print(
            f"  {x['item']:26s} {x['last_modified']}  servidor {x['bytes_servidor']}  local {x['bytes_local']}"
        )
    print(f"  metadados da trilha: {meta}")
    grava("s37_versoes_ucsc.tsv", linhas)


def kmers(seq: str) -> dict[str, int]:
    return {seq[i : i + K]: i for i in range(len(seq) - K + 1) if "N" not in seq[i : i + K]}


def p6_reliquia_galgal3() -> None:
    print("\n=== 6. reliquia do GRCg7b no galGal3 e o bin de Elferink 2010 ===")
    # galGal3 em torno da janela NRG1..LPL do s28 (52,84-53,41 Mb)
    ini3, fim3 = 52_500_000, 53_800_000
    r = requests.get(
        f"{UCSC}/getData/sequence",
        params={"genome": "galGal3", "chrom": "chrZ", "start": ini3, "end": fim3},
        timeout=300,
    )
    r.raise_for_status()
    alvo = r.json()["dna"].upper()
    indice: dict[str, list[int]] = {}
    for i in range(len(alvo) - K + 1):
        indice.setdefault(alvo[i : i + K], []).append(i)
    z = s28.tabela()
    z = z.loc[z["crom"] == "GGZ"]
    reliquias = {"ENAM": (55_155_480, 55_157_241), "AMBN": (55_168_855, 55_168_992)}
    flanco = 2_000  # a reliquia de AMBN tem 138 pb; o flanco da ancoras unicas de sobra
    linhas = []
    for gene, (a, b) in reliquias.items():
        fasta = s27.fatia("NC_052572.1", a - flanco, b + flanco, gene)
        if not fasta:
            raise SystemExit(f"efetch falhou para a reliquia de {gene}")
        q = "".join(fasta.split("\n")[1:]).upper()
        diagonais = {"+": [], "-": []}
        for fita, s in (("+", q), ("-", str(Seq(q).reverse_complement()))):
            for kmer, i in kmers(s).items():
                pos = indice.get(kmer, [])
                if len(pos) == 1:  # so ancora unica
                    diagonais[fita].append(pos[0] - i)
        fita = max(diagonais, key=lambda f: len(diagonais[f]))
        d = diagonais[fita]
        if len(d) < 5:
            print(f"  {gene}: {len(d)} ancoras — sem mapeamento")
            continue
        diag = statistics.median(d)
        concordantes = sum(abs(x - diag) <= 50 for x in d)
        # centro da reliquia: posicao (flanco + meio) na consulta, convertida pela diagonal
        meio_q = flanco + (b - a) // 2
        if fita == "-":
            meio_q = len(q) - 1 - meio_q
        centro3 = ini3 + int(diag) + meio_q + 1
        bin_ = z.loc[(z["de"] <= centro3) & (z["ate"] > centro3)].iloc[0]
        borda = min(centro3 - int(bin_["de"]), int(bin_["ate"]) - centro3)
        linha = {
            "reliquia": gene,
            "grcg7b": f"NC_052572.1:{a}-{b}",
            "fita": fita,
            "ancoras_32mer": len(d),
            "concordantes": concordantes,
            "centro_galgal3": centro3,
            "bin": int(bin_["bin"]),
            "bin_de": int(bin_["de"]),
            "bin_ate": int(bin_["ate"]),
            "dist_borda_kb": round(borda / 1e3, 1),
            "macho_cmmb": round(float(bin_["macho_cmmb"]), 2),
            "pop1_macho_cmmb": round(float(bin_["pop1_macho_cmmb"]), 2),
            "pop2_macho_cmmb": round(float(bin_["pop2_macho_cmmb"]), 2),
        }
        linhas.append(linha)
        print(
            f"  {gene}: fita {fita}, {concordantes}/{len(d)} ancoras na mesma diagonal -> galGal3 chrZ:"
            f"{centro3:,}; bin {linha['bin']} ({linha['bin_de']:,}-{linha['bin_ate']:,}), "
            f"{linha['dist_borda_kb']} kb da borda; macho {linha['macho_cmmb']} cM/Mb "
            f"(pop1 {linha['pop1_macho_cmmb']}, pop2 {linha['pop2_macho_cmmb']})"
        )
    if linhas:
        grava("s37_reliquia_galgal3.tsv", linhas)


def main() -> int:
    p2_anolis()
    p3_genbank_refseq()
    p4_datas()
    p5_versoes()
    p6_reliquia_galgal3()
    print(f"\ntabelas em {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
