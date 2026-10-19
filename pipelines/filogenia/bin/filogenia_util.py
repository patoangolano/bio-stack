#!/usr/bin/env python3
"""Utilitarios do pipeline de filogenia da bio-stack.

Roda dentro de container com Python 3.12 e so biblioteca padrao. O mesmo
arquivo atende tres processos do pipeline, e nenhum deles precisa de
dependencia externa — o que mantem o container minimo e a versao fixada.

Subcomandos:
  prepara    valida e saneia o FASTA de um loco e decide se ele segue
  concatena  monta a supermatriz e o arquivo de particoes NEXUS
  resumo     junta relatorios do IQ-TREE, concordancia e ASTER em TSV/JSON/MD
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import statistics
import sys
from collections.abc import Iterable
from pathlib import Path

# O IQ-TREE (medido no 3.1.4) troca por "_" todo caractere fora de
# [A-Za-z0-9_.\-|/] e so avisa no log; espaco no cabecalho FASTA ele nao
# troca: corta o nome e so junta palavras com "_" enquanto houver nome
# repetido no arquivo. Sanear aqui, antes do alinhamento, faz o nome do taxon
# ser o mesmo em todos os locos — e a concatenacao casa linhas pelo nome. Este
# conjunto e mais estrito ('|' e '/' tambem viram "_").
_ILEGAL = re.compile(r"[^A-Za-z0-9_.\-]+")
_SUBLINHADOS = re.compile(r"_+")

_DNA = set("ACGTUN")
_IUPAC_DNA = set("ACGTUNRYSWKMBDHV")
_AA = set("ACDEFGHIKLMNPQRSTVWYBJZXUO")
_GAPS = set("-.?")


class ErroEntrada(Exception):
    """Entrada invalida que o usuario precisa corrigir; vira mensagem, nao traceback."""


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------


def ler_fasta(caminho: Path) -> list[tuple[str, str]]:
    """Le FASTA (opcionalmente .gz), tolerando CRLF, linhas quebradas e espacos."""
    abrir = gzip.open if caminho.suffix.lower() == ".gz" else open
    registros: list[tuple[str, str]] = []
    nome: str | None = None
    partes: list[str] = []
    # utf-8-sig: o Set-Content -Encoding UTF8 do PowerShell 5.1 grava BOM, e com
    # "utf-8" a linha 1 virava "\ufeff>nome" e o arquivo era recusado.
    with abrir(caminho, "rt", encoding="utf-8-sig", errors="replace") as fh:
        for n, bruta in enumerate(fh, start=1):
            linha = bruta.strip()
            if not linha:
                continue
            if linha.startswith(">"):
                if nome is not None:
                    registros.append((nome, "".join(partes)))
                nome = linha[1:].strip()
                partes = []
            else:
                if nome is None:
                    raise ErroEntrada(
                        f"{caminho.name}: linha {n} tem sequencia antes do primeiro '>'"
                    )
                partes.append("".join(linha.split()))
    if nome is not None:
        registros.append((nome, "".join(partes)))
    return registros


def escrever_fasta(caminho: Path, registros: Iterable[tuple[str, str]], largura: int = 80) -> None:
    with caminho.open("w", encoding="utf-8", newline="\n") as fh:
        for nome, seq in registros:
            fh.write(f">{nome}\n")
            for i in range(0, len(seq), largura):
                fh.write(seq[i : i + largura] + "\n")


def sanear_nome(cabecalho: str, modo: str) -> str:
    """Extrai o nome do taxon do cabecalho e troca o que o IQ-TREE recusaria."""
    if modo == "colchete":
        # Dois padroes do NCBI. O classico (protein.faa do RefSeq):
        #   ">XP_015134201.2 enamelin [Gallus gallus]"
        # e o do NCBI Datasets, com chave=valor e o organismo no MEIO:
        #   ">NP_114095.2 ENAM [organism=Homo sapiens] [GeneID=10117] [isoform=1]"
        # Pegar o ultimo colchete daria "isoform_1" neste exemplo e "GeneID_<n>"
        # na proteina sem isoforma, que e a maioria. Nenhum dos dois e taxon: o
        # GeneID e um nome por gene (nenhum taxon casaria entre locos), e o
        # rotulo de isoforma se repete entre especies (nome repetido no loco).
        # Colchete so com chave que nao e organism nao serve: vira cabecalho sem
        # nome, e da erro.
        achados = re.findall(r"\[([^\[\]]+)\]", cabecalho)
        organismo = [
            a.split("=", 1)[1] for a in achados if a.strip().lower().startswith("organism=")
        ]
        sem_chave = [a for a in achados if "=" not in a]
        if achados:
            bruto = (organismo or sem_chave or [""])[-1]
        else:
            bruto = cabecalho.split()[0] if cabecalho.split() else ""
    elif modo == "completo":
        bruto = cabecalho
    else:  # primeiro
        bruto = cabecalho.split()[0] if cabecalho.split() else ""
    limpo = _SUBLINHADOS.sub("_", _ILEGAL.sub("_", bruto.strip())).strip("_")
    return limpo


def detectar_tipo(seqs: Iterable[str]) -> str:
    """DNA se >= 90% dos caracteres informativos forem ACGTUN; senao proteina."""
    total = dna = 0
    for s in seqs:
        # upper(): o MAFFT grava DNA em minusculas; sem isto todo loco alinhado
        # contava como proteina e a trava de tipos misturados nunca disparava.
        for c in s.upper():
            if c in _GAPS:
                continue
            total += 1
            if c in _DNA:
                dna += 1
    if total == 0:
        return "DNA"
    return "DNA" if dna / total >= 0.9 else "AA"


def _lista(texto: str) -> list[str]:
    return [t.strip() for t in texto.split(",") if t.strip()]


# --------------------------------------------------------------------------
# prepara
# --------------------------------------------------------------------------


def cmd_prepara(a: argparse.Namespace) -> int:
    registros = ler_fasta(Path(a.fasta))
    if not registros:
        raise ErroEntrada(f"{a.fasta}: nenhum registro FASTA encontrado")

    renomeados: dict[str, str] = {}
    avisos: list[str] = []
    descartados: list[str] = []
    vistos: dict[str, str] = {}
    limpos: list[tuple[str, str]] = []

    for cab, seq_bruta in registros:
        nome = sanear_nome(cab, a.nome_taxon)
        if not nome:
            raise ErroEntrada(f"cabecalho sem nome utilizavel: '>{cab}'")
        if nome != cab:
            renomeados[cab] = nome
        if nome in vistos:
            raise ErroEntrada(
                f"nome de taxon repetido depois do saneamento: '{nome}' "
                f"(cabecalhos '>{vistos[nome]}' e '>{cab}'). A concatenacao casa "
                "linhas pelo nome; deixe um registro por taxon por loco."
            )
        vistos[nome] = cab

        seq = seq_bruta.upper()
        if a.desalinhar:
            seq = seq.replace("-", "").replace(".", "")
        if "*" in seq:
            # rstrip com gap: num alinhamento o '*' final vem seguido de "---".
            internos = seq.rstrip("*-.").count("*")
            if internos:
                avisos.append(
                    f"{nome}: {internos} '*' interno(s) removido(s) — pseudogene ou erro de traducao?"
                )
            # Desalinhado, o '*' sai; ja alinhado, vira gap — apagar deslocaria as
            # colunas e o loco seria recusado como "nao alinhado", a causa errada.
            seq = seq.replace("*", "" if a.desalinhar else "-")
        if not seq.strip("-.?"):
            descartados.append(nome)
            avisos.append(f"{nome}: sequencia vazia depois da limpeza; registro descartado")
            continue
        limpos.append((nome, seq))

    detectado = detectar_tipo(s for _, s in limpos)
    tipo = a.tipo.upper() if a.tipo else detectado
    if a.tipo and tipo != detectado:
        # A, C, G, T e N tambem sao letras de aminoacido: sem esta trava o IQ-TREE
        # ajustaria um modelo de proteina a nucleotideos sem reclamar.
        raise ErroEntrada(
            f"--tipo_seq {tipo}, mas as sequencias parecem {detectado}. Forcar o tipo nao "
            "traduz nada: para a arvore de proteina, traduza antes (seqkit translate -T "
            "<codigo genetico>) e rode sobre os FASTA de proteina."
        )
    permitidos = (_IUPAC_DNA if tipo == "DNA" else _AA) | _GAPS
    for nome, seq in limpos:
        estranhos = sorted(set(seq) - permitidos)
        if estranhos:
            raise ErroEntrada(
                f"{nome}: caracteres invalidos para {tipo}: {''.join(estranhos)} "
                "(se o tipo detectado estiver errado, passe --tipo_seq)"
            )

    if not a.desalinhar:
        comprimentos = {len(s) for _, s in limpos}
        if len(comprimentos) > 1:
            raise ErroEntrada(
                f"alinhador = nenhum, mas as sequencias tem comprimentos diferentes "
                f"({min(comprimentos)} a {max(comprimentos)}): a entrada nao esta alinhada"
            )

    nomes = [n for n, _ in limpos]
    outgroup = [sanear_nome(o, "completo") for o in _lista(a.outgroup or "")]
    og_presentes = [o for o in outgroup if o in vistos and o not in descartados]
    status = "ok" if len(limpos) >= a.min_taxa else "poucos_taxa"

    prefixo = a.prefixo
    escrever_fasta(Path(f"{prefixo}.fa"), limpos)
    comprimentos = [len(s.replace("-", "")) for _, s in limpos]
    relatorio = {
        "loco": a.loco,
        "status": status,
        "n_taxa": len(limpos),
        "min_taxa": a.min_taxa,
        "tipo": tipo,
        "comprimento_min": min(comprimentos) if comprimentos else 0,
        "comprimento_mediano": statistics.median(comprimentos) if comprimentos else 0,
        "comprimento_max": max(comprimentos) if comprimentos else 0,
        "taxa": nomes,
        "outgroup_presentes": og_presentes,
        "renomeados": renomeados,
        "descartados": descartados,
        "avisos": avisos,
    }
    Path(f"{prefixo}.json").write_text(
        json.dumps(relatorio, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    # Lido pelo processo com `source`: e assim que o Nextflow recebe as
    # decisoes do script como saidas `env` sem precisar abrir JSON no Groovy.
    Path(f"{prefixo}.env").write_text(
        f"N_TAXA={len(limpos)}\nTIPO={tipo}\nSTATUS={status}\n"
        f"OUTGROUP_PRESENTES={','.join(og_presentes)}\n",
        encoding="utf-8",
    )
    for aviso in avisos:
        print(f"[aviso] {a.loco}: {aviso}", file=sys.stderr)
    if status != "ok":
        print(
            f"[aviso] {a.loco}: {len(limpos)} taxa (< {a.min_taxa}); loco fica fora da inferencia",
            file=sys.stderr,
        )
    return 0


# --------------------------------------------------------------------------
# concatena
# --------------------------------------------------------------------------


def _charset(nome: str) -> str:
    """Identificador NEXUS valido: so letra, digito e sublinhado, sem comecar por digito."""
    limpo = re.sub(r"[^A-Za-z0-9_]", "_", nome)
    return limpo if limpo[:1].isalpha() else f"l_{limpo}"


def cmd_concatena(a: argparse.Namespace) -> int:
    locos = a.locos
    arquivos = [Path(p) for p in a.alinhamentos]
    if len(locos) != len(arquivos):
        raise ErroEntrada(f"{len(locos)} nomes de loco para {len(arquivos)} alinhamentos")

    blocos: list[tuple[str, int, dict[str, str]]] = []
    tipos: set[str] = set()
    for loco, arq in zip(locos, arquivos, strict=True):
        regs = ler_fasta(arq)
        tamanhos = {len(s) for _, s in regs}
        if len(tamanhos) != 1:
            raise ErroEntrada(
                f"{arq.name}: alinhamento com sequencias de comprimentos diferentes {sorted(tamanhos)}"
            )
        tipos.add(detectar_tipo(s for _, s in regs))
        blocos.append((loco, tamanhos.pop(), {n: s for n, s in regs}))
    if len(tipos) > 1:
        raise ErroEntrada("locos de DNA e de proteina misturados; rode cada tipo em separado")

    taxa = sorted({t for _, _, seqs in blocos for t in seqs})
    supermatriz = {t: [] for t in taxa}
    charsets: list[str] = []
    inicio = 1
    for loco, comp, seqs in blocos:
        for t in taxa:
            supermatriz[t].append(seqs.get(t, "-" * comp))
        fim = inicio + comp - 1
        charsets.append(f"  charset {_charset(loco)} = {inicio}-{fim};")
        inicio = fim + 1

    escrever_fasta(Path(a.saida_fasta), ((t, "".join(supermatriz[t])) for t in taxa))
    Path(a.saida_particoes).write_text(
        "#nexus\nbegin sets;\n" + "\n".join(charsets) + "\nend;\n", encoding="utf-8", newline="\n"
    )

    with Path(a.saida_ocupacao).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["taxon", *locos, "n_locos"])
        for t in taxa:
            pres = [int(t in seqs) for _, _, seqs in blocos]
            w.writerow([t, *pres, sum(pres)])

    outgroup = [sanear_nome(o, "completo") for o in _lista(a.outgroup or "")]
    presentes = [o for o in outgroup if o in supermatriz]
    celulas = len(taxa) * len(blocos)
    ocupadas = sum(1 for _, _, seqs in blocos for t in taxa if t in seqs)
    Path(a.saida_env).write_text(
        f"OUTGROUP_PRESENTES={','.join(presentes)}\n"
        f"N_TAXA={len(taxa)}\nN_LOCOS={len(blocos)}\nCOMPRIMENTO={inicio - 1}\n"
        f"OCUPACAO={ocupadas / celulas:.4f}\n",
        encoding="utf-8",
    )
    return 0


# --------------------------------------------------------------------------
# resumo
# --------------------------------------------------------------------------

_RX = {
    # Supermatriz: "Input data: 12 taxa with 4 partitions and 4053 total sites (2.1% missing data)"
    "particionado": re.compile(
        r"Input data:\s+(\d+)\s+taxa\s+with\s+(\d+)\s+partitions\s+and\s+(\d+)\s+total sites"
        r"(?:\s+\(([\d.]+)% missing data\))?"
    ),
    # Loco unico: "Input data: 12 sequences with 936 nucleotide sites"
    "entrada": re.compile(r"Input data:\s+(\d+)\s+(?:taxa|sequences)\s+with\s+(\d+)\s+([\w-]+)"),
    "constantes": re.compile(r"Number of constant sites:\s+(\d+)"),
    "informativos": re.compile(r"Number of parsimony informative sites:\s+(\d+)"),
    "padroes": re.compile(r"Number of distinct site patterns:\s+(\d+)"),
    "modelo": re.compile(r"Best-fit model according to (\w+):\s+(\S+)"),
    "lnl": re.compile(r"Log-likelihood of the tree:\s+(-?[\d.]+)"),
    "comprimento": re.compile(r"Total tree length \(sum of branch lengths\):\s+([\d.]+)"),
    "internos": re.compile(r"Sum of internal branch lengths:\s+([\d.]+)"),
}


def ler_relatorio_iqtree(caminho: Path) -> dict[str, object]:
    texto = caminho.read_text(encoding="utf-8", errors="replace")
    out: dict[str, object] = {}
    if m := _RX["particionado"].search(texto):
        # Com MFP+MERGE esta linha descreve o esquema JA FUNDIDO, nao a entrada:
        # quatro locos fundidos em um viram "1 partitions ... (0% missing data)".
        # Contagem de entrada e dado faltante vem de particoes.nex e ocupacao.tsv.
        out["n_seqs"], out["n_particoes_relatorio"], out["n_sitios"] = (
            int(m[1]),
            int(m[2]),
            int(m[3]),
        )
        if m[4] is not None:
            out["perc_faltantes_relatorio"] = float(m[4])
    elif m := _RX["entrada"].search(texto):
        out["n_seqs"], out["n_sitios"], out["tipo_sitio"] = int(m[1]), int(m[2]), m[3]
    for chave in ("constantes", "informativos", "padroes"):
        if m := _RX[chave].search(texto):
            out[chave] = int(m[1])
    if m := _RX["modelo"].search(texto):
        out["criterio"], out["modelo"] = m[1], m[2]
    if m := _RX["lnl"].search(texto):
        out["lnL"] = float(m[1])
    if m := _RX["comprimento"].search(texto):
        out["comprimento_arvore"] = float(m[1])
    if m := _RX["internos"].search(texto):
        out["comprimento_internos"] = float(m[1])
    return out


_ROTULO = re.compile(r"\)([^:,;()]+)")


def ordem_rotulos(params: dict[str, object]) -> list[str]:
    """
    Campos do rotulo interno do IQ-TREE, na ordem em que ele os escreve:
    SH-aLRT/aBayes/UFBoot, so os que foram pedidos. Sem a ordem, "95/100" e
    ambiguo — pode ser SH/UFBoot ou aBayes/UFBoot.
    """
    return _ordem_por_params(params)


def suportes(newick: str, ordem: list[str]) -> list[dict[str, float]]:
    """Um dicionario {campo: valor} por no interno; rotulo com outro numero de campos e ignorado."""
    nos: list[dict[str, float]] = []
    for rot in _ROTULO.findall(newick):
        partes = rot.split("/")
        # aBayes sem SH-aLRT: o IQ-TREE deixa a casa do SH-aLRT vazia ("/0.98").
        if len(partes) > 1 and partes[0] == "":
            partes = partes[1:]
        if not ordem or len(partes) != len(ordem):
            continue
        try:
            nos.append({campo: float(v) for campo, v in zip(ordem, partes, strict=True)})
        except ValueError:
            continue
    return nos


def estat_suporte(newick: str, ordem: list[str]) -> dict[str, object]:
    nos = suportes(newick, ordem)
    out: dict[str, object] = {"n_ramos_internos": len(nos)}
    ufb = [n["ufboot"] for n in nos if "ufboot" in n]
    ab = [n["abayes"] for n in nos if "abayes" in n]
    if ufb:
        out["ufboot_medio"] = round(statistics.fmean(ufb), 2)
        out["frac_ufboot_95"] = round(sum(u >= 95 for u in ufb) / len(ufb), 3)
    if ab:
        out["abayes_medio"] = round(statistics.fmean(ab), 3)
        out["frac_abayes_095"] = round(sum(v >= 0.95 for v in ab) / len(ab), 3)
    ambos = [(n["alrt"], n["ufboot"]) for n in nos if "alrt" in n and "ufboot" in n]
    if ambos:
        # Criterio conjunto da documentacao do IQ-TREE: ramo confiavel quando
        # SH-aLRT >= 80% E UFBoot >= 95%.
        out["frac_sh80_ufb95"] = round(sum(s >= 80 and u >= 95 for s, u in ambos) / len(ambos), 3)
    return out


def ler_cf(caminho: Path) -> dict[str, object]:
    linhas = [
        ln
        for ln in caminho.read_text(encoding="utf-8").splitlines()
        if ln and not ln.startswith("#")
    ]
    if not linhas:
        return {}
    leitor = csv.DictReader(linhas, delimiter="\t")
    out: dict[str, object] = {}
    valores: dict[str, list[float]] = {}
    for linha in leitor:
        for col in ("gCF", "sCF"):
            v = linha.get(col)
            if v not in (None, "", "NA"):
                valores.setdefault(col, []).append(float(v))
    for col, vs in valores.items():
        out[f"{col}_medio"] = round(statistics.fmean(vs), 2)
        out[f"{col}_min"] = round(min(vs), 2)
        out[f"n_ramos_{col}"] = len(vs)
    return out


def ler_versoes(caminho: Path | None) -> dict[str, str]:
    """versoes.yml no formato do pipeline: 'PROCESSO:' e '  ferramenta: versao'."""
    versoes: dict[str, str] = {}
    if caminho is None or not caminho.is_file():
        return versoes
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        if linha.startswith("  ") and ":" in linha:
            ferramenta, versao = (p.strip() for p in linha.strip().split(":", 1))
            versoes.setdefault(ferramenta, versao)
    return versoes


def _comprimento_alinhamento(caminho: Path) -> int | None:
    regs = ler_fasta(caminho)
    return len(regs[0][1]) if regs else None


_REFERENCIAS = {
    "iqtree": "Wong TKF, Ly-Trong N, Ren H, et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Molecular Biology and Evolution 43(5):msag117. doi:10.1093/molbev/msag117",
    "modelfinder": "Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14:587-589. doi:10.1038/nmeth.4285",
    "ufboot": "Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS (2018) UFBoot2: improving the ultrafast bootstrap approximation. Molecular Biology and Evolution 35:518-522. doi:10.1093/molbev/msx281",
    "shalrt": "Guindon S, Dufayard JF, Lefort V, et al. (2010) New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology 59:307-321. doi:10.1093/sysbio/syq010",
    "abayes": "Anisimova M, Gil M, Dufayard JF, Dessimoz C, Gascuel O (2011) Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. Systematic Biology 60:685-699. doi:10.1093/sysbio/syr041",
    "particao": "Chernomor O, von Haeseler A, Minh BQ (2016) Terrace aware data structure for phylogenomic inference from supermatrices. Systematic Biology 65:997-1008. doi:10.1093/sysbio/syw037",
    "rcluster": "Lanfear R, Calcott B, Kainer D, Mayer C, Stamatakis A (2014) Selecting optimal partitioning schemes for phylogenomic datasets. BMC Evolutionary Biology 14:82. doi:10.1186/1471-2148-14-82",
    "partitionfinder2": "Lanfear R, Frandsen PB, Wright AM, Senfeld T, Calcott B (2017) PartitionFinder 2: new methods for selecting partitioned models of evolution for molecular and morphological phylogenetic analyses. Molecular Biology and Evolution 34:772-773. doi:10.1093/molbev/msw260",
    "gcf": "Minh BQ, Hahn MW, Lanfear R (2020) New methods to calculate concordance factors for phylogenomic datasets. Molecular Biology and Evolution 37:2727-2733. doi:10.1093/molbev/msaa106",
    "scfl": "Mo YK, Lanfear R, Hahn MW, Minh BQ (2023) Updated site concordance factors minimize effects of homoplasy and taxon sampling. Bioinformatics 39(1):btac741. doi:10.1093/bioinformatics/btac741",
    "mafft": "Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution 30:772-780. doi:10.1093/molbev/mst010",
    "famsa": "Deorowicz S, Debudaj-Grabysz A, Gudys A (2016) FAMSA: fast and accurate multiple sequence alignment of huge protein families. Scientific Reports 6:33964. doi:10.1038/srep33964",
    "clipkit": "Steenwyk JL, Buida TJ, Li Y, Shen XX, Rokas A (2020) ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. PLoS Biology 18:e3001007. doi:10.1371/journal.pbio.3001007",
    "trimal": "Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T (2009) trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics 25:1972-1973. doi:10.1093/bioinformatics/btp348",
    "aster": "Zhang C, Nielsen R, Mirarab S (2025) ASTER: a package for large-scale phylogenomic reconstructions. Molecular Biology and Evolution 42(8):msaf172. doi:10.1093/molbev/msaf172",
    "tan2015": "Tan G, Muffato M, Ledergerber C, et al. (2015) Current methods for automated filtering of multiple sequence alignments frequently worsen single-gene phylogenetic inference. Systematic Biology 64:778-791. doi:10.1093/sysbio/syv033",
    "wastral": "Zhang C, Mirarab S (2022) Weighting by gene tree uncertainty improves accuracy of quartet-based species trees. Molecular Biology and Evolution 39:msac215. doi:10.1093/molbev/msac215",
    "nextflow": "Di Tommaso P, Chatzou M, Floden EW, et al. (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319. doi:10.1038/nbt.3820",
    "nfcore": "Ewels PA, Peltzer A, Fillinger S, et al. (2020) The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology 38:276-278. doi:10.1038/s41587-020-0439-x",
}

#: Marcador de versao que nao veio do versoes.yml. O teste procura esta string.
PENDENTE = "[VERSÃO A CONFIRMAR]"


def _v(versoes: dict[str, str], ferramenta: str) -> str:
    return versoes.get(ferramenta, PENDENTE)


def _ligado(valor: object) -> bool:
    """
    Verdade de parametro vindo de JSON, tolerante a texto. O Nextflow 26.04 sem
    bloco params tipado entregava "false" e "0" como String — verdadeiras. O
    pipeline agora tipa os parametros, mas execucoes antigas ficaram assim.
    """
    if isinstance(valor, str):
        return valor.strip().lower() not in ("", "0", "false", "no", "nao", "não", "null", "none")
    return bool(valor)


def _ordem_por_params(params: dict[str, object]) -> list[str]:
    pedidos = (
        ("alrt", params.get("alrt")),
        ("abayes", params.get("abayes")),
        ("ufboot", params.get("ufboot")),
    )
    return [campo for campo, ligado in pedidos if _ligado(ligado)]


def _usa_modelfinder(modelo: str) -> bool:
    """ModelFinder so roda com -m MF, MFP, TEST e variantes (MFP+MERGE, TESTONLY...)."""
    m = modelo.upper()
    return m.startswith(("MF", "TEST"))


def _texto_limpeza(modo: str, versao: str) -> str:
    if modo in ("smart-gap", "gappy"):
        return (
            f"Colunas com excesso de gaps foram removidas com o ClipKIT v{versao} (modo {modo}"
            + (", limiar de gaps determinado dinamicamente)" if modo == "smart-gap" else ")")
            + "."
        )
    if modo.startswith("kpic"):
        return (
            f"Os alinhamentos foram aparados com o ClipKIT v{versao} (modo {modo}), que mantém "
            "só os sítios parcimônia-informativos e os constantes."
        )
    if modo.startswith("kpi"):
        return (
            f"Os alinhamentos foram aparados com o ClipKIT v{versao} (modo {modo}), que mantém "
            "só os sítios parcimônia-informativos."
        )
    return f"Os alinhamentos foram aparados com o ClipKIT v{versao} (modo {modo})."


def montar_metodos(
    params: dict[str, object], versoes: dict[str, str], resumo: dict[str, object]
) -> str:
    """
    Paragrafo de Metodos em portugues. Cada frase depende do que DE FATO rodou:
    versoes lidas do versoes.yml (nunca inventadas), ModelFinder so com -m
    MF/MFP/TEST, suporte so os pedidos, peso do wASTRAL pela mesma regra do
    main.nf (rotuloSuporte), texto da limpeza pelo modo do ClipKIT.
    """
    alinhador = str(params.get("alinhador", "mafft"))
    limpeza = str(params.get("limpeza", "clipkit"))
    modelo = str(params.get("modelo", "MFP"))
    amostragem = str(params.get("amostragem_concat") or "")
    n_locos = resumo.get("n_locos_inferidos", 0)
    usados = ["iqtree", "nextflow", "nfcore"]
    partes: list[str] = []

    if alinhador == "mafft":
        partes.append(
            f"As sequências de cada loco foram alinhadas com o MAFFT v{_v(versoes, 'mafft')} "
            "(L-INS-i até 200 sequências e estratégia automática acima disso; refinamento "
            "iterativo em uma thread, para resultado determinístico)."
        )
        usados.append("mafft")
    elif alinhador == "famsa":
        partes.append(
            f"As sequências de cada loco foram alinhadas com o FAMSA v{_v(versoes, 'famsa')}."
        )
        usados.append("famsa")
    else:
        partes.append("Os alinhamentos fornecidos foram usados sem realinhamento.")

    if limpeza == "clipkit":
        partes.append(
            _texto_limpeza(str(params.get("clipkit_modo", "smart-gap")), _v(versoes, "clipkit"))
        )
        usados.append("clipkit")
    elif limpeza == "trimal":
        partes.append(
            f"Os alinhamentos foram aparados com o trimAl v{_v(versoes, 'trimal')} (-automated1)."
        )
        usados.append("trimal")

    suporte = []
    if _ligado(params.get("alrt")):
        suporte.append(f"SH-aLRT ({params['alrt']} réplicas)")
        usados.append("shalrt")
    if _ligado(params.get("abayes")):
        suporte.append("aBayes")
        usados.append("abayes")
    if _ligado(params.get("ufboot")):
        suporte.append(
            f"UFBoot2 ({params['ufboot']} réplicas"
            + (", com otimização --bnni)" if _ligado(params.get("bnni")) else ")")
        )
        usados.append("ufboot")
    lista_suporte = (
        ", ".join(suporte[:-1]) + (" e " if len(suporte) > 1 else "") + suporte[-1]
        if suporte
        else ""
    )
    if _usa_modelfinder(modelo):
        escolha_modelo = (
            f"com o modelo de substituição escolhido pelo ModelFinder (-m {modelo}) segundo o BIC"
        )
        usados.append("modelfinder")
    else:
        escolha_modelo = f"sob o modelo {modelo}, definido pelo usuário"
    partes.append(
        ("Para o único loco" if str(n_locos) == "1" else f"Para cada um dos {n_locos} locos")
        + ", a árvore de máxima verossimilhança foi inferida no "
        f"IQ-TREE v{_v(versoes, 'iqtree')}, {escolha_modelo}"
        + (f"; o suporte dos ramos foi estimado por {lista_suporte}" if suporte else "")
        + "."
    )

    concat = resumo.get("concatenado")
    if isinstance(concat, dict):
        n_part = concat.get("n_particoes_entrada", n_locos)
        fundidas = concat.get("n_particoes_finais")
        ocupacao = concat.get("ocupacao_perc")
        frase = (
            ("Os alinhamentos limpos" if limpeza in ("clipkit", "trimal") else "Os alinhamentos")
            + " foram concatenados numa supermatriz de "
            f"{concat.get('n_taxa', concat.get('n_seqs', '?'))} táxons e {concat.get('n_sitios', '?')} "
            f"sítios, particionada por loco ({_plural(n_part, 'partição', 'partições')}"
            + (
                f"; {ocupacao}% das células táxon × loco preenchidas"
                if ocupacao is not None
                else ""
            )
            + "), e analisada com modelo de partição de comprimentos de ramo proporcionalmente ligados"
        )
        if _usa_modelfinder(modelo) and amostragem:
            # Com reamostragem por loco o pipeline roda -m MFP (sem MERGE): cada
            # loco fica com o seu modelo e a unidade sorteada e o loco.
            frase += (
                ", com um modelo escolhido por loco (sem fusão de partições) e UFBoot com "
                f"reamostragem {amostragem} sobre os {n_part} locos"
            )
        elif _usa_modelfinder(modelo):
            frase += " e fusão de partições pelo ModelFinder (MFP+MERGE, algoritmo rclusterf)" + (
                f"; o esquema final teve {_plural(fundidas, 'partição', 'partições')}"
                if fundidas
                else ""
            )
            usados += ["rcluster", "partitionfinder2"]
        partes.append(frase + ", com os mesmos critérios de suporte.")
        usados.append("particao")
    if resumo.get("concordancia"):
        partes.append(
            "Os fatores de concordância de genes (gCF) e de sítios por verossimilhança (sCFL, "
            f"{params.get('scfl', 100)} quartetos por ramo) foram calculados no IQ-TREE sobre a árvore concatenada."
        )
        usados += ["gcf", "scfl"]
    if resumo.get("arvore_especies"):
        # Mesma regra do rotuloSuporte() do main.nf.
        if _ligado(params.get("abayes")):
            peso = "pelo suporte aBayes e pelo comprimento de ramo (modo híbrido)"
        elif _ligado(params.get("ufboot")):
            peso = "pelo suporte UFBoot e pelo comprimento de ramo (modo híbrido)"
        elif _ligado(params.get("alrt")):
            peso = "pelo suporte SH-aLRT e pelo comprimento de ramo (modo híbrido)"
        else:
            peso = "só pelo comprimento de ramo (--mode 3)"
        partes.append(
            f"A árvore de espécies foi estimada com o wASTRAL (ASTER v{_v(versoes, 'aster')}) a partir das "
            f"árvores de genes, ponderando os quartetos {peso}."
        )
        usados += ["wastral", "aster"]
    threads_concat = concat.get("threads") if isinstance(concat, dict) else None
    threads_texto = (
        f"{threads_concat}" if threads_concat else str(params.get("cpus_concatenado", 1))
    )
    partes.append(
        f"O fluxo foi executado no Nextflow v{_v(versoes, 'nextflow')} com módulos nf-core em containers "
        f"BioContainers, semente aleatória {params.get('seed')} e IQ-TREE em uma thread por loco"
        # Sem supermatriz (um loco so, ou --concatenar false) nao ha thread dela a declarar.
        + (f" (supermatriz: {threads_texto})." if isinstance(concat, dict) else ".")
    )
    pendencias = []
    if PENDENTE in " ".join(partes):
        pendencias.append(f"- Há versões marcadas {PENDENTE}: o versoes.yml não as registrou.")
    if threads_concat and int(str(threads_concat)) > 1:
        pendencias.append(
            "- A supermatriz rodou com mais de uma thread: a topologia e o lnL se reproduzem, mas "
            "UFBoot e SH-aLRT podem variar 1 a 2 pontos entre execuções e os comprimentos de ramo "
            "nas últimas casas decimais. Declare isso no texto."
        )
    avisos = resumo.get("avisos_preparo") or {}
    if isinstance(avisos, dict) and avisos:
        pendencias.append(
            f"- O preparo alterou registros em {_plural(len(avisos), 'loco', 'locos')} "
            f"({', '.join(map(str, avisos))}): '*' interno removido ou registro vazio "
            "descartado. Veja resumo.json → avisos_preparo antes de usar o resultado — "
            "stop interno numa proteína traduzida é pseudogene ou erro de tradução."
        )
    if resumo.get("outgroup_ausente"):
        pendencias.append(
            f"- O outgroup pedido ({params.get('outgroup')}) não aparece em nenhum loco: as árvores "
            "ficaram sem raiz definida e as figuras, no ponto médio. Confira a grafia e o "
            "--nome_taxon."
        )
    if isinstance(concat, dict):
        ocupacao = concat.get("ocupacao_perc")
        if isinstance(ocupacao, int | float) and ocupacao < 50:
            pendencias.append(
                f"- Ocupação táxon × loco de {ocupacao}% (abaixo de 50%): a supermatriz é mais "
                "gap que dado. Se os nomes de táxon são acessos (XP_…, NP_…), o --nome_taxon "
                "está errado para esse cabeçalho; senão, declare a ocupação no texto."
            )
        if concat.get("n_taxa_multiloco") == 0:
            pendencias.append(
                "- Nenhum táxon aparece em mais de um loco: a supermatriz não liga os locos, e "
                "a árvore concatenada e a de espécies não têm sentido. Quase sempre é nome de "
                "táxon diferente entre locos — veja o --nome_taxon."
            )
    agressivos = resumo.get("locos_limpeza_agressiva") or []
    if isinstance(agressivos, list) and agressivos:
        pendencias.append(
            f"- Locos com mais de 20% dos sítios removidos na limpeza: {', '.join(map(str, agressivos))}. "
            "Limpeza agressiva tende a piorar árvores de gene único (Tan et al. 2015); confira se "
            "a conclusão sobrevive a --limpeza nenhuma."
        )
        usados.append("tan2015")
    refs = "\n".join(f"- {_REFERENCIAS[k]}" for k in dict.fromkeys(usados))
    return (
        "# Métodos — filogenia\n\n"
        + " ".join(partes)
        + "\n\n## Referências\n\n"
        + refs
        + "\n\n## Pendências\n\n"
        + ("\n".join(pendencias) if pendencias else "- Nenhuma.")
        + "\n"
    )


def _plural(n: object, singular: str, plural: str) -> str:
    return f"{n} {singular if str(n) == '1' else plural}"


def _primeiro(diretorio: Path | None, padrao: str) -> Path | None:
    if diretorio is None or not diretorio.is_dir():
        return None
    return next(iter(sorted(diretorio.glob(padrao))), None)


def _contar_charsets(caminho: Path) -> int:
    return sum(
        ln.strip().lower().startswith("charset")
        for ln in caminho.read_text(encoding="utf-8").splitlines()
    )


def ler_ocupacao(caminho: Path) -> dict[str, object]:
    """ocupacao.tsv do CONCATENA: taxon x loco com 1/0 e a coluna n_locos no fim."""
    with caminho.open(encoding="utf-8", newline="") as fh:
        linhas = list(csv.reader(fh, delimiter="\t"))
    if len(linhas) < 2:
        return {}
    n_locos = len(linhas[0]) - 2
    celulas = presentes = multiloco = 0
    for linha in linhas[1:]:
        valores = [int(v) for v in linha[1 : 1 + n_locos]]
        celulas += len(valores)
        presentes += sum(valores)
        multiloco += sum(valores) >= 2
    return {
        "n_taxa": len(linhas) - 1,
        "ocupacao_perc": round(100 * presentes / celulas, 1) if celulas else None,
        # Taxon em um loco so nao liga nada na supermatriz. Zero aqui quase
        # sempre e nome de taxon diferente entre locos (acesso em vez de especie).
        "n_taxa_multiloco": multiloco,
    }


def cmd_resumo(a: argparse.Namespace) -> int:
    params = json.loads(Path(a.params).read_text(encoding="utf-8")) if a.params else {}
    versoes = ler_versoes(Path(a.versoes) if a.versoes else None)
    if a.versao_nextflow:
        versoes.setdefault("nextflow", a.versao_nextflow)

    ordem = ordem_rotulos(params)
    entrada = {}
    for p in a.preparos or []:
        d = json.loads(Path(p).read_text(encoding="utf-8"))
        entrada[d["loco"]] = d

    # O nome do loco nunca tem ponto (o workflow saneia), entao o que vem antes
    # do primeiro ponto e o loco em qualquer etapa: X.fas, X.limpo.fa, X.entrada.fa.
    alinhados = {Path(p).name.split(".")[0]: Path(p) for p in (a.alinhados or [])}
    limpos = {Path(p).name.split(".")[0]: Path(p) for p in (a.limpos or [])}

    linhas: list[dict[str, object]] = []
    for rel in sorted(a.relatorios_genes or []):
        loco = Path(rel).name.removesuffix(".iqtree")
        info: dict[str, object] = {"loco": loco}
        info.update(ler_relatorio_iqtree(Path(rel)))
        arv = Path(rel).with_suffix(".treefile")
        if arv.is_file():
            info.update(estat_suporte(arv.read_text(encoding="utf-8"), ordem))
        antes = _comprimento_alinhamento(alinhados[loco]) if loco in alinhados else None
        depois = _comprimento_alinhamento(limpos[loco]) if loco in limpos else None
        if antes is not None:
            info["comprimento_alinhado"] = antes
        if depois is not None:
            info["comprimento_limpo"] = depois
        if antes and depois is not None:
            info["frac_retida"] = round(depois / antes, 3)
        if loco in entrada:
            info["renomeados"] = len(entrada[loco].get("renomeados", {}))
        linhas.append(info)

    fora = sorted(n for n, d in entrada.items() if d.get("status") != "ok")
    resumo: dict[str, object] = {
        "n_locos_entrada": len(entrada),
        "n_locos_inferidos": len(linhas),
        "locos_fora": {n: entrada[n]["status"] for n in fora},
        # O PREPARA_LOCO altera dados (tira '*' interno, descarta registro vazio)
        # e so avisava no stderr da tarefa, que o Nextflow nao mostra quando a
        # tarefa da certo. Aqui o aviso vira dado — e pendencia no metodos.md.
        "avisos_preparo": {n: d["avisos"] for n, d in sorted(entrada.items()) if d.get("avisos")},
        # Outgroup pedido e ausente de todos os locos: arvores sem raiz definida e
        # figuras no ponto medio, sem que nada falhe.
        "outgroup_ausente": bool(params.get("outgroup"))
        and bool(entrada)
        and not any(d.get("outgroup_presentes") for d in entrada.values()),
        "rotulos_suporte": ordem,
        # Mais de 20% removido: Tan et al. (2015) mostram que limpeza automatica
        # agressiva piora arvore de gene unico; vira pendencia no metodos.md.
        "locos_limpeza_agressiva": sorted(
            str(x["loco"]) for x in linhas if float(str(x.get("frac_retida", 1))) < 0.8
        ),
        "locos": linhas,
    }

    dir_concat = Path(a.dir_concat) if a.dir_concat else None
    if rel_concat := _primeiro(dir_concat, "*.iqtree"):
        concat = ler_relatorio_iqtree(rel_concat)
        if arv_concat := _primeiro(dir_concat, "*.treefile"):
            concat.update(estat_suporte(arv_concat.read_text(encoding="utf-8"), ordem))
        if particoes := _primeiro(dir_concat, "particoes.nex"):
            concat["n_particoes_entrada"] = _contar_charsets(particoes)
        if esquema := _primeiro(dir_concat, "*.best_scheme.nex"):
            concat["n_particoes_finais"] = _contar_charsets(esquema)
        if ocupacao := _primeiro(dir_concat, "ocupacao.tsv"):
            concat.update(ler_ocupacao(ocupacao))
        # O -T EFETIVO, nao o pedido: o resourceLimits corta --cpus_concatenado
        # acima do teto sem avisar, e o Metodos tem de dizer o que rodou.
        if log_concat := _primeiro(dir_concat, "*.log"):
            if m := re.search(
                r"^Command:.*?\s-T\s+(\d+)",
                log_concat.read_text(encoding="utf-8", errors="replace"),
                flags=re.MULTILINE,
            ):
                concat["threads"] = int(m[1])
        resumo["concatenado"] = concat

    cf: dict[str, object] = {}
    for p in a.cf or []:
        cf.update(ler_cf(Path(p)))
    if cf:
        resumo["concordancia"] = cf

    dir_especies = Path(a.dir_especies) if a.dir_especies else None
    if arv_esp := _primeiro(dir_especies, "*.tre"):
        especies: dict[str, object] = {"newick": arv_esp.read_text(encoding="utf-8").strip()}
        if log_esp := _primeiro(dir_especies, "*.log"):
            log = log_esp.read_text(encoding="utf-8", errors="replace")
            if m := re.findall(r"^Score:\s*(-?[\d.eE+-]+)", log, flags=re.MULTILINE):
                especies["score_quartetos"] = float(m[-1])
        resumo["arvore_especies"] = especies

    resumo["versoes"] = versoes
    resumo["parametros"] = params

    prefixo = a.saida_prefixo
    Path(f"{prefixo}.json").write_text(
        json.dumps(resumo, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    colunas = [
        "loco",
        "n_seqs",
        "comprimento_alinhado",
        "comprimento_limpo",
        "frac_retida",
        "n_sitios",
        "informativos",
        "padroes",
        "modelo",
        "lnL",
        "comprimento_arvore",
        "ufboot_medio",
        "frac_ufboot_95",
        "frac_sh80_ufb95",
        "renomeados",
    ]
    with Path(f"{prefixo}.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=colunas, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        for linha in linhas:
            w.writerow(linha)
    Path("metodos.md").write_text(
        montar_metodos(params, versoes, resumo), encoding="utf-8", newline="\n"
    )
    return 0


# --------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepara", help="valida e saneia o FASTA de um loco")
    p.add_argument("--loco", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--prefixo", required=True)
    p.add_argument("--min-taxa", type=int, default=4)
    p.add_argument("--outgroup", default="")
    p.add_argument("--tipo", default="", choices=["", "DNA", "AA", "dna", "aa"])
    p.add_argument("--nome-taxon", default="primeiro", choices=["primeiro", "completo", "colchete"])
    p.add_argument("--desalinhar", action="store_true", help="remove gaps para realinhar")
    p.set_defaults(func=cmd_prepara)

    c = sub.add_parser("concatena", help="supermatriz + particoes NEXUS")
    c.add_argument("--locos", nargs="+", required=True)
    c.add_argument("--alinhamentos", nargs="+", required=True)
    c.add_argument("--outgroup", default="")
    c.add_argument("--saida-fasta", default="supermatriz.fa")
    c.add_argument("--saida-particoes", default="particoes.nex")
    c.add_argument("--saida-ocupacao", default="ocupacao.tsv")
    c.add_argument("--saida-env", default="concatena.env")
    c.set_defaults(func=cmd_concatena)

    r = sub.add_parser("resumo", help="tabela, JSON e paragrafo de Metodos")
    r.add_argument("--params")
    r.add_argument("--versoes")
    r.add_argument("--versao-nextflow", default="")
    r.add_argument("--preparos", nargs="*")
    r.add_argument("--alinhados", nargs="*")
    r.add_argument("--limpos", nargs="*")
    r.add_argument("--relatorios-genes", nargs="*")
    r.add_argument(
        "--dir-concat", help="diretorio com .iqtree, .treefile e .best_scheme.nex da supermatriz"
    )
    r.add_argument("--cf", nargs="*")
    r.add_argument("--dir-especies", help="diretorio com a arvore (.tre) e o log do wASTRAL")
    r.add_argument("--saida-prefixo", default="resumo")
    r.set_defaults(func=cmd_resumo)

    a = ap.parse_args(argv)
    try:
        return a.func(a)
    except ErroEntrada as exc:
        print(f"ERRO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
