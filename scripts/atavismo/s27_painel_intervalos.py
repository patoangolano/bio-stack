"""O resto do painel: AMELX, MMP20, DSPP e ACP4, cada um no seu endereco.

O s26 mediu ENAM, AMBN e AMTN dentro do intervalo NRG1..LPL. Os outros genes do
painel nao moram ali — no jacare estao em tres outros cromossomos:

    AMELX  chr1  NC_081824.1   dentro/junto de ARHGAP6 nos mamiferos
    MMP20  chr1  NC_081824.1   cluster de metaloproteinases
    DSPP   chr2  NC_081825.1   cluster SCPP/SIBLING (dentina)
    ACP4   chr5  NC_081828.1

Rodar o s26 com eles daria zero em tudo, e zero que nao significa nada. Entao
cada gene ganha o seu proprio intervalo sintenico:

    ANCORAS   o vizinho com simbolo oficial mais proximo de cada lado do alvo no
              jacare, exigido anotado nas especies sem dente e no Anolis. Sao o
              controle POSITIVO de cada intervalo: se a ancora nao alinhar, o
              intervalo esta errado ou o metodo nao tem sensibilidade ali.
    COMPETIDORES  toda proteina anotada no intervalo do jacare entre as ancoras
              (inclusive os LOC). O MMP20 vive cercado de colagenases parálogas;
              sem competicao, o tblastn do MMP20 "acharia" MMP13 nas aves. Um HSP
              so conta para o alvo se nenhum competidor alinhar melhor no mesmo
              trecho (s27_indice.py).
    NEGATIVO  LPL do jacare (chr3), o mesmo controle positivo do s24/s26 — aqui
              fora de qualquer intervalo do painel.

ODAM e ODAPH ficam de fora por um motivo medido, nao por pressa: nao ha um unico
gene anotado com esses simbolos em Sauropsida no NCBI Gene (busca de
2026-10-17). Consultar com a proteina humana, a ~320 Myr, daria negativo
ininterpretavel. O passo deles e outro: localizar primeiro no proprio jacare.

Este script so BAIXA e escreve FASTA/JSON. O alinhamento roda no WSL:

    uv run python scripts/atavismo/s27_painel_intervalos.py
    wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s27_alinhar.sh
    uv run python scripts/atavismo/s27_indice.py
"""

from __future__ import annotations

import json
import os
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import requests

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
CM: dict[str, str] = {"tool": "bio-stack", "email": "limatheus458@gmail.com"}
if os.environ.get("NCBI_API_KEY"):
    CM["api_key"] = os.environ["NCBI_API_KEY"]
PAUSA = 0.15 if "api_key" in CM else 0.4

RAIZ = Path(__file__).resolve().parents[2]
OUT = RAIZ / "runs" / "atavismo" / "painel"
LPL_FAA = RAIZ / "runs" / "atavismo" / "lncrna" / "jacare_proteinas.faa"

S = requests.Session()
S.headers["User-Agent"] = "bio-stack/s27"

# taxid, rotulo, nome, estado dentario — mesmo conjunto do s26
ESPECIES = [
    (8496, "Amis", "Alligator mississippiensis", "COM dente"),
    (94835, "Ggan", "Gavialis gangeticus", "COM dente"),
    (28377, "Acar", "Anolis carolinensis", "COM dente"),
    (8478, "Cpic", "Chrysemys picta", "sem dente"),
    (9031, "Ggal", "Gallus gallus", "sem dente"),
    (60713, "Mund", "Melopsittacus undulatus", "sem dente"),
    (2489341, "Shab", "Strigops habroptila", "sem dente"),
]
# Ancora tem de existir onde a pergunta e decidida: nas sem dente e no Anolis
# (o controle distante). Gavialis tem o painel em scaffolds nao posicionados e
# nao veta ancora. Melopsittacus tambem nao: a API de genes nao devolve
# localizacao para nenhum vizinho dele (medido em 2026-10-17; o NRG1 ja faltava
# no s26), entao exigi-lo vetaria toda ancora.
EXIGIDAS = {"Acar", "Cpic", "Ggal", "Shab"}

ALVOS = ["AMELX", "MMP20", "DSPP", "ACP4"]
JANELA_FT = 1_500_000  # bp de cada lado do alvo na tabela de features do jacare
# 25 e nao 8: a vizinhanca do ACP4 (equivalente ao 19q13 humano) e rica em GC,
# e nas aves boa parte desses genes falta no assembly. Com 8 candidatas por lado
# nenhuma estava no frango E no kakapo.
MAX_CANDIDATAS = 25
FOLGA = 60_000  # bp alem das ancoras em cada especie
TETO_INTERVALO = 3_000_000  # acima disso a "sintenia" e outra coisa
# Quando o par falha numa especie (rearranjo, sequencias diferentes, ancora
# ausente), janela em volta da ancora que existir, a mais proxima do alvo no
# jacare. No Anolis o DSPP fica a 179 Mb do NUDT9, mas colado no IBSP.
JANELA_UMA = 750_000


@dataclass
class Gene:
    sym: str
    desc: str
    ini: int
    fim: int
    proteina: str | None = None


@dataclass
class Intervalo:
    rotulo: str
    situacao: str
    acesso: str | None = None
    sequencia: str | None = None
    janela_ini: int | None = None
    janela_fim: int | None = None
    bp: int | None = None
    arquivo: str | None = None
    ancoras: list[str] = field(default_factory=list)  # as que servem de controle aqui


@dataclass
class Painel:
    alvo: str
    jacare: dict
    ancora_esq: str | None = None
    ancora_dir: str | None = None
    competidores: list[str] = field(default_factory=list)
    intervalos: list[Intervalo] = field(default_factory=list)


def gene_pos(sym: str, taxon: int) -> dict | None:
    """Primeira localizacao genomica anotada do gene; None se nao anotado."""
    try:
        d = S.get(f"{API}/gene/symbol/{sym}/taxon/{taxon}",
                  params={"returned_content": "COMPLETE"}, timeout=90).json()
    except (requests.RequestException, json.JSONDecodeError):
        return None
    finally:
        time.sleep(PAUSA)
    for rep in d.get("reports", []):
        g = rep.get("gene", {})
        if g.get("symbol", "").upper() != sym.upper():
            continue
        for a in g.get("annotations", []):
            for gl in a.get("genomic_locations", []):
                gr = gl.get("genomic_range", {})
                return {"acc": gl.get("genomic_accession_version"), "seq": gl.get("sequence_name"),
                        "ini": int(gr.get("begin", 0)), "fim": int(gr.get("end", 0))}
    return None


def tabela_features(acesso: str, ini: int, fim: int) -> list[Gene]:
    """Genes e a primeira proteina de cada, de uma fatia do cromossomo.

    O efetch rettype=ft devolve coordenadas RELATIVAS ao recorte; aqui elas
    voltam a ser absolutas.
    """
    r = S.get(f"{EU}/efetch.fcgi",
              params={**CM, "db": "nuccore", "id": acesso, "rettype": "ft", "retmode": "text",
                      "seq_start": ini, "seq_stop": fim}, timeout=300)
    r.raise_for_status()
    time.sleep(PAUSA)
    genes: dict[str, Gene] = {}
    tipo, coords, quals = None, (0, 0), {}

    def fecha() -> None:
        if tipo == "gene" and "gene" in quals:
            a, b = sorted(coords)
            genes.setdefault(quals["gene"], Gene(quals["gene"], quals.get("gene_desc", ""),
                                                 ini - 1 + a, ini - 1 + b))
        elif tipo == "CDS" and "gene" in quals and "protein_id" in quals:
            g = genes.get(quals["gene"])
            if g and not g.proteina:
                g.proteina = quals["protein_id"].strip("|").split("|")[-1]

    for linha in r.text.splitlines():
        if m := re.match(r"^<?(\d+)\t>?(\d+)\t(\w+)$", linha):
            fecha()
            tipo, coords, quals = m.group(3), (int(m.group(1)), int(m.group(2))), {}
        elif linha.startswith("\t\t\t"):
            partes = linha.strip("\t").split("\t", 1)
            if len(partes) == 2:
                quals.setdefault(partes[0], partes[1])
    fecha()
    return sorted(genes.values(), key=lambda g: g.ini)


def escolhe_ancora(candidatas: list[Gene], cache: dict) -> str | None:
    """Primeiro vizinho nomeado presente em todas as especies EXIGIDAS."""
    melhor, melhor_n = None, -1
    for g in candidatas[:MAX_CANDIDATAS]:
        pos = {rot: cache.setdefault((g.sym, rot), gene_pos(g.sym, tx)) for tx, rot, *_ in ESPECIES}
        presentes = {rot for rot, p in pos.items() if p}
        print(f"      candidata {g.sym:10s} anotada em {len(presentes)}/7 "
              f"(faltam: {', '.join(sorted(EXIGIDAS - presentes)) or '-'})")
        if EXIGIDAS <= presentes:
            return g.sym
        if len(presentes & EXIGIDAS) > melhor_n:
            melhor, melhor_n = g.sym, len(presentes & EXIGIDAS)
    return melhor


def fatia(acesso: str, ini: int, fim: int, rotulo: str) -> str | None:
    r = S.get(f"{EU}/efetch.fcgi",
              params={**CM, "db": "nuccore", "id": acesso, "rettype": "fasta", "retmode": "text",
                      "seq_start": max(1, ini), "seq_stop": fim}, timeout=600)
    time.sleep(PAUSA)
    if r.status_code != 200 or not r.text.strip().startswith(">"):
        return None
    corpo = r.text.strip().split("\n", 1)[1]
    return f">{rotulo} {acesso}:{max(1, ini)}-{fim}\n{corpo}\n"


def proteinas(acessos: dict[str, str]) -> str:
    """FASTA das proteinas do jacare, cabecalho >SIMBOLO_Amis acesso."""
    r = S.get(f"{EU}/efetch.fcgi",
              params={**CM, "db": "protein", "id": ",".join(acessos.values()),
                      "rettype": "fasta", "retmode": "text"}, timeout=300)
    r.raise_for_status()
    time.sleep(PAUSA)
    por_acesso = {v: k for k, v in acessos.items()}
    saida = []
    for bloco in r.text.strip().split(">")[1:]:
        cab, corpo = bloco.split("\n", 1)
        acc = cab.split(" ", 1)[0]
        sym = por_acesso.get(acc)
        if sym:
            saida.append(f">{sym}_Amis {acc}\n{corpo.strip()}\n")
    return "".join(saida)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cache: dict = {}
    paineis: list[Painel] = []
    lpl = LPL_FAA.read_text(encoding="utf-8").split(">LPL_Amis", 1)[1].split(">", 1)[0]
    lpl_bloco = ">LPL_Amis" + lpl

    for alvo in ALVOS:
        print(f"\n=== {alvo} ===")
        pj = gene_pos(alvo, 8496)
        if not pj:
            print("  nao anotado no jacare — pulado")
            continue
        genes = tabela_features(pj["acc"], max(1, pj["ini"] - JANELA_FT), pj["fim"] + JANELA_FT)
        nomeados = [g for g in genes if g.proteina and not g.sym.startswith("LOC") and g.sym != alvo]
        a_ini, a_fim = int(pj["ini"]), int(pj["fim"])
        esq = sorted((g for g in nomeados if g.fim < a_ini), key=lambda g: a_ini - g.fim)
        dir_ = sorted((g for g in nomeados if g.ini > a_fim), key=lambda g: g.ini - a_fim)
        print(f"  jacare {pj['acc']}:{pj['ini']:,}-{pj['fim']:,}  |  {len(genes)} genes na janela, "
              f"{len(esq)} nomeados a esquerda, {len(dir_)} a direita")
        print("  ancora esquerda:")
        a_esq = escolhe_ancora(esq, cache)
        print("  ancora direita:")
        a_dir = escolhe_ancora(dir_, cache)
        painel = Painel(alvo, pj, a_esq, a_dir)
        if not (a_esq and a_dir):
            print("  sem par de ancoras — pulado")
            paineis.append(painel)
            continue

        # Consultas: alvo, ancoras e todo competidor do intervalo do jacare.
        gi = {g.sym: g for g in genes}
        lo = min(gi[a_esq].ini, gi[a_dir].ini)
        hi = max(gi[a_esq].fim, gi[a_dir].fim)
        consulta = {g.sym: g.proteina for g in genes if g.proteina and g.ini <= hi and g.fim >= lo}
        painel.competidores = sorted(s for s in consulta if s not in (alvo, a_esq, a_dir))
        faa = proteinas(consulta) + lpl_bloco
        (OUT / f"painel_{alvo}_proteinas.faa").write_text(faa, encoding="utf-8")
        print(f"  consultas: alvo + 2 ancoras ({a_esq}, {a_dir}) + "
              f"{len(painel.competidores)} competidores + LPL (negativo)")

        dist = {s: min(abs(gi[s].fim - a_ini), abs(gi[s].ini - a_fim)) for s in (a_esq, a_dir)}
        por_proximidade = sorted(dist, key=dist.__getitem__)
        for _tx, rot, _nome, _estado in ESPECIES:
            pe, pd = cache.get((a_esq, rot)), cache.get((a_dir, rot))
            janela: tuple[str, str, int, int] | None = None
            ancoras_ok: list[str] = []
            if pe and pd and pe["acc"] == pd["acc"]:
                ini = min(pe["ini"], pe["fim"], pd["ini"], pd["fim"]) - FOLGA
                fim = max(pe["ini"], pe["fim"], pd["ini"], pd["fim"]) + FOLGA
                if fim - ini <= TETO_INTERVALO:
                    janela, ancoras_ok, nota = (pe["acc"], pe["seq"], ini, fim), [a_esq, a_dir], "ok"
                else:
                    nota = f"ancoras a {fim - ini:,} bp"
            elif pe and pd:
                nota = f"ancoras em sequencias diferentes ({pe['seq']} / {pd['seq']})"
            else:
                nota = f"sem ancora {a_esq if not pe else a_dir}"
            if janela is None:
                for sym in por_proximidade:
                    p = cache.get((sym, rot))
                    if p:
                        meio = (p["ini"] + p["fim"]) // 2
                        janela = (p["acc"], p["seq"], meio - JANELA_UMA, meio + JANELA_UMA)
                        ancoras_ok, nota = [sym], f"uma ancora ({sym}); {nota}"
                        break
            if janela is None:
                painel.intervalos.append(Intervalo(rot, nota))
                print(f"    {rot}: {nota} — pulado")
                continue
            acc, seq, ini, fim = janela
            arq = OUT / f"intervalo_{alvo}_{rot}.fna"
            fasta = fatia(acc, ini, fim, f"{rot}_{alvo}")
            if not fasta:
                painel.intervalos.append(Intervalo(rot, "efetch falhou", acc, seq))
                continue
            arq.write_text(fasta, encoding="utf-8")
            bp = len(fasta.split("\n", 1)[1].replace("\n", ""))
            painel.intervalos.append(
                Intervalo(rot, nota, acc, seq, max(1, ini), fim, bp, arq.name, ancoras_ok))
            print(f"    {rot}: seq {seq:>6}  {bp:>9,} bp  {nota}")
        paineis.append(painel)

    (OUT / "s27_resumo.json").write_text(
        json.dumps([asdict(p) for p in paineis], indent=2, ensure_ascii=False), encoding="utf-8")
    prontos = sum(1 for p in paineis for i in p.intervalos if i.situacao == "ok")
    print(f"\n{prontos} intervalos extraidos. Agora, no WSL:")
    print("  wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s27_alinhar.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
