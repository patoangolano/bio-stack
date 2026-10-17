"""ODAM e ODAPH (= C4orf26): localizar no jacare, depois medir nas outras especies.

Os dois genes ficaram fora do s27 por um motivo medido: nao ha gene anotado com
esses simbolos em Sauropsida no NCBI Gene (busca de 2026-10-17). Ausencia de
anotacao e o nivel 1 da escada — fraco. Aqui o caminho e o inverso do s26/s27:
primeiro achar o gene no PROPRIO jacare (com dente, onde ele tem de estar
funcional), depois usar a proteina do jacare como consulta nas demais.

ODAM   No humano fica no cluster SCPP de 4q13, ao lado de FDCSP/CSN3. Nos
       sauropsidas o cluster SCPP (ENAM-AMBN-AMTN) mora no intervalo NRG1..LPL do
       s26, e e la que se procura. Consultas: ODAM humano e de dois anfibios
       (Xenopus, Rhinatrema) — mamifero e anfibio de lados opostos de
       Sauropsida, para nao depender de uma so distancia.
ODAPH  Fora dos clusters SCPP (Liang et al. 2021: introns de fase 1, nao e
       SCPP). No humano: RCHY1-THAP6-ODAPH-CDKL2-G3BP2 em 4q21.1. A janela de
       cada especie sai das ancoras CDKL2/THAP6/RCHY1/G3BP2 anotadas (gene_pos do
       s27). Consultas: ODAPH humano e do gamba (Monodelphis).

Controles (o "mesmo limiar" da escada, nivel 2):
  ancoras    ENAM/AMBN/AMTN do jacare no intervalo SCPP; CDKL2/THAP6/RCHY1/G3BP2
             do jacare na janela do ODAPH. Se a ancora nao alinha, a janela esta
             errada ou o metodo nao tem sensibilidade ali — negativo nao vale.
  com dente  Crocodylus, Gavialis, Anolis, Python, Podarcis: a proteina do jacare
             TEM de achar o gene nelas no mesmo limiar, e sem lesao (controle
             negativo de lesao, como no s29).

Ordem:
    uv run python scripts/atavismo/s32_odam_odaph.py baixar
    export MSYS_NO_PATHCONV=1
    wsl.exe -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s32_odam_odaph.sh etapa1
    uv run python scripts/atavismo/s32_odam_odaph.py modelo
    wsl.exe -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s32_odam_odaph.sh etapa2
    uv run python scripts/atavismo/s32_odam_odaph.py tabela
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
OUT = RAIZ / "runs" / "atavismo" / "odam_odaph"
LNC = RAIZ / "runs" / "atavismo" / "lncrna"
PROTEOMA = RAIZ / "runs" / "atavismo" / "painel" / "proteoma_Amis" / "protein.faa"


def _importa(nome: str, arquivo: str):
    spec = importlib.util.spec_from_file_location(nome, AQUI / arquivo)
    if spec is None or spec.loader is None:
        raise SystemExit(f"nao carreguei {arquivo}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[nome] = mod  # dataclass precisa do modulo registrado
    spec.loader.exec_module(mod)
    return mod


s27 = _importa("s27", "s27_painel_intervalos.py")
s29 = _importa("s29", "s29_lesoes.py")

# rotulo, fonte, acesso
CONSULTAS = [
    ("ODAM_Hsap", "uniprot", "A1E959"),
    ("ODAM_Xtro", "ncbi", "XP_012813187"),
    ("ODAM_Rbiv", "ncbi", "XP_029453991"),
    ("ODAPH_Hsap", "uniprot", "Q17RF5"),
    ("ODAPH_Mdom", "ncbi", "XP_001375528"),
]
# taxid, rotulo, nome, dente, grupo — as 18 especies com intervalo NRG1..LPL
ESPECIES = [
    (8496, "Amis", "Alligator mississippiensis", "COM", "Crocodylia"),
    (8502, "Cpor", "Crocodylus porosus", "COM", "Crocodylia"),
    (94835, "Ggan", "Gavialis gangeticus", "COM", "Crocodylia"),
    (28377, "Acar", "Anolis carolinensis", "COM", "Squamata"),
    (176946, "Pbiv", "Python bivittatus", "COM", "Squamata"),
    (64176, "Pmur", "Podarcis muralis", "COM", "Squamata"),
    (8478, "Cpic", "Chrysemys picta", "sem", "Testudines"),
    (8469, "Cmyd", "Chelonia mydas", "sem", "Testudines"),
    (1825980, "Gevg", "Gopherus evgoodei", "sem", "Testudines"),
    (13735, "Psin", "Pelodiscus sinensis", "sem", "Testudines"),
    (9031, "Ggal", "Gallus gallus", "sem", "Galloanserae"),
    (8839, "Apla", "Anas platyrhynchos", "sem", "Galloanserae"),
    (2489341, "Shab", "Strigops habroptila", "sem", "Neoaves"),
    (59729, "Tgut", "Taeniopygia guttata", "sem", "Neoaves"),
    (8790, "Dnov", "Dromaius novaehollandiae", "sem", "Palaeognathae"),
    (8801, "Scam", "Struthio camelus", "sem", "Palaeognathae"),
    (308060, "Arow", "Apteryx rowi", "sem", "Palaeognathae"),
    (94827, "Tgua", "Tinamus guttatus", "sem", "Palaeognathae"),
]
# CDKL2 primeiro: no humano e no jacare o ODAPH esta colado nele.
ANCORAS_ODAPH = ["CDKL2", "THAP6", "RCHY1", "G3BP2"]
VIZINHANCA = 1_000_000  # ancora a mais que isso do centro e outro bloco (frango: THAP6 a 12 Mb)
FOLGA = 100_000  # o ODAPH pode estar do lado de fora se o bloco inverteu
JANELA_UMA = 300_000
THAP6_JANELA = 150_000
CTRL_XTRO = 650_000  # meia janela: 1,3 Mb no total, o tamanho do intervalo do jacare
CONTROLES_ODAM = ["ENAM_Amis", "AMBN_Amis", "AMTN_Amis"]
# Limiar unico de "achado" para alvo e controle (tblastn/blastp). Num intervalo
# de ~0,3-1,3 Mb, E <= 1e-3 separa os hits de ancora do fundo: o fundo de
# consultas aleatorias fica em E 0,01-10 (medido no e1_*: melhor ODAPH humano
# fora do gene, E >= 0,1).
LIMIAR_E = 1e-3


def fasta(texto: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    for bloco in texto.strip().split(">")[1:]:
        cab, corpo = bloco.split("\n", 1)
        seqs[cab.split()[0]] = re.sub(r"\s", "", corpo)
    return seqs


def escreve_fasta(caminho: Path, seqs: dict[str, str]) -> None:
    with caminho.open("w", encoding="utf-8", newline="\n") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n" + "\n".join(v[i:i + 60] for i in range(0, len(v), 60)) + "\n")


def proteina(fonte: str, acesso: str) -> str:
    if fonte == "uniprot":
        r = s27.S.get(f"https://rest.uniprot.org/uniprotkb/{acesso}.fasta", timeout=60)
    else:
        r = s27.S.get(f"{s27.EU}/efetch.fcgi", timeout=60,
                      params={**s27.CM, "db": "protein", "id": acesso, "rettype": "fasta", "retmode": "text"})
    r.raise_for_status()
    return next(iter(fasta(r.text).values()))


def janela_odaph(taxon: int) -> dict | None:
    pos = {g: s27.gene_pos(g, taxon) for g in ANCORAS_ODAPH}
    centro = next((g for g in ANCORAS_ODAPH if pos[g]), None)
    if centro is None:
        return None
    c = pos[centro]
    meio = (c["ini"] + c["fim"]) // 2
    usadas = {g: p for g, p in pos.items()
              if p and p["acc"] == c["acc"] and abs((p["ini"] + p["fim"]) // 2 - meio) <= VIZINHANCA}
    if len(usadas) > 1:
        lo = min(min(p["ini"], p["fim"]) for p in usadas.values()) - FOLGA
        hi = max(max(p["ini"], p["fim"]) for p in usadas.values()) + FOLGA
    else:
        lo, hi = meio - JANELA_UMA, meio + JANELA_UMA
    return {"acesso": c["acc"], "sequencia": c["seq"], "ini": max(1, lo), "fim": hi,
            "ancoras": {g: [p["ini"], p["fim"]] for g, p in usadas.items()},
            "ancoras_fora": {g: [p["acc"], p["ini"], p["fim"]] for g, p in pos.items()
                             if p and g not in usadas}}


def baixar() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    # 1. consultas
    consultas = {rot: proteina(f, a) for rot, f, a in CONSULTAS}
    escreve_fasta(OUT / "consultas_odam.faa", {k: v for k, v in consultas.items() if k.startswith("ODAM")})
    escreve_fasta(OUT / "consultas_odaph.faa", {k: v for k, v in consultas.items() if k.startswith("ODAPH")})
    for k, v in consultas.items():
        print(f"  consulta {k:11s} {len(v):4d} aa")

    # 2. controles do intervalo SCPP: as mesmas proteinas do jacare do s26
    jac = fasta((LNC / "jacare_proteinas.faa").read_text(encoding="utf-8"))
    escreve_fasta(OUT / "controles_odam.faa", {k: jac[k] for k in CONTROLES_ODAM})

    # 3. janela do ODAPH por especie + tabela de genes anotados nela
    janelas, genes_tsv = {}, []
    for taxon, rot, nome, dente, grupo in ESPECIES:
        j = janela_odaph(taxon)
        if j is None:
            print(f"  {rot} {nome:28s} nenhuma ancora anotada")
            janelas[rot] = {"situacao": "sem ancora"}
            continue
        arq = OUT / f"janela_odaph_{rot}.fna"
        if not arq.exists():
            fa = s27.fatia(j["acesso"], j["ini"], j["fim"], f"{rot}_ODAPH")
            if fa is None:
                janelas[rot] = {"situacao": "efetch falhou", **j}
                continue
            arq.write_text(fa, encoding="utf-8")
        try:
            genes = s27.tabela_features(j["acesso"], j["ini"], j["fim"])
        except Exception as e:  # noqa: BLE001 — tabela e informativa, a janela ja esta no disco
            print(f"    tabela de features falhou: {e}")
            genes = []
        for g in genes:
            genes_tsv.append({"especie": rot, "simbolo": g.sym, "ini": g.ini, "fim": g.fim,
                              "proteina": g.proteina or "", "descricao": g.desc})
        janelas[rot] = {"situacao": "ok", "nome": nome, "dente": dente, "grupo": grupo, **j,
                        "bp": j["fim"] - j["ini"] + 1}
        print(f"  {rot} {nome:28s} {j['acesso']}:{j['ini']}-{j['fim']} ({j['fim'] - j['ini'] + 1:,} bp) "
              f"ancoras {','.join(j['ancoras'])}  genes anotados {len(genes)}")
        # Nas aves CDKL2 encosta direto em RCHY1: o trecho ODAPH-THAP6 do jacare
        # saiu dali. Onde o THAP6 esta anotado em outro lugar (frango: chr4
        # 45,9 Mb, entre LIN54 e SEC31A), o ODAPH pode ter ido junto — segunda
        # janela em volta dele.
        if "THAP6" in j["ancoras_fora"]:
            acc2, a, b = j["ancoras_fora"]["THAP6"]
            lo2, hi2 = max(1, (a + b) // 2 - THAP6_JANELA), (a + b) // 2 + THAP6_JANELA
            arq2 = OUT / f"janela_odaph_{rot}_THAP6.fna"
            if not arq2.exists():
                fa = s27.fatia(acc2, lo2, hi2, f"{rot}_THAP6")
                if fa:
                    arq2.write_text(fa, encoding="utf-8")
            for g in s27.tabela_features(acc2, lo2, hi2):
                genes_tsv.append({"especie": f"{rot}_THAP6", "simbolo": g.sym, "ini": g.ini, "fim": g.fim,
                                  "proteina": g.proteina or "", "descricao": g.desc})
            janelas[f"{rot}_THAP6"] = {"situacao": "ok" if arq2.exists() else "efetch falhou",
                                       "nome": nome, "dente": dente, "grupo": grupo, "acesso": acc2,
                                       "ini": lo2, "fim": hi2, "bp": hi2 - lo2 + 1, "ancoras": {"THAP6": [a, b]}}
            print(f"  {rot}_THAP6 janela secundaria {acc2}:{lo2}-{hi2}")
    (OUT / "s32_janelas_odaph.json").write_text(json.dumps(janelas, indent=2, ensure_ascii=False),
                                                encoding="utf-8")
    with (OUT / "s32_janela_genes.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["especie", "simbolo", "ini", "fim", "proteina", "descricao"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(genes_tsv)

    # 4. controles da janela do ODAPH: as ancoras do jacare
    gj = {g["simbolo"]: g["proteina"] for g in genes_tsv if g["especie"] == "Amis"}
    acc = {g: gj[g] for g in ANCORAS_ODAPH if gj.get(g)}
    ctrl = fasta(s27.proteinas(acc))
    escreve_fasta(OUT / "controles_odaph.faa", ctrl)
    print(f"  controles ODAPH: {', '.join(ctrl)}")

    # 5. Controle positivo de SENSIBILIDADE para o ODAM: a regiao genomica do
    # ODAM de Xenopus (anotado), do mesmo tamanho do intervalo do jacare
    # (~1,3 Mb, para o E-value ser comparavel). Humano x Xenopus ~350 Ma, mais
    # longe que humano x jacare (~320 Ma): se a consulta humana acha o ODAM do
    # sapo no genoma, nao achar no jacare deixa de ser falta de sensibilidade.
    px = s27.gene_pos("ODAM", 8364)
    arq = OUT / "controle_odam_Xtro.fna"
    if px and not arq.exists():
        meio = (px["ini"] + px["fim"]) // 2
        fa = s27.fatia(px["acc"], meio - CTRL_XTRO, meio + CTRL_XTRO, "Xtro_ODAM")
        if fa:
            arq.write_text(fa, encoding="utf-8")
    if px:
        janelas["controle_Xtro"] = {"acesso": px["acc"], "odam": [px["ini"], px["fim"]],
                                    "ini": (px["ini"] + px["fim"]) // 2 - CTRL_XTRO,
                                    "fim": (px["ini"] + px["fim"]) // 2 + CTRL_XTRO}
        print(f"  controle Xenopus: ODAM {px['acc']}:{px['ini']}-{px['fim']}")

    # 6. HMM do Pfam para ODAM (PF15424). Perfil e mais sensivel que uma
    # consulta so para proteina rica em P/Q de exons curtos, onde o tblastn
    # perde ate o ODAM do Xenopus.
    hmm = OUT / "pf15424_odam.hmm"
    if not hmm.exists():
        r = s27.S.get("https://www.ebi.ac.uk/interpro/api/entry/pfam/PF15424",
                      params={"annotation": "hmm"}, timeout=120)
        r.raise_for_status()
        hmm.write_bytes(gzip.decompress(r.content))
    print(f"  {hmm.read_text(encoding='utf-8').splitlines()[2]}")

    # 7. Swiss-Prot humano para o segundo lado do melhor hit reciproco
    sp = OUT / "sprot_humano.faa"
    if not sp.exists():
        r = s27.S.get("https://rest.uniprot.org/uniprotkb/stream", timeout=600,
                      params={"query": "organism_id:9606 AND reviewed:true", "format": "fasta"})
        r.raise_for_status()
        sp.write_text(r.text, encoding="utf-8", newline="\n")
    print(f"  Swiss-Prot humano: {sp.read_text(encoding='utf-8').count('>')} proteinas")
    (OUT / "s32_janelas_odaph.json").write_text(json.dumps(janelas, indent=2, ensure_ascii=False),
                                                encoding="utf-8")
    return vaos()


def _nomeado(g: dict) -> bool:
    """Gene com simbolo oficial e proteina: serve de parede do vao."""
    return bool(g["proteina"]) and not g["simbolo"].upper().startswith(("LOC", "MIR", "TRNA"))


def _distancia(a: dict, b: dict) -> int:
    return max(0, max(a["ini"], b["ini"]) - min(a["fim"], b["fim"]))


def lado_de_fora(cd: dict, gs: list[dict]) -> int:
    """+1 se o ODAPH ancestral fica a direita do CDKL2, -1 a esquerda, 0 sem referencia.

    Primeiro: o lado do vizinho de fora mais proximo (THAP6, RCHY1 ou PARM1) a
    ate 100 kb. Se nenhum (frango: RCHY1 a 530 kb), o lado oposto ao vizinho de
    dentro mais proximo (G3BP2, SH2D4A ou AP1AR — o que encosta no CDKL2 do lado
    do G3BP2 nas aves) a ate 100 kb.
    """
    for grupo, sinal in ((("THAP6", "RCHY1", "PARM1"), 1), (("G3BP2", "SH2D4A", "AP1AR"), -1)):
        perto = [g for g in gs if g["simbolo"].upper() in grupo and _distancia(g, cd) <= 100_000]
        if perto:
            g = min(perto, key=lambda g: _distancia(g, cd))
            return sinal * (1 if g["ini"] > cd["ini"] else -1)
    return 0


def vaos() -> int:
    """O vao do ODAPH: entre CDKL2 e o primeiro gene nomeado do lado de fora.

    Ordem ancestral (humano e jacare): G3BP2-CDKL2-ODAPH-THAP6-RCHY1-PARM1. O
    ODAPH, se existe ali, esta no vao entre CDKL2 e o vizinho nomeado seguinte
    do lado de fora (THAP6 em crocodilo e tartaruga, PARM1 em escamado, RCHY1
    em ave). Genes LOC dentro do vao sao candidatos, nao parede: no Anolis e no
    Podarcis o ODAPH esta anotado como lncRNA. Vao de 0 bp = o vizinho encosta
    ou sobrepoe o CDKL2.
    """
    janelas = json.loads((OUT / "s32_janelas_odaph.json").read_text(encoding="utf-8"))
    genes = list(csv.DictReader((OUT / "s32_janela_genes.tsv").open(encoding="utf-8"), delimiter="\t"))
    resumo: dict[str, dict] = {}
    transcritos: dict[str, str] = {}
    for rot, j in janelas.items():
        if j.get("situacao") != "ok" or rot.endswith("_THAP6") or "nome" not in j:
            continue
        gs: list[dict[str, Any]] = [{**g, "ini": int(g["ini"]), "fim": int(g["fim"])}
                                    for g in genes if g["especie"] == rot]
        cd = next((g for g in gs if g["simbolo"].upper() == "CDKL2"), None)
        if cd is None:
            resumo[rot] = {"situacao": "sem CDKL2"}
            continue
        lado = lado_de_fora(cd, gs)
        if lado == 0:
            resumo[rot] = {"situacao": "sem referencia de lado"}
            continue
        if lado > 0:
            viz = min((g for g in gs if _nomeado(g) and g["ini"] > cd["ini"] and g["fim"] > cd["fim"]),
                      key=lambda g: g["ini"], default=None)
            lo, hi = cd["fim"] + 1, (viz["ini"] - 1) if viz else j["fim"]
        else:
            viz = max((g for g in gs if _nomeado(g) and g["ini"] < cd["ini"] and g["fim"] < cd["fim"]),
                      key=lambda g: g["fim"], default=None)
            lo, hi = (viz["fim"] + 1) if viz else j["ini"], cd["ini"] - 1
        corpo = (OUT / f"janela_odaph_{rot}.fna").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
        seq = corpo[max(0, lo - j["ini"]):max(0, hi - j["ini"] + 1)]
        locs = [g for g in gs if g["simbolo"].upper().startswith("LOC") and g["fim"] >= lo and g["ini"] <= hi]
        resumo[rot] = {"situacao": "ok", "acesso": j["acesso"], "ini": lo, "fim": hi, "bp": len(seq),
                       "N": seq.upper().count("N"), "parede": viz["simbolo"] if viz else None,
                       "LOC_no_vao": [g["simbolo"] for g in locs]}
        if seq:
            escreve_fasta(OUT / f"vao_odaph_{rot}.fna", {f"{rot}_vao {j['acesso']}:{lo}-{hi}": seq})
        # transcritos dos LOC sem proteina no vao (o ODAPH anotado como lncRNA?)
        for g in locs:
            if g["proteina"]:
                continue
            r = s27.S.get(f"{s27.EU}/efetch.fcgi", timeout=120,
                          params={**s27.CM, "db": "nuccore", "id": j["acesso"], "rettype": "ft",
                                  "retmode": "text", "seq_start": g["ini"], "seq_stop": g["fim"]})
            for tid in re.findall(r"transcript_id\tref\|([^|]+)\|", r.text):
                transcritos[tid] = f"{rot}_{g['simbolo']}"
        print(f"  vao {rot:5s} {j['acesso']}:{lo}-{hi} {len(seq):>7,} bp  N={resumo[rot]['N']:<5d} "
              f"parede={resumo[rot]['parede']}  LOC={','.join(resumo[rot]['LOC_no_vao']) or '-'}")
    if transcritos:
        r = s27.S.get(f"{s27.EU}/efetch.fcgi", timeout=120,
                      params={**s27.CM, "db": "nuccore", "id": ",".join(transcritos), "rettype": "fasta",
                              "retmode": "text"})
        r.raise_for_status()
        tx = {f"{transcritos.get(k, 'LOC')}_{k}": v for k, v in fasta(r.text).items()}
        escreve_fasta(OUT / "transcritos_vao.fna", tx)
        print(f"  transcritos de LOC no vao: {', '.join(tx)}")
    (OUT / "s32_vaos.json").write_text(json.dumps(resumo, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


def le_tsv(caminho: Path) -> list[dict]:
    """Saida -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'."""
    if not caminho.exists():
        return []
    linhas = []
    for c in csv.reader(caminho.open(encoding="utf-8"), delimiter="\t"):
        if len(c) < 12:
            continue
        linhas.append({"q": c[0], "s": c[1], "pid": float(c[2]), "L": int(c[3]), "qi": int(c[4]),
                       "qf": int(c[5]), "si": int(c[6]), "sf": int(c[7]), "e": float(c[8]),
                       "bits": float(c[9]), "qlen": int(c[10]), "slen": int(c[11])})
    return linhas


class _Texto:
    """Uma linha de PAF com cara de Path, para reaproveitar s29.le_paf linha a linha."""

    def __init__(self, texto: str) -> None:
        self.texto = texto

    def exists(self) -> bool:
        return True

    def read_text(self, encoding: str = "utf-8") -> str:
        return self.texto


def le_paf_nomeado(caminho: Path) -> list[tuple[str, Any]]:
    """(consulta, Aln do s29 com fs/st/lesoes) por linha de alinhamento; ignora ##STA."""
    if not caminho.exists():
        return []
    saida = []
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        if linha.startswith("##") or not linha.strip():
            continue
        alns = s29.le_paf(_Texto(linha))
        if alns:
            saida.append((linha.split("\t", 1)[0], alns[0]))
    return saida


def melhor(hsps: list[dict], prefixo: str) -> dict | None:
    sel = [h for h in hsps if h["q"].startswith(prefixo)]
    return min(sel, key=lambda h: (h["e"], -h["bits"])) if sel else None


def estrutura_cds(acesso: str, ini: int, fim: int, proteina_id: str) -> list[tuple[int, int]]:
    """Exons codificantes (coordenadas absolutas, na ordem da traducao) da proteina."""
    r = s27.S.get(f"{s27.EU}/efetch.fcgi", timeout=120,
                  params={**s27.CM, "db": "nuccore", "id": acesso, "rettype": "ft", "retmode": "text",
                          "seq_start": ini, "seq_stop": fim})
    r.raise_for_status()
    atual: list[tuple[int, int]] | None = None
    for linha in r.text.splitlines():
        if m := re.match(r"^<?(\d+)\t>?(\d+)\tCDS$", linha):
            atual = [(int(m.group(1)), int(m.group(2)))]
        elif atual is not None and (m := re.match(r"^<?(\d+)\t>?(\d+)$", linha)):
            atual.append((int(m.group(1)), int(m.group(2))))
        elif atual is not None and "protein_id" in linha:
            if proteina_id in linha:
                return [(ini - 1 + a, ini - 1 + b) for a, b in atual]
            atual = None
    return []


def modelo() -> int:
    """Decide se ha proteina do jacare para cada gene e escreve modelo_<gene>_Amis.faa."""
    A = OUT / "aln"  # noqa: N806
    bp = le_tsv(A / "e1_bp_consultas_vs_Amis.tsv")
    janela_genes = [r for r in csv.DictReader((OUT / "s32_janela_genes.tsv").open(encoding="utf-8"),
                                              delimiter="\t") if r["especie"] == "Amis"]
    proteoma = fasta(PROTEOMA.read_text(encoding="utf-8"))
    decisao: dict[str, dict] = {}

    # ODAM: precisa de hit no proteoma OU alinhamento do miniprot com >= 30% da consulta no intervalo
    h = melhor(bp, "ODAM_")
    cob = max((a.cobertura for n, a in le_paf_nomeado(A / "mp_e1_odam_Amis.paf") if n.startswith("ODAM_")),
              default=0.0)
    decisao["ODAM"] = {"melhor_blastp": h, "miniprot_cobertura_max": round(cob, 3),
                       "modelo": None}
    if h and h["e"] <= LIMIAR_E or cob >= 0.3:
        decisao["ODAM"]["situacao"] = "candidato — revisar a mao antes de seguir"
    else:
        decisao["ODAM"]["situacao"] = (f"nao localizado: melhor blastp E={h['e'] if h else 'nenhum'} "
                                       f"(> {LIMIAR_E}), miniprot cobre no maximo {cob:.0%} da consulta")
    print(f"  ODAM: {decisao['ODAM']['situacao']}")

    # ODAPH: melhor hit no proteoma, exigido dentro da janela CDKL2..RCHY1
    h = melhor(bp, "ODAPH_")
    na_janela = {g["proteina"]: g for g in janela_genes if g["proteina"]}
    d: dict = {"melhor_blastp": h, "modelo": None}
    if h and h["e"] <= LIMIAR_E and h["s"] in na_janela:
        g = na_janela[h["s"]]
        seq = proteoma[h["s"]]
        escreve_fasta(OUT / "modelo_odaph_Amis.faa", {"ODAPH_Amis": seq})
        info = json.loads((OUT / "s32_janelas_odaph.json").read_text(encoding="utf-8"))["Amis"]
        exons = estrutura_cds(info["acesso"], int(g["ini"]), int(g["fim"]), h["s"])
        tam = [abs(b - a) + 1 for a, b in exons]
        fases = [sum(tam[:i + 1]) % 3 for i in range(len(tam) - 1)]
        # a traducao que o miniprot tirou do GENOMA tem de estar dentro do modelo RefSeq
        sta, ultima = [], ""
        for ln in (A / "mp_e1_odaph_Amis.paf").read_text(encoding="utf-8").splitlines():
            if ln.startswith("##STA") and ultima.startswith("ODAPH_"):
                sta.append(ln.split("\t")[1])
            elif not ln.startswith("##"):
                ultima = ln
        d.update({"modelo": h["s"], "gene": g["simbolo"], "acesso": info["acesso"],
                  "gene_ini": int(g["ini"]), "gene_fim": int(g["fim"]), "aa": len(seq),
                  "exons_cds": exons, "fases_intron": fases,
                  "motivo_SSEE": [m.start() + 1 for m in re.finditer("SSEE", seq)],
                  "cisteinas": [i + 1 for i, c in enumerate(seq) if c == "C"],
                  "miniprot_no_modelo": [s.rstrip("*") in seq for s in sta],
                  "situacao": "localizado"})
        # o gene inteiro (exons + introns) do jacare, para blastn contra os vaos:
        # reliquia nucleotidica sobrevive onde a traducao ja nao alinha
        corpo = (OUT / "janela_odaph_Amis.fna").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
        a, b = min(e for ex in exons for e in ex), max(e for ex in exons for e in ex)
        escreve_fasta(OUT / "gene_odaph_Amis.fna",
                      {f"ODAPH_Amis_gene {info['acesso']}:{a}-{b}": corpo[a - info["ini"]:b - info["ini"] + 1]})
        print(f"  ODAPH: {h['s']} ({g['simbolo']}) {info['acesso']}:{g['ini']}-{g['fim']}, {len(seq)} aa, "
              f"{len(exons)} exons codificantes, fases dos introns {fases}, SSEE em {d['motivo_SSEE']}")
    else:
        d["situacao"] = "sem hit significativo na janela"
        print(f"  ODAPH: {d['situacao']}")
    decisao["ODAPH"] = d
    (OUT / "s32_modelo.json").write_text(json.dumps(decisao, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


def _melhor_paf(caminho: Path, prefixo: str) -> Any:
    alns = [a for n, a in le_paf_nomeado(caminho) if n.startswith(prefixo)]
    return max(alns, key=lambda a: a.escore) if alns else None


def _genomico(base: int, a: int, b: int) -> str:
    return f"{base + min(a, b) - 1}-{base + max(a, b) - 1}"


def tabela() -> int:
    """Escada de evidencia por gene e especie -> s32_odam_odaph.tsv (+ detalhe)."""
    A = OUT / "aln"  # noqa: N806
    janelas = json.loads((OUT / "s32_janelas_odaph.json").read_text(encoding="utf-8"))
    vao = json.loads((OUT / "s32_vaos.json").read_text(encoding="utf-8"))
    decisao = json.loads((OUT / "s32_modelo.json").read_text(encoding="utf-8"))
    gene_bp = 0
    if decisao["ODAPH"].get("exons_cds"):
        pos = [e for ex in decisao["ODAPH"]["exons_cds"] for e in ex]
        gene_bp = max(pos) - min(pos) + 1
    linhas, detalhe = [], []

    # ---- controle positivo das consultas de ODAPH: achou nas especies com dente?
    # (a consulta nunca conta contra a propria especie: Amis x modelo do jacare e
    # Pmur x ORF do Podarcis ficam de fora)
    achou_com_dente = {}
    for _tx, rot, _nome, dente, _grupo in ESPECIES:
        if dente != "COM":
            continue
        es = []
        for q, origem in (("e2_vao", "Amis"), ("e2_vaoPmur", "Pmur")):
            if rot != origem:
                es += [h["e"] for h in le_tsv(A / f"tbn_noseg_{q}_{rot}.tsv")]
        if es:
            achou_com_dente[rot] = min(es)
    ctrl_dente = ", ".join(f"{r} E={e:.0e}" for r, e in achou_com_dente.items() if e <= LIMIAR_E)

    # ---- controle de sensibilidade do ODAM: consultas de fora contra a regiao do ODAM de Xenopus
    tx = le_tsv(A / "tbn_noseg_e1_odam_ctrlXtro.tsv")
    mx = [(n, a) for n, a in le_paf_nomeado(A / "mp_e1_odam_ctrlXtro.paf") if n != "ODAM_Xtro"]
    sens = "; ".join(
        f"{q}: tblastn E={h['e']:.2g}" if (h := melhor(tx, q)) else f"{q}: tblastn nada"
        for q in ("ODAM_Hsap", "ODAM_Rbiv"))
    sens += "; miniprot " + (", ".join(f"{n} {a.qini + 1}-{a.qfim}/{a.qlen} aa" for n, a in mx) or "nada")
    print(f"  controle ODAM no Xenopus: {sens}")

    for _tx, rot, nome, dente, grupo in ESPECIES:
        base = {"especie": rot, "nome": nome, "grupo": grupo, "dente": dente}
        # ================= ODAPH
        v = vao.get(rot, {})
        if v.get("situacao") != "ok":
            linhas.append({**base, "gene": "ODAPH", "estado": "nao testado",
                           "nivel": "1 (sem anotacao)", "coordenada": "-",
                           "controle": "sem CDKL2/THAP6/RCHY1/G3BP2 anotados: sem janela sintenica",
                           "nota": janelas.get(rot, {}).get("situacao", "")})
        else:
            # ancoras na janela, mesmo comando do modelo (e2_q_odaph.faa)
            jan = le_tsv(A / f"tbn_noseg_e2_odaph_{rot}.tsv")
            anc = {}
            for g in ANCORAS_ODAPH:
                h = melhor(jan, f"{g}_Amis")
                if h and h["e"] <= LIMIAR_E:
                    anc[g] = h["e"]
            # alvo, so no vao: duas consultas (jacare e Podarcis), tblastn + miniprot + blastn
            res = {}
            for q, pref in (("e2_vao", "ODAPH_Amis"), ("e2_vaoPmur", "ODAPH_Pmur")):
                h = melhor(le_tsv(A / f"tbn_noseg_{q}_{rot}.tsv"), pref)
                a = _melhor_paf(A / f"mp_{q}_{rot}.paf", pref)
                res[pref] = (h, a)
                detalhe.append({**base, "consulta": pref, "vao_bp": v["bp"],
                                "tblastn_E": h["e"] if h else "", "tblastn_bits": h["bits"] if h else "",
                                "tblastn_q": f"{h['qi']}-{h['qf']}" if h else "",
                                "tblastn_coord": f"{v['acesso']}:{_genomico(v['ini'], h['si'], h['sf'])}" if h else "",
                                "miniprot_cob": round(a.cobertura, 3) if a else "",
                                "miniprot_ident": round(a.identidade, 3) if a else "",
                                "miniprot_fs": a.fs if a else "", "miniprot_stop": a.st if a else "",
                                "miniprot_coord": f"{v['acesso']}:{_genomico(v['ini'], a.sini + 1, a.sfim)}" if a else "",
                                "lesoes": ";".join(f"{t}@{p}" for t, p in a.lesoes) if a else ""})
            bn = [h for h in le_tsv(A / f"bn_vao_{rot}.tsv") if h["e"] <= LIMIAR_E]
            hs = [(p, h) for p, (h, _) in res.items() if h and h["e"] <= LIMIAR_E]
            ctrl_anc = ", ".join(f"{g} E={e:.0e}" for g, e in anc.items()) or "nenhuma ancora a E<=1e-3"
            if v["bp"] == 0:
                nota = f"vao de 0 bp: {v['parede']} sobrepoe o CDKL2 ({v['acesso']})"
            else:
                nota = f"vao {v['acesso']}:{v['ini']}-{v['fim']} = {v['bp']} bp (N={v['N']}), parede {v['parede']}"
            # janela secundaria em volta do THAP6 deslocado (frango)
            sec = le_tsv(A / f"tbn_noseg_e2_odaph_{rot}_THAP6.tsv")
            if sec:
                hm, ht = melhor(sec, "ODAPH_Amis"), melhor(sec, "THAP6_Amis")
                j2 = janelas.get(f"{rot}_THAP6", {})
                nota += (f"; janela do THAP6 deslocado {j2.get('acesso')}:{j2.get('ini')}-{j2.get('fim')}: "
                         f"THAP6 E={ht['e']:.0e}, ODAPH melhor E={hm['e']:.2g}" if hm and ht else "")
            if v.get("LOC_no_vao"):
                nota += f"; LOC no vao: {','.join(v['LOC_no_vao'])}"
            if hs:
                p, h = min(hs, key=lambda x: x[1]["e"])
                a = res[p][1]
                if a is not None and a.cobertura >= 0.5:
                    coord = f"{v['acesso']}:{_genomico(v['ini'], a.sini + 1, a.sfim)}"
                    les = a.fs + a.st
                    estado = (f"relíquia com lesão ({a.fs} fs, {a.st} stop)" if les else
                              f"ORF sem lesão no trecho alinhado ({p} {a.qini + 1}-{a.qfim}/{a.qlen}, "
                              f"ident {a.identidade:.2f})")
                    nivel = "3" if les else "presente"
                else:
                    coord = f"{v['acesso']}:{_genomico(v['ini'], h['si'], h['sf'])}"
                    estado = (f"presente: tblastn {p} E={h['e']:.1e} (q{h['qi']}-{h['qf']}); "
                              "integridade do ORF nao medida (miniprot < 50% da consulta)")
                    nivel = "presente"
                if bn:
                    nota += f"; blastn do gene do jacare: {len(bn)} HSP"
                linhas.append({**base, "gene": "ODAPH", "estado": estado, "nivel": nivel,
                               "coordenada": coord, "controle": ctrl_anc, "nota": nota})
            else:
                melhor_e = min((h["e"] for h, _ in res.values() if h), default=None)
                sem_espaco = gene_bp and v["bp"] < gene_bp
                estado = (f"nao achado no vao (melhor E={melhor_e:.2g})" if melhor_e is not None
                          else "nao achado no vao (nenhum HSP)")
                if sem_espaco:
                    estado += f"; vao de {v['bp']} bp < {gene_bp} bp do gene no jacare"
                if bn:
                    estado += f"; blastn: {len(bn)} HSP"
                nivel = "2" if anc and ctrl_dente else "1"
                linhas.append({**base, "gene": "ODAPH", "estado": estado, "nivel": nivel,
                               "coordenada": (f"{v['acesso']}:{v['ini']}-{v['fim']} (vao)" if v["bp"]
                                              else f"{v['acesso']}: CDKL2 sem vao"),
                               "controle": f"ancoras: {ctrl_anc}; mesmas consultas acham ODAPH em {ctrl_dente}",
                               "nota": nota})
        # ================= ODAM
        tb = le_tsv(A / f"tbn_noseg_e1_odam_{rot}.tsv")
        h = melhor(tb, "ODAM_")
        cob = max((a.cobertura for n, a in le_paf_nomeado(A / f"mp_e1_odam_{rot}.paf") if n.startswith("ODAM_")),
                  default=0.0)
        ctrl = [g for g in CONTROLES_ODAM if (x := melhor(tb, g)) and x["e"] <= LIMIAR_E]
        linhas.append({**base, "gene": "ODAM",
                       "estado": f"nao achado no intervalo NRG1..LPL (melhor tblastn E={h['e']:.2g}; "
                                 f"miniprot cobre {cob:.0%} da consulta)" if h else "nao achado (nenhum HSP)",
                       "nivel": "1 (negativo sem controle positivo em Sauropsida)",
                       "coordenada": "-",
                       "controle": f"ancoras E<=1e-3: {','.join(ctrl) or 'nenhuma'}; ODAM nao localizado em "
                                   f"nenhuma especie COM dente; sensibilidade no ODAM do Xenopus: {sens}",
                       "nota": ""})

    campos = ["gene", "especie", "nome", "grupo", "dente", "estado", "nivel", "coordenada", "controle", "nota"]
    for nome_arq, dados, cab in (("s32_odam_odaph.tsv", linhas, campos),
                                 ("s32_odaph_detalhe.tsv", detalhe, list(detalhe[0]) if detalhe else [])):
        with (OUT / nome_arq).open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cab, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(sorted(dados, key=lambda d: (d.get("gene", ""), d["especie"])) if nome_arq.startswith(
                "s32_odam") else dados)
    for d in linhas:
        if d["gene"] == "ODAPH":
            print(f"  ODAPH {d['especie']:5s} {d['dente']:3s} nivel {d['nivel']:9s} {d['estado'][:90]}")
    print(f"\n  tabelas: {OUT / 's32_odam_odaph.tsv'}  |  {OUT / 's32_odaph_detalhe.tsv'}")
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    acoes = {"baixar": baixar, "vaos": vaos, "modelo": modelo, "tabela": tabela}
    if cmd not in acoes:
        raise SystemExit(f"uso: s32_odam_odaph.py {{{'|'.join(acoes)}}}")
    raise SystemExit(acoes[cmd]())
