"""Cliente BLAST corrigido: o retorno Tabular vem embrulhado em HTML (<PRE>).
O parser antigo filtrava por contagem de TAB e descartava tudo silenciosamente."""

from __future__ import annotations

import html
import re
import sys
import time
from pathlib import Path

import requests

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
SEQ = Path(r"C:\Users\Manec\bio-stack\runs\atavismo\seqs")
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
S = requests.Session()


def getfa(lab):
    for b in (SEQ / "ref_proteinas.fa").read_text(encoding="utf-8").split(">"):
        if b.startswith(lab):
            return ">" + b
    raise KeyError(lab)


def parse(txt: str):
    """Extrai as linhas tabulares de dentro do <PRE>, se houver."""
    blocks = re.findall(r"<PRE>(.*?)</PRE>", txt, re.S | re.I)
    body = "\n".join(html.unescape(b) for b in blocks) if blocks else txt
    rows = []
    for linha in body.splitlines():
        linha = linha.strip()
        if not linha or linha.startswith("#") or linha.startswith("<"):
            continue
        f = linha.split("\t")
        if len(f) >= 12:
            rows.append(f)
    return rows, body


def run(nome, query, db, program="tblastn", entrez=None, evalue="1e-3", ws=None):
    p = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": db,
        "QUERY": query,
        "EXPECT": evalue,
        "HITLIST_SIZE": "100",
        "FILTER": "F",
        "EMAIL": "limatheus458@gmail.com",
        "TOOL": "bio-stack",
    }
    if entrez:
        p["ENTREZ_QUERY"] = entrez
    if ws:
        p["WORD_SIZE"] = ws
    r = S.post(URL, data=p, timeout=200)
    m = re.search(r"RID = (\S+)", r.text)
    if not m:
        print(f"{nome:26} SEM RID")
        return None
    rid = m.group(1)
    st = "?"
    info = ""
    for _ in range(60):
        time.sleep(20)
        info = S.get(
            URL, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}, timeout=150
        ).text
        mm = re.search(r"Status=(\w+)", info)
        st = mm.group(1) if mm else "?"
        if st in ("READY", "UNKNOWN"):
            break
    hits = re.search(r"ThereAreHits=(\w+)", info)
    t = S.get(
        URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "Tabular",
            "ALIGNMENTS": "100",
            "DESCRIPTIONS": "100",
        },
        timeout=200,
    ).text
    rows, body = parse(t)
    (OUT / f"blastfix_{nome}.tsv").write_text(
        "\n".join("\t".join(r_) for r_ in rows), encoding="utf-8"
    )
    print(
        f"{nome:26} RID={rid} status={st} ThereAreHits={hits.group(1) if hits else '?'} -> {len(rows)} HSPs",
        flush=True,
    )
    for f in rows[:8]:
        print(
            f"      {f[1][:24]:24} pid={f[2]:>6} len={f[3]:>5} sbjct={f[8]}-{f[9]} E={f[10]:>10} bits={f[11]}"
        )
    return rows


DB = "GCF_016699485.2"  # genoma de Gallus gallus GRCg7b
print(f"### banco {DB} (genoma de Gallus gallus GRCg7b), parser corrigido\n")
ctrl = run("CTRL_SULT1E1_Ggal", getfa("Ggal_SULT1E1"), DB)
if not ctrl:
    print("\n!! CONTROLE POSITIVO FALHOU DE NOVO -> nao interpretar os negativos")
    sys.exit(0)
print("\n-> controle positivo passou; rodando as consultas de esmalte\n")
for lab in [
    "Amis_ENAM",
    "Amis_AMBN",
    "Amis_AMTN",
    "Amis_AMELX",
    "Amis_DSPP",
    "Cmac_ENAM",
    "Scam_AMELX",
]:
    run(f"TESTE_{lab}", getfa(lab), DB, evalue="1e-2")
