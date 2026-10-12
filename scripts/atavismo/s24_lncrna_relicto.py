"""Os lncRNAs do chrZ do frango sao os restos do cluster de esmalte do jacare?

Continua o s18_chrZ.py. Aquele script estabeleceu a ancora sintenica
(LPL/PSD3 arcossauria, nao a SULT1E1/IGJ mamifera de Sire 2008). A medida de
2026-10-12 estendeu a ancora para quatro genes e fechou o endereco:

    NRG1 -> [ENAM . AMBN . AMTN] -> LPL -> PSD3     mesma ordem, mesma fita

    jacare  NC_081826.1 (chr3)   NRG1..LPL = 893.413 bp, tres genes codificantes
    frango  NC_052572.1 (chrZ)   NRG1..LPL = 105.209 bp, SEIS lncRNAs, zero
                                 gene codificante

Nenhum dos nove genes do painel esta anotado em Gallus gallus, nem no NCBI nem
no Ensembl. Entao "reliquia" aqui nao pode significar pseudogene anotado.

A PERGUNTA: os lncRNAs que ocupam o endereco tem homologia com o cluster?

    Se SIM  -> a reliquia e um lncRNA transcrito, e a persistencia dela ganha
               mecanismo candidato. Liga em Crombez et al. 2025
               (10.1038/s41467-025-...), que reporta sobreposicao de pseudogene
               com lncRNA e Ka/Ks nao-neutro em parte deles.
    Se NAO  -> os 105 kb sao sequencia nao relacionada que ocupou o endereco
               depois da delecao. Tambem e resultado, e e o controle que a
               afirmacao exige.

DESENHO DO TESTE — sem os dois controles, um hit nao significa nada:

    alvo              ENAM, AMBN, AMTN do jacare  -> intervalo do chrZ do frango
    controle POSITIVO LPL do jacare               -> mesmo intervalo
                      LPL esta anotado nas DUAS especies. Se o alinhamento nao
                      recuperar LPL, o metodo nao tem sensibilidade para 116 Myr
                      e nenhuma conclusao negativa sobre o cluster se sustenta.
    controle NEGATIVO DSPP do jacare              -> mesmo intervalo
                      DSPP e dentina, esta no chr2 do jacare, e nao deveria ter
                      nada a ver com este endereco. DSPP tem repeticao
                      Ser-Ser-Asp de baixa complexidade: se ele "acertar", o
                      filtro de complexidade esta desligado e todo hit e suspeito.

Este script so BAIXA e escreve FASTA. O alinhamento roda no WSL, onde vivem o
blast e o miniprot: scripts/atavismo/s24_alinhar.sh

    uv run python scripts/atavismo/s24_lncrna_relicto.py
    wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s24_alinhar.sh
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import requests

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENS = "https://rest.ensembl.org"
CM = {"tool": "bio-stack", "email": "limatheus458@gmail.com"}

OUT = Path(__file__).resolve().parents[2] / "runs" / "atavismo" / "lncrna"

S = requests.Session()
S.headers["User-Agent"] = "bio-stack/s24"

# --- Coordenadas medidas em 2026-10-12, nao inventadas. ----------------------
# Janela do frango: abre 40 kb antes do fim de NRG1 (55.109.784) e fecha 70 kb
# depois do inicio de PSD3 (55.231.641). Pega o intervalo inteiro, LPL (o
# controle positivo) e folga nas duas pontas.
FRANGO_CHRZ = "NC_052572.1"
FRANGO_INI, FRANGO_FIM = 55_070_000, 55_300_000

# Janela do jacare: abre 45 kb antes de ENAM (157.293.890) e fecha depois de
# LPL (157.512.736). Alvo do blastn dos lncRNAs.
JACARE_CHR3 = "NC_081826.1"
JACARE_INI, JACARE_FIM = 157_250_000, 157_550_000

# UID de proteina no NCBI, resolvido por elink a partir do gene_id.
# Primeira isoforma de cada, que basta para tblastn/miniprot.
PROTEINAS = {
    "ENAM": "1928163124",   # alvo
    "AMBN": "1113763746",   # alvo
    "AMTN": "1011589666",   # alvo
    "LPL": "1113763590",    # CONTROLE POSITIVO — anotado nas duas especies
    "DSPP": "2580375102",   # CONTROLE NEGATIVO — chr2 do jacare, baixa complexidade
}

# lncRNAs do Ensembl no intervalo NRG1..LPL do chrZ (GRCg7b). Os cinco
# primeiros estao DENTRO do intervalo; o ultimo fica alem de LPL, do lado de
# PSD3, e entra como sexto para nao perder nada por recorte.
LNCRNAS = [
    "ENSGALG00010014552",
    "ENSGALG00010014554",
    "ENSGALG00010014556",
    "ENSGALG00010014557",
    "ENSGALG00010014761",
    "ENSGALG00010014771",
]


def fatia_genomica(acesso: str, ini: int, fim: int, rotulo: str) -> str:
    """Baixa uma fatia de cromossomo por efetch. Falha alto: sem a fatia nao ha teste."""
    r = S.get(
        f"{EU}/efetch.fcgi",
        params={**CM, "db": "nuccore", "id": acesso, "rettype": "fasta",
                "retmode": "text", "seq_start": ini, "seq_stop": fim},
        timeout=300,
    )
    r.raise_for_status()
    texto = r.text.strip()
    if not texto.startswith(">"):
        raise RuntimeError(f"{rotulo}: efetch nao devolveu FASTA ({texto[:120]!r})")
    # Reescreve o cabecalho para algo que o blast nao trunque e que diga a origem.
    corpo = texto.split("\n", 1)[1]
    return f">{rotulo} {acesso}:{ini}-{fim}\n{corpo}\n"


def proteina(uid: str, sym: str) -> str:
    r = S.get(
        f"{EU}/efetch.fcgi",
        params={**CM, "db": "protein", "id": uid, "rettype": "fasta", "retmode": "text"},
        timeout=120,
    )
    r.raise_for_status()
    texto = r.text.strip()
    if not texto.startswith(">"):
        raise RuntimeError(f"{sym}: efetch de proteina falhou ({texto[:120]!r})")
    cabecalho, corpo = texto.split("\n", 1)
    original = cabecalho[1:].split(" ", 1)[0]
    return f">{sym}_Amis {original}\n{corpo}\n"


def lncrna(ident: str) -> str | None:
    """Sequencia genomica do lncRNA pelo Ensembl. Devolve None se o id nao existir."""
    r = S.get(f"{ENS}/sequence/id/{ident}",
              params={"type": "genomic"},
              headers={"Content-Type": "text/x-fasta"}, timeout=120)
    if r.status_code != 200:
        print(f"    {ident}: HTTP {r.status_code} — pulado", file=sys.stderr)
        return None
    texto = r.text.strip()
    if not texto.startswith(">"):
        return None
    corpo = texto.split("\n", 1)[1]
    return f">{ident}_Ggal lncRNA chrZ GRCg7b\n{corpo}\n"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"saida: {OUT}")

    # 1. Intervalo do frango — o SUJEITO de tblastn e miniprot.
    print("\n[1/4] intervalo do chrZ do frango")
    alvo_frango = OUT / "frango_chrZ_intervalo.fna"
    fasta = fatia_genomica(FRANGO_CHRZ, FRANGO_INI, FRANGO_FIM, "Ggal_chrZ_NRG1_PSD3")
    alvo_frango.write_text(fasta, encoding="utf-8")
    pares = len(fasta.split("\n", 1)[1].replace("\n", ""))
    print(f"      {alvo_frango.name}  {pares:,} bp")

    # 2. Regiao do jacare — o SUJEITO do blastn dos lncRNAs.
    print("\n[2/4] regiao do cluster no jacare")
    time.sleep(0.4)
    alvo_jacare = OUT / "jacare_chr3_cluster.fna"
    fasta = fatia_genomica(JACARE_CHR3, JACARE_INI, JACARE_FIM, "Amis_chr3_cluster_esmalte")
    alvo_jacare.write_text(fasta, encoding="utf-8")
    pares = len(fasta.split("\n", 1)[1].replace("\n", ""))
    print(f"      {alvo_jacare.name}  {pares:,} bp")

    # 3. Proteinas do jacare — as CONSULTAS, alvo e os dois controles.
    print("\n[3/4] proteinas do jacare (3 alvos + 2 controles)")
    blocos = []
    for sym, uid in PROTEINAS.items():
        time.sleep(0.4)
        bloco = proteina(uid, sym)
        blocos.append(bloco)
        n = len(bloco.split("\n", 1)[1].replace("\n", ""))
        papel = {"LPL": "CONTROLE POSITIVO", "DSPP": "CONTROLE NEGATIVO"}.get(sym, "alvo")
        print(f"      {sym:6s} {n:>5} aa   {papel}")
    (OUT / "jacare_proteinas.faa").write_text("".join(blocos), encoding="utf-8")

    # 4. lncRNAs do frango — as CONSULTAS do blastn reciproco.
    print("\n[4/4] lncRNAs do frango no endereco sintenico")
    blocos = []
    for ident in LNCRNAS:
        time.sleep(0.3)
        bloco = lncrna(ident)
        if not bloco:
            continue
        blocos.append(bloco)
        n = len(bloco.split("\n", 1)[1].replace("\n", ""))
        print(f"      {ident}  {n:>6,} bp")
    if not blocos:
        print("      NENHUM lncRNA baixado — o teste reciproco nao roda", file=sys.stderr)
        return 1
    (OUT / "frango_lncrnas.fna").write_text("".join(blocos), encoding="utf-8")

    print("\npronto. agora o alinhamento, no WSL:")
    print("  wsl -d Ubuntu-24.04 -- bash -l "
          "/mnt/c/Users/Manec/bio-stack/scripts/atavismo/s24_alinhar.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
