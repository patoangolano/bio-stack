#!/usr/bin/env python
"""Segunda passada nos artigos que o baixar_pdfs.py nao conseguiu.

Diagnostico da primeira passada (95 artigos, 42 baixados):

    25x HTTP 403   academic.oup.com, pmc.ncbi, europepmc
     9x nao e PDF  paywall servindo HTML com 200
     1x HTTP 202   Europe PMC ainda gerando o render
     1x erro de rede
    17x sem URL de acesso aberto

Os 403 do PMC e do Europe PMC NAO sao paywall — sao bloqueio de cliente
automatizado na URL de navegador. A correcao nao e insistir na mesma URL com
outro cabecalho, e sim perguntar a quem sabe onde a copia aberta esta.

Ordem de tentativa, da melhor para a pior:

  1. Unpaywall (api.unpaywall.org) -> best_oa_location e depois TODAS as
     oa_locations. Ele indexa copias em repositorio (PMC, institucional,
     arXiv/bioRxiv), e repositorio nao bloqueia cliente scriptado como a
     editora bloqueia.
  2. Europe PMC REST fullTextXML -> texto completo em XML para quem tem
     PMCID. Medido: nem todo PMCID tem (devolve 500 quando o artigo esta
     so como resumo). Quando tem, e MELHOR que PDF para indexar: ja vem
     parseado, sem erro de extracao de coluna dupla, legenda ou tabela.

O que sobrar e paywall de verdade, e isso fica declarado no manifesto em vez
de virar silencio.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import httpx

RAIZ = Path(__file__).resolve().parents[2]
PDF_DIR = RAIZ / "data" / "biblio" / "pdf"
TXT_DIR = RAIZ / "data" / "biblio" / "fulltext"
MANIFESTO = RAIZ / "data" / "biblio" / "manifesto.jsonl"
PROV = Path(os.environ.get("BIO_STACK_PROVENANCE", RAIZ / "runs" / "provenance.jsonl"))

EMAIL = os.environ.get("BIO_STACK_EMAIL", "limatheus458@gmail.com")
UA = f"bio-stack/1.0 (academico; mailto:{EMAIL})"
PISO_PDF = 20_000
PISO_TXT = 4_000
INTERVALO_S = 0.2


def _sanitiza(texto: str, limite: int) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", texto).strip("_")[:limite] or "sem_titulo"


def nome_base(art: dict) -> str:
    autor = _sanitiza(str(art.get("autores") or "anon").split()[0], 24)
    titulo = _sanitiza("_".join(str(art.get("titulo", "")).split()[:6]), 56)
    return f"{autor}_{art.get('ano', '0000')}_{titulo}"


def urls_unpaywall(cliente: httpx.Client, doi: str) -> list[str]:
    """Todas as URLs de PDF que o Unpaywall conhece, da melhor para a pior."""
    try:
        r = cliente.get(f"https://api.unpaywall.org/v2/{doi}", params={"email": EMAIL})
        if r.status_code != 200:
            return []
        d = r.json()
    except (httpx.HTTPError, json.JSONDecodeError):
        return []

    vistos: list[str] = []
    locais = [d.get("best_oa_location")] + list(d.get("oa_locations") or [])
    for loc in locais:
        if not loc:
            continue
        # url_for_pdf primeiro: url pode apontar para a landing page.
        for chave in ("url_for_pdf", "url"):
            u = loc.get(chave)
            if u and u not in vistos:
                vistos.append(u)
    return vistos


def baixa_pdf(cliente: httpx.Client, url: str) -> tuple[bytes | None, str]:
    try:
        r = cliente.get(url)
    except httpx.HTTPError as exc:
        return None, f"rede: {type(exc).__name__}"
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}"
    if not r.content.startswith(b"%PDF"):
        return None, "nao e PDF"
    if len(r.content) < PISO_PDF:
        return None, f"PDF pequeno ({len(r.content)} bytes)"
    return r.content, "ok"


def texto_ncbi_efetch(cliente: httpx.Client, pmcid: str) -> tuple[str | None, str]:
    """Texto completo pelo E-utilities do NCBI, a rota oficial com chave.

    Por que esta rota existe alem do Europe PMC: as URLs de navegador do PMC
    (pmc.ncbi.nlm.nih.gov/articles/.../pdf) devolvem 403 para cliente
    scriptado, e a pagina /pmc/articles/PMC... devolve HTML de destino, nao o
    PDF. O efetch e a interface sancionada, aceita a NCBI_API_KEY e nao
    bloqueia. Medido em 2026-10-10: PMC4386873 devolveu 222 KB com corpo;
    PMC2464697 devolveu 200 com 6 KB e SEM corpo, porque esta no PMC so como
    resumo. Por isso a checagem de corpo abaixo nao e opcional — sem ela eu
    gravaria 6 KB de metadado e contaria como texto completo.
    """
    params = {"db": "pmc", "id": pmcid.replace("PMC", ""), "retmode": "xml"}
    chave = os.environ.get("NCBI_API_KEY")
    if chave:
        params["api_key"] = chave
    try:
        r = cliente.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
    except httpx.HTTPError as exc:
        return None, f"rede: {type(exc).__name__}"
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}"
    bruto = r.text
    if "<body" not in bruto and "<sec" not in bruto:
        return None, f"sem corpo de texto ({len(bruto)} bytes) — so resumo no PMC"
    limpo = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", bruto)).strip()
    if len(limpo) < PISO_TXT:
        return None, f"texto curto ({len(limpo)} chars)"
    return limpo, "ok"


def texto_europepmc(cliente: httpx.Client, pmcid: str) -> tuple[str | None, str]:
    """Texto completo em XML. Melhor que PDF para indexar quando existe."""
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    try:
        r = cliente.get(url)
    except httpx.HTTPError as exc:
        return None, f"rede: {type(exc).__name__}"
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}"
    bruto = r.text
    # Tira marcacao, mantendo o texto. Suficiente para busca semantica; nao e
    # parser de JATS e nao pretende ser.
    sem_tag = re.sub(r"<[^>]+>", " ", bruto)
    limpo = re.sub(r"\s+", " ", sem_tag).strip()
    if len(limpo) < PISO_TXT:
        return None, f"texto curto ({len(limpo)} chars) — provavelmente so resumo"
    return limpo, "ok"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--so-centrais", action="store_true")
    args = ap.parse_args()

    registros = [
        json.loads(linha)
        for linha in MANIFESTO.read_text(encoding="utf-8").splitlines()
        if linha.strip()
    ]
    pendentes = [r for r in registros if r.get("estado") in ("falhou", "sem_url_aberta", "fechado")]
    if args.so_centrais:
        pendentes = [r for r in pendentes if r.get("central")]

    PDF_DIR.mkdir(parents=True, exist_ok=True)
    TXT_DIR.mkdir(parents=True, exist_ok=True)
    por_doi = {r.get("doi"): r for r in registros if r.get("doi")}
    contagem = {"pdf": 0, "texto": 0, "fechado": 0}

    with httpx.Client(timeout=60, follow_redirects=True, headers={"User-Agent": UA}) as cli:
        for i, art in enumerate(pendentes, 1):
            base = nome_base(art)
            doi = art.get("doi")
            resolvido = False
            tentativas: list[str] = []

            for url in urls_unpaywall(cli, doi) if doi else []:
                corpo, motivo = baixa_pdf(cli, url)
                time.sleep(INTERVALO_S)
                tentativas.append(f"{url[:60]} -> {motivo}")
                if corpo:
                    alvo = PDF_DIR / f"{base}.pdf"
                    alvo.write_bytes(corpo)
                    art |= {
                        "estado": "baixado_2a",
                        "arquivo": alvo.name,
                        "sha256": hashlib.sha256(corpo).hexdigest(),
                        "bytes": len(corpo),
                        "pdf_url": url,
                        "via": "unpaywall",
                    }
                    contagem["pdf"] += 1
                    resolvido = True
                    print(f"[{i:3d}/{len(pendentes)}] PDF     {base[:58]}")
                    break

            # efetch antes do Europe PMC: medido como mais eficaz, e a chave do
            # NCBI ja esta no ambiente.
            for rota, fn in (("ncbi efetch", texto_ncbi_efetch), ("europepmc", texto_europepmc)):
                if resolvido or not art.get("pmcid"):
                    break
                txt, motivo = fn(cli, art["pmcid"])
                time.sleep(INTERVALO_S)
                tentativas.append(f"{rota} -> {motivo}")
                if txt:
                    alvo = TXT_DIR / f"{base}.txt"
                    alvo.write_text(txt, encoding="utf-8")
                    art |= {
                        "estado": "texto_2a",
                        "arquivo": alvo.name,
                        "sha256": hashlib.sha256(txt.encode()).hexdigest(),
                        "bytes": len(txt),
                        "via": rota,
                    }
                    contagem["texto"] += 1
                    resolvido = True
                    print(f"[{i:3d}/{len(pendentes)}] TEXTO   {base[:58]}")

            if not resolvido:
                art |= {"estado": "fechado", "tentativas": tentativas[:6]}
                contagem["fechado"] += 1
                print(f"[{i:3d}/{len(pendentes)}] fechado {base[:58]}")

            if doi:
                por_doi[doi] = art

    with MANIFESTO.open("w", encoding="utf-8") as fh:
        for r in registros:
            fh.write(json.dumps(por_doi.get(r.get("doi"), r), ensure_ascii=False) + "\n")

    print(f"\npdf={contagem['pdf']}  texto={contagem['texto']}  fechado={contagem['fechado']}")
    try:
        PROV.parent.mkdir(parents=True, exist_ok=True)
        with PROV.open("a", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    {
                        "ts": datetime.now(UTC).isoformat(),
                        "server": "biblio-segunda-passada",
                        "tool": "segunda_passada",
                        "mode": "local",
                        "params": {"pendentes": len(pendentes)},
                        "outputs": [str(MANIFESTO)],
                        "resumo": contagem,
                    },
                    ensure_ascii=False,
                )
                + "\n"
            )
    except OSError as exc:
        print(f"AVISO: procedencia nao gravada: {exc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
