#!/usr/bin/env python
"""Baixa os PDFs de acesso aberto da bibliografia e grava o manifesto.

Separacao deliberada, seguindo o padrao que o .gitignore do projeto ja usa para
data/raw: o PDF e dado bruto, fica FORA do git (direito autoral e tamanho); o
manifesto e o artefato citavel e E versionado. Quem clonar o repositorio
reconstroi a biblioteca rodando este script com o manifesto.

Incremental: um PDF ja baixado com o mesmo sha256 e pulado, entao rodar de novo
depois de acrescentar cinco artigos custa cinco downloads.

Uso:
    python scripts/biblio/baixar_pdfs.py entrada.json
    python scripts/biblio/baixar_pdfs.py entrada.json --so-centrais
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
DESTINO = RAIZ / "data" / "biblio" / "pdf"
MANIFESTO = RAIZ / "data" / "biblio" / "manifesto.jsonl"
PROV = Path(os.environ.get("BIO_STACK_PROVENANCE", RAIZ / "runs" / "provenance.jsonl"))

# Etiqueta de API: NCBI e Unpaywall pedem identificacao. Sem isso o NCBI
# estrangula a taxa e o Unpaywall recusa.
EMAIL = os.environ.get("BIO_STACK_EMAIL", "limatheus458@gmail.com")
UA = f"bio-stack/1.0 (academico; mailto:{EMAIL})"

# 3 req/s e o teto do NCBI sem chave; com chave sobe para 10. Fico em 3 porque
# os downloads vao para hosts variados e o gargalo nao e o NCBI.
INTERVALO_S = 0.34
PISO_PDF_BYTES = 20_000  # abaixo disso quase sempre e pagina de erro, nao artigo


def _sanitiza(texto: str, limite: int) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", texto).strip("_")[:limite] or "sem_titulo"


def nome_arquivo(art: dict) -> str:
    """Primeiro autor, ano e titulo curto — legivel na pasta e estavel entre execucoes."""
    autor = _sanitiza(str(art.get("autores", "")).split()[0] if art.get("autores") else "anon", 24)
    ano = str(art.get("ano", "0000"))
    titulo = _sanitiza("_".join(str(art.get("titulo", "")).split()[:6]), 56)
    return f"{autor}_{ano}_{titulo}.pdf"


def ja_baixado(caminho: Path) -> str | None:
    """sha256 se o arquivo existe e parece PDF de verdade; None caso contrario."""
    if not caminho.is_file() or caminho.stat().st_size < PISO_PDF_BYTES:
        return None
    dados = caminho.read_bytes()
    if not dados.startswith(b"%PDF"):
        return None
    return hashlib.sha256(dados).hexdigest()


def baixa_um(cliente: httpx.Client, url: str) -> tuple[bytes | None, str]:
    """Devolve (conteudo, motivo). conteudo None quando nao veio PDF utilizavel."""
    try:
        r = cliente.get(url)
    except httpx.HTTPError as exc:
        return None, f"rede: {type(exc).__name__}: {exc}"
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}"
    corpo = r.content
    # Checar o conteudo, nao o Content-Type: editora que serve paywall costuma
    # devolver 200 com text/html, e as vezes devolve HTML com rotulo de PDF.
    if not corpo.startswith(b"%PDF"):
        trecho = corpo[:200].decode("utf-8", errors="replace").replace("\n", " ")
        return None, f"nao e PDF (comeca com {trecho[:80]!r})"
    if len(corpo) < PISO_PDF_BYTES:
        return None, f"PDF pequeno demais ({len(corpo)} bytes), provavelmente pagina de erro"
    return corpo, "ok"


def registra_procedencia(resumo: dict) -> None:
    """Append-only, mesmo formato dos servidores MCP, para o /metodos enxergar."""
    try:
        PROV.parent.mkdir(parents=True, exist_ok=True)
        with PROV.open("a", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    {
                        "ts": datetime.now(UTC).isoformat(),
                        "server": "biblio-baixar-pdfs",
                        "tool": "baixar_pdfs",
                        "mode": "local",
                        "params": {"entrada": resumo.get("entrada")},
                        "outputs": [str(MANIFESTO)],
                        "resumo": {k: v for k, v in resumo.items() if k != "entrada"},
                    },
                    ensure_ascii=False,
                )
                + "\n"
            )
    except OSError as exc:
        print(f"AVISO: procedencia nao gravada em {PROV}: {exc}", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entrada", help="JSON com a lista de artigos (campo 'artigos')")
    ap.add_argument(
        "--so-centrais", action="store_true", help="baixa apenas os marcados central=true"
    )
    args = ap.parse_args()

    dados = json.loads(Path(args.entrada).read_text(encoding="utf-8"))
    artigos = dados.get("artigos", dados if isinstance(dados, list) else [])
    if args.so_centrais:
        artigos = [a for a in artigos if a.get("central")]

    DESTINO.mkdir(parents=True, exist_ok=True)
    linhas: list[dict] = []
    contagem = {"baixado": 0, "ja_tinha": 0, "sem_url": 0, "falhou": 0}

    with httpx.Client(timeout=60, follow_redirects=True, headers={"User-Agent": UA}) as cliente:
        for i, art in enumerate(artigos, 1):
            alvo = DESTINO / nome_arquivo(art)
            url = (art.get("pdf_final") or art.get("pdf_aberto_url") or "").strip()
            registro = {
                "titulo": art.get("titulo"),
                "autores": art.get("autores"),
                "ano": art.get("ano"),
                "revista": art.get("revista"),
                "doi": art.get("doi") or None,
                "pmid": art.get("pmid") or None,
                "pmcid": art.get("pmcid") or None,
                "lacuna": art.get("lacuna"),
                "central": bool(art.get("central")),
                "por_que_importa": art.get("por_que_importa"),
                "pdf_url": url or None,
                "arquivo": None,
                "sha256": None,
                "estado": None,
            }

            existente = ja_baixado(alvo)
            if existente:
                registro |= {"arquivo": alvo.name, "sha256": existente, "estado": "ja_tinha"}
                contagem["ja_tinha"] += 1
            elif not url:
                registro["estado"] = "sem_url_aberta"
                contagem["sem_url"] += 1
            else:
                corpo, motivo = baixa_um(cliente, url)
                time.sleep(INTERVALO_S)
                if corpo is None:
                    registro |= {"estado": "falhou", "motivo": motivo}
                    contagem["falhou"] += 1
                else:
                    alvo.write_bytes(corpo)
                    registro |= {
                        "arquivo": alvo.name,
                        "sha256": hashlib.sha256(corpo).hexdigest(),
                        "bytes": len(corpo),
                        "estado": "baixado",
                    }
                    contagem["baixado"] += 1

            linhas.append(registro)
            print(f"[{i:3d}/{len(artigos)}] {registro['estado']:14s} {alvo.name[:62]}")

    MANIFESTO.parent.mkdir(parents=True, exist_ok=True)
    with MANIFESTO.open("w", encoding="utf-8") as fh:
        for linha in linhas:
            fh.write(json.dumps(linha, ensure_ascii=False) + "\n")

    print(f"\nmanifesto: {MANIFESTO}  ({len(linhas)} registros)")
    print("  " + "  ".join(f"{k}={v}" for k, v in contagem.items()))
    registra_procedencia({"entrada": args.entrada, **contagem, "total": len(linhas)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
