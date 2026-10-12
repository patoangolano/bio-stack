#!/usr/bin/env python
"""Recolhe o pacote de um job do Colab de Downloads para runs/colab/, conferindo sha256.

O MCP do Colab roda numa sessao de navegador na propria maquina, entao
`files.download()` deposita o arquivo em Downloads. Este script fecha a ponte.

Por que UM zip e nao varios arquivos soltos: medido em 2026-10-10, cinco
`files.download()` em sequencia resultaram em UM arquivo entregue. O navegador
bloqueia downloads multiplos automaticos de uma mesma pagina, e o Colab nao
reclama — ele reporta "5 arquivos enviados" e quatro nunca saem. Falha
silenciosa. Com um zip unico isso nao acontece, e o manifesto viaja DENTRO do
pacote, entao nao ha como o manifesto chegar sem os arquivos que descreve.

A conferencia de sha256 nao e cerimonia: download de navegador pode truncar em
silencio, e um PDB truncado ainda parece um PDB.

Uso:
    python scripts/biblio/recolher_colab.py mmp20_allmi
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import zipfile
from datetime import UTC, datetime
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
DOWNLOADS = Path.home() / "Downloads"
PROV = Path(os.environ.get("BIO_STACK_PROVENANCE", RAIZ / "runs" / "provenance.jsonl"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("job", help="nome do job; espera Downloads/<job>.zip e cria runs/colab/<job>/")
    args = ap.parse_args()

    pacote = DOWNLOADS / f"{args.job}.zip"
    if not pacote.is_file():
        print(f"pacote nao encontrado: {pacote}", file=sys.stderr)
        return 1

    destino = RAIZ / "runs" / "colab" / args.job
    destino.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(pacote) as z:
        ruins = z.testzip()
        if ruins is not None:
            print(f"zip corrompido na entrada {ruins} — download truncado", file=sys.stderr)
            return 1
        manifesto = json.loads(z.read("manifesto_transporte.json"))
        z.extractall(destino)

    ok, recusados = [], []
    for item in manifesto:
        alvo = destino / item["nome"]
        if not alvo.is_file():
            recusados.append({"nome": item["nome"], "motivo": "ausente no zip"})
            print(f"AUSENTE   {item['nome']}")
            continue
        obtido = hashlib.sha256(alvo.read_bytes()).hexdigest()
        if obtido != item["sha256"]:
            recusados.append({"nome": item["nome"], "motivo": "sha256 diferente"})
            alvo.unlink()  # nao deixa arquivo suspeito entrar na analise
            print(f"RECUSADO  {item['nome']}  (sha256 nao bate; arquivo removido)")
            continue
        ok.append({"nome": item["nome"], "bytes": item["bytes"], "sha256": obtido})
        print(f"OK        {item['bytes']:>9,}  {item['nome']}")

    pacote.unlink()

    try:
        PROV.parent.mkdir(parents=True, exist_ok=True)
        with PROV.open("a", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    {
                        "ts": datetime.now(UTC).isoformat(),
                        "server": "colab-recolher",
                        "tool": "recolher_colab",
                        "mode": "colab",
                        "params": {"job": args.job},
                        "outputs": [str(destino / m["nome"]) for m in ok],
                        "resumo": {"ok": len(ok), "recusados": len(recusados), "arquivos": ok},
                    },
                    ensure_ascii=False,
                )
                + "\n"
            )
    except OSError as exc:
        print(f"AVISO: procedencia nao gravada: {exc}", file=sys.stderr)

    print(f"\ndestino: {destino}")
    print(f"  ok={len(ok)} recusados={len(recusados)}")
    return 1 if recusados else 0


if __name__ == "__main__":
    raise SystemExit(main())
