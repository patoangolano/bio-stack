"""Copia uma resposta JSON pública de MCP para o diretório desta pesquisa."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    """Parse and save tool output without changing the original managed file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    data = json.loads(args.source.read_text(encoding="utf-8"))
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    args.destination.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"saved": str(args.destination), "keys": list(data)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
