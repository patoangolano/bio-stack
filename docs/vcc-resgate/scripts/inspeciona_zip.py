"""
Inventario de um zip baixado: membros, caminhos suspeitos e CRC.

O zip e tratado como dado nao confiavel: so e lido, nada dele e executado, e
este script fica fora do diretorio do arquivo e roda com python -I.
"""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path, PurePosixPath


def suspeito(nome: str) -> bool:
    """Absoluto, com letra de unidade ou com '..' escaparia do destino (zip slip)."""
    normal = nome.replace("\\", "/")
    partes = PurePosixPath(normal).parts
    return normal.startswith("/") or ":" in normal or ".." in partes


def main() -> None:
    caminho = Path(sys.argv[1])
    print(f"zip: {caminho.name} | {caminho.stat().st_size / 2**20:.1f} MiB")
    with zipfile.ZipFile(caminho) as zf:
        membros = zf.infolist()
        for m in membros:
            print(f"  {m.file_size / 2**20:10.1f} MiB  {m.filename}")
        ruins = [m.filename for m in membros if suspeito(m.filename)]
        print(f"caminhos suspeitos: {ruins or 'nenhum'}")
        corrompido = zf.testzip()
        print(f"CRC de {len(membros)} membros: {'OK' if corrompido is None else 'FALHOU em ' + corrompido}")


if __name__ == "__main__":
    main()
