"""
Chama uma ferramenta do vcc-mcp por JSON-RPC, sem depender do cliente MCP.

Uso: python -I chama_ferramenta.py <ferramenta> '<argumentos em JSON>'
"""

from __future__ import annotations

import json
import subprocess
import sys

SERVIDOR = "C:/Users/Manec/bio-stack/servers/vcc_mcp/server.py"


def main() -> None:
    ferramenta = sys.argv[1]
    argumentos = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    proc = subprocess.Popen(
        ["uv", "run", "--script", SERVIDOR],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        encoding="utf-8",
    )
    # Popen tipa os dois como opcionais; com PIPE eles sempre existem.
    assert proc.stdin is not None and proc.stdout is not None
    entrada, saida = proc.stdin, proc.stdout

    def enviar(mensagem: dict) -> None:
        entrada.write(json.dumps(mensagem) + "\n")
        entrada.flush()

    def receber(id_esperado: int) -> dict:
        for linha in saida:
            linha = linha.strip()
            if linha.startswith("{"):
                dados = json.loads(linha)
                if dados.get("id") == id_esperado:
                    return dados
        raise RuntimeError(f"servidor fechou antes da resposta {id_esperado}")

    try:
        enviar({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                           "clientInfo": {"name": "chamada-direta", "version": "1"}}})
        receber(1)
        enviar({"jsonrpc": "2.0", "method": "notifications/initialized"})
        enviar({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                "params": {"name": ferramenta, "arguments": argumentos}})
        resultado = receber(2).get("result", {})
        texto = " ".join(c.get("text", "") for c in resultado.get("content", []))
        if resultado.get("isError"):
            print(f"ERRO: {texto}")
            sys.exit(1)
        print(json.dumps(json.loads(texto), indent=2, ensure_ascii=False))
    finally:
        entrada.close()
        proc.terminate()


if __name__ == "__main__":
    main()
