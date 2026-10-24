"""
Driver JSON-RPC para o servidor vcc-mcp: handshake, inventario e guardrails.

Nenhuma chamada aqui passa confirmar_final=True. Todo vcc_submit deste teste
precisa ser recusado ANTES do upload, e e isso que se mede.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

SERVIDOR = "C:/Users/Manec/bio-stack/servers/vcc_mcp/server.py"
RASCUNHO = Path(sys.argv[1])
AMOSTRA_FALSA = Path("C:/Users/Manec/bio-stack/runs/vcc/sample/guardrail_teste.vcc")
REGISTRO = Path("C:/Users/Manec/bio-stack/runs/vcc/sample/registro.jsonl")


def sha256(caminho: Path) -> str:
    return hashlib.sha256(caminho.read_bytes()).hexdigest()

h5ad_falso = RASCUNHO / "pred_falsa.h5ad"
vcc_falso = RASCUNHO / "pred_falsa.vcc"
for caminho in (h5ad_falso, vcc_falso, AMOSTRA_FALSA):
    caminho.parent.mkdir(parents=True, exist_ok=True)
    caminho.write_bytes(b"nao e uma predicao")

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
ENTRADA = proc.stdin
SAIDA = proc.stdout


def enviar(mensagem: dict) -> None:
    ENTRADA.write(json.dumps(mensagem) + "\n")
    ENTRADA.flush()


def receber(id_esperado: int) -> dict:
    for linha in SAIDA:
        linha = linha.strip()
        if not linha.startswith("{"):
            continue
        dados = json.loads(linha)
        if dados.get("id") == id_esperado:
            return dados
    raise RuntimeError(f"servidor fechou antes da resposta {id_esperado}")


def chamar(id_: int, ferramenta: str, argumentos: dict) -> tuple[bool, str]:
    enviar({"jsonrpc": "2.0", "id": id_, "method": "tools/call",
            "params": {"name": ferramenta, "arguments": argumentos}})
    resposta = receber(id_)
    resultado = resposta.get("result", {})
    texto = " ".join(c.get("text", "") for c in resultado.get("content", []))
    return bool(resultado.get("isError")), texto


try:
    enviar({"jsonrpc": "2.0", "id": 1, "method": "initialize",
            "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                       "clientInfo": {"name": "teste-guardrail", "version": "1"}}})
    info = receber(1)["result"]["serverInfo"]
    print(f"initialize   -> {info['name']}")
    enviar({"jsonrpc": "2.0", "method": "notifications/initialized"})

    enviar({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
    nomes = sorted(t["name"] for t in receber(2)["result"]["tools"])
    print(f"tools/list   -> {len(nomes)}: {', '.join(nomes)}")

    erro, texto = chamar(3, "vcc_whoami", {})
    who = json.loads(texto) if not erro else {}
    print(f"vcc_whoami   -> erro={erro} fase={who.get('fase')} contextos={who.get('contextos')} "
          f"horas_restantes={who.get('horas_restantes')} pacote={who.get('pacote_de_controles')}")

    # Amostra renomeada e movida para fora de runs/vcc/sample/: so o sha256 no
    # registro pode reconhece-la. O registro nao existia antes deste teste.
    amostra_movida = RASCUNHO / "modelo_final_v3.vcc"
    amostra_movida.write_bytes(b"amostra disfarcada de predicao real")
    REGISTRO.parent.mkdir(parents=True, exist_ok=True)
    REGISTRO.write_text(
        json.dumps({"arquivo": "renomeado", "sha256": sha256(amostra_movida)}) + "\n",
        encoding="utf-8",
    )

    casos = [
        ("h5ad cru", {"arquivo_vcc": str(h5ad_falso), "nome_modelo": "teste"}),
        ("amostra na fase final", {"arquivo_vcc": str(AMOSTRA_FALSA), "nome_modelo": "teste"}),
        ("amostra renomeada", {"arquivo_vcc": str(amostra_movida), "nome_modelo": "teste"}),
        ("final sem confirmar", {"arquivo_vcc": str(vcc_falso), "nome_modelo": "teste"}),
    ]
    for i, (rotulo, argumentos) in enumerate(casos, start=10):
        erro, texto = chamar(i, "vcc_submit", argumentos)
        estado = "RECUSADO" if erro else "!! PASSOU"
        print(f"submit [{rotulo:22s}] -> {estado}: {texto[:150]}")
finally:
    ENTRADA.close()
    proc.terminate()
    for caminho in (h5ad_falso, vcc_falso, AMOSTRA_FALSA, RASCUNHO / "modelo_final_v3.vcc", REGISTRO):
        caminho.unlink(missing_ok=True)
