# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
SessionStart hook - injeta o estado real do ambiente no contexto.

Sem isto, o modelo adivinha: propoe rodar AlphaFold local numa GPU de 8 GB,
sugere conda que nao existe, ou tenta um NIM sem chave configurada. Com isto,
ele comeca a sessao sabendo exatamente o que esta disponivel e o que nao esta.

Saida vai para stdout e entra no contexto como informacao de sistema.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# O console do Windows usa cp1252 por padrao. Forcar UTF-8 aqui e manter a
# saida em ASCII puro (ver abaixo) garante que o texto chegue intacto ao
# modelo, independente de como o consumidor decodifica o stdout.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, OSError):
        pass


def probe_gpu() -> str:
    if not shutil.which("nvidia-smi"):
        return "GPU: nvidia-smi indisponivel."
    try:
        p = subprocess.run(
            ["nvidia-smi",
             "--query-gpu=name,memory.total,memory.used,compute_cap,driver_version",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=15, encoding="utf-8", errors="replace",
        )
        if p.returncode != 0 or not p.stdout.strip():
            return "GPU: nvidia-smi falhou."

        name, total, used, cap, driver = [f.strip() for f in p.stdout.strip().splitlines()[0].split(",")]
        total_gb = int(total) / 1024
        free_gb = (int(total) - int(used)) / 1024

        lines = [
            f"GPU: {name} | {total_gb:.1f} GB total, {free_gb:.1f} GB livres "
            f"| compute {cap} | driver {driver}"
        ]
        if total_gb < 16:
            lines.append(
                "  Orcamento de VRAM: NAO tente rodar ESMFold, AlphaFold2, OpenFold3, "
                "Boltz-2 ou Evo 2 localmente - nao cabem. Use as ferramentas mcp__nvidia-bio__* "
                "(endpoints hospedados). Cabem local: ESM-2 650M, DNABERT-2, scanpy, RDKit, nvMolKit."
            )
        return "\n".join(lines)
    except (subprocess.TimeoutExpired, OSError, ValueError, IndexError):
        return "GPU: leitura falhou."


def probe_tools() -> str:
    wanted = {
        "git": "plugins do Claude Code",
        "uv": "servidores MCP Python",
        "node": "npx skills add (BioNeMo Toolkit)",
        "python": "scripts de analise",
        "nextflow": "pipelines nf-core",
        "docker": "containers NIM auto-hospedados",
        "samtools": "manipulacao de BAM",
        "R": "Bioconductor",
    }
    present: list[str] = []
    missing: list[str] = []
    for tool, why in wanted.items():
        found = shutil.which(tool)
        # O stub da Microsoft Store aparece no PATH mas so abre a loja -
        # reportar "python presente" por causa dele e pior que reportar ausente.
        if tool.startswith("python") and found and "WindowsApps" in found:
            found = None
        if found:
            present.append(tool)
        else:
            missing.append(f"{tool} ({why})")

    out = [f"Ferramentas presentes: {', '.join(present) or 'nenhuma'}"]
    if missing:
        out.append(f"Ausentes: {', '.join(missing)}")
    return "\n".join(out)


def probe_keys() -> str:
    keys = {
        "NVIDIA_API_KEY": "NIM de biologia",
        "NCBI_API_KEY": "rate limit do E-utilities",
        "S2_API_KEY": "Semantic Scholar",
    }
    lines = []
    for k, why in keys.items():
        lines.append(f"  {k}: {'configurada' if os.environ.get(k) else f'AUSENTE ({why})'}")
    return "Chaves:\n" + "\n".join(lines)


def probe_runs() -> str:
    log = Path(os.environ.get("BIO_STACK_PROVENANCE", "runs/provenance.jsonl"))
    if not log.is_file():
        return "Historico: nenhuma analise registrada ainda."
    try:
        lines = log.read_text(encoding="utf-8", errors="replace").strip().splitlines()
        if not lines:
            return "Historico: vazio."
        last = json.loads(lines[-1])
        return (
            f"Historico: {len(lines)} operacoes registradas em {log}. "
            f"Ultima: {last.get('tool')} em {last.get('ts')}."
        )
    except (OSError, json.JSONDecodeError, IndexError):
        return "Historico: log ilegivel."


def main() -> int:
    print("=== bio-stack: estado do ambiente ===")
    print(probe_gpu())
    print(probe_tools())
    print(probe_keys())
    print(probe_runs())
    print(
        "Escopo deste projeto: analise computacional academica - genomica, "
        "transcriptomica, estrutura de proteinas, engenharia de pipelines. "
        "Toda operacao fica registrada em runs/provenance.jsonl para a secao de Metodos."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
