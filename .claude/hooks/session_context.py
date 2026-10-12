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


#: Distro Linux onde vivem as ferramentas que nao tem build nativo Windows.
WSL_DISTRO = "Ubuntu-24.04"

#: Ferramentas que so existem dentro do WSL. Procura-las no PATH do Windows
#: sempre devolveria "ausente" — foi o que a versao anterior deste hook fazia,
#: e por isso o relatorio de sessao continuou dizendo "nextflow ausente" depois
#: de o nextflow estar instalado e funcionando.
WSL_TOOLS = {
    "nextflow": "pipelines nf-core",
    "nf-core": "lint e scaffold de pipeline",
    "samtools": "manipulacao de BAM",
    "bcftools": "manipulacao de VCF",
    "bedtools": "operacoes de intervalo",
    "seqkit": "manipulacao de FASTA/FASTQ",
    # Entrou em 2026-10-10: o docker agora vive aqui, nao no Windows. E o
    # `-profile docker` do Nextflow depende dele, entao a ausencia dele
    # importa tanto quanto a do nextflow.
    "docker": "-profile docker do Nextflow",
    # Camada pixi global (bootstrap/05_wsl_ferramentas.sh, 2026-10-11). Um
    # representante por familia: se ele falta, a familia inteira falta, e
    # listar as 30 ferramentas so alongaria o relatorio de toda sessao.
    "pixi": "gerenciador das CLIs (bootstrap/pixi-global.toml)",
    "iqtree3": "filogenia interativa (o pipeline usa container)",
    "mafft": "alinhamento multiplo",
    "hyphy": "testes de selecao (RELAX, aBSREL, BUSTED)",
    "codeml": "PAML: dN/dS por ramo e sitio",
    "blastn": "BLAST+ local (tblastn incluso)",
    "miniprot": "proteina contra genoma, com frameshift",
    "orthofinder": "ortologia entre proteomas (entrada do pipeline de filogenia)",
    "datasets": "NCBI Datasets: genomas e ortologos",
    "nf-test": "testes de pipeline Nextflow",
}


def _wsl_presentes() -> set[str] | None:
    """
    Quais ferramentas respondem dentro da distro. `None` = WSL indisponivel.

    Uma unica chamada resolve todas: abrir o WSL custa ~1s, e fazer isso seis
    vezes tornaria o hook de inicio de sessao perceptivelmente lento.
    """
    if not shutil.which("wsl"):
        return None

    # Sem cifrao no comando, de proposito. Passado por subprocess, o wsl.exe
    # engole sequencias `$var` na traducao da linha de comando: um laco
    # `for t in ...; do echo "$t"; done` executa o numero certo de vezes e
    # imprime vazio em todas. Medido nesta maquina. `command -v` aceita varios
    # nomes de uma vez e imprime o caminho de cada um que existe, o que resolve
    # sem precisar de variavel nenhuma.
    try:
        proc = subprocess.run(
            ["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc",
             "command -v " + " ".join(WSL_TOOLS)],
            capture_output=True, text=True, timeout=40, check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0 and not proc.stdout.strip():
        return None

    # `command -v` devolve caminho completo; o nome da ferramenta e o basename.
    return {
        ln.strip().rsplit("/", 1)[-1]
        for ln in proc.stdout.replace("\x00", "").splitlines()
        if ln.strip()
    }


def probe_tools() -> str:
    wanted = {
        "git": "plugins do Claude Code",
        "uv": "servidores MCP Python",
        "node": "npx skills add (BioNeMo Toolkit)",
        "python": "scripts de analise",
        "R": "Bioconductor",
    }
    # `docker` saiu desta lista em 2026-10-10. A partir da reconstrucao do WSL
    # o daemon e o docker-ce nativo DENTRO da Ubuntu-24.04, por systemd, e nao
    # o Docker Desktop. Logo nao existe `docker` no PATH do Windows e nunca vai
    # existir — reportar isso como ausencia descrevia a maquina errada. A
    # presenca real do docker e checada no bloco do WSL, abaixo.
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

    out = [f"Ferramentas presentes (Windows): {', '.join(present) or 'nenhuma'}"]
    if missing:
        out.append(f"Ausentes no Windows: {', '.join(missing)}")

    no_wsl = _wsl_presentes()
    if no_wsl is None:
        out.append(f"WSL ({WSL_DISTRO}): indisponivel — nenhuma ferramenta Linux acessivel.")
    else:
        faltando = [f"{t} ({w})" for t, w in WSL_TOOLS.items() if t not in no_wsl]
        out.append(
            f"Ferramentas no WSL ({WSL_DISTRO}): "
            f"{', '.join(sorted(no_wsl)) or 'nenhuma'}"
        )
        if faltando:
            out.append(f"  ausentes no WSL: {', '.join(faltando)}")
        out.append(
            f'  Como chamar: wsl -d {WSL_DISTRO} -- bash -lc "samtools view ..." '
            "— o projeto fica em /mnt/c/Users/Manec/bio-stack."
        )
    return "\n".join(out)


def probe_keys() -> str:
    # Apenas chaves com consumidor real. O S2_API_KEY saiu daqui em 2026-10-10:
    # era repassado ao `biomcp`, e o biomcp-python 0.7.3 nao tem uma unica
    # ocorrencia de S2_API_KEY nem de semanticscholar. Reportar AUSENTE para
    # uma variavel que ninguem le produz uma lacuna falsa no inicio de toda
    # sessao, e uma lacuna falsa custa mais que nenhuma informacao. Ver
    # docs/CHAVES.md secao 2.
    keys = {
        "NVIDIA_API_KEY": "NIM de biologia",
        "NCBI_API_KEY": "rate limit do E-utilities",
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
