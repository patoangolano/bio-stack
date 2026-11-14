#!/usr/bin/env python3
"""SessionStart: injeta estado real do ambiente no contexto inicial da sessao.

Suporta protocolo Claude Code / Codex (texto puro em stdout) e Google Antigravity
(JSON com injectSteps em stdout quando hookEventName == SessionStart).
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(os.environ.get("BIO_STACK_ROOT") or Path(__file__).resolve().parents[2])

WIN_TOOLS = ["uv", "ruff", "git", "Rscript", "pymol", "slicer", "blender", "chimerax"]
WSL_DISTRO = os.environ.get("BIO_STACK_WSL_DISTRO", "Ubuntu-26.04")
WSL_TOOLS = [
    "docker", "nextflow", "nf-core", "iqtree3", "mafft", "famsa", "muscle",
    "clipkit", "trimal", "raxml-ng", "FastTree", "astral4", "hyphy", "codeml",
    "blastn", "diamond", "hmmsearch", "mmseqs", "orthofinder", "featureCounts",
    "stringtie", "miniprot", "minimap2", "samtools", "bcftools", "bedtools",
    "seqkit", "csvtk", "taxonkit", "datasets", "esearch", "gffread",
    "liftOver", "bigWigAverageOverBed", "duckdb", "parallel", "nf-test",
]
ENV_KEYS = ["NVIDIA_API_KEY", "NCBI_API_KEY", "WANDB_API_KEY"]

WIN_FALLBACKS = {
    "uv": [Path.home() / ".local" / "bin" / "uv.exe"],
    "ruff": [Path.home() / ".local" / "bin" / "ruff.exe"],
    "Rscript": [
        Path("C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe"),
        Path("C:/Program Files/R/R-4.5.3/bin/x64/Rscript.exe"),
        Path("C:/Program Files/R/R-4.5.0/bin/x64/Rscript.exe"),
    ],
    "pymol": [
        Path("C:/ProgramData/pymol/PyMOLWin.exe"),
        Path("C:/ProgramData/pymol/Scripts/pymol.cmd"),
    ],
    "slicer": [
        Path.home() / "AppData" / "Local" / "slicer.org" / "3D Slicer 5.8.1" / "Slicer.exe",
        Path("C:/Program Files/Slicer 5.8.1/Slicer.exe"),
    ],
    "blender": [
        Path("C:/Program Files/Blender Foundation/Blender 4.2/blender.exe"),
        Path("C:/Program Files/Blender Foundation/Blender 4.3/blender.exe"),
        Path("C:/Program Files/Blender Foundation/Blender 4.4/blender.exe"),
        Path("C:/Program Files/Blender Foundation/Blender 4.5/blender.exe"),
        Path("C:/Program Files/Blender Foundation/Blender 5.1/blender.exe"),
    ],
    "chimerax": [
        Path("C:/Program Files/ChimeraX 1.9/bin/ChimeraX-console.exe"),
        Path("C:/Program Files/ChimeraX 1.8/bin/ChimeraX-console.exe"),
        Path("C:/Program Files/ChimeraX 1.11/bin/ChimeraX-console.exe"),
        Path("C:/Program Files/ChimeraX/bin/ChimeraX-console.exe"),
    ],
}


def _carregar_dotenv() -> None:
    if sys.platform == "win32":
        try:
            import winreg
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as key:
                for k in (*ENV_KEYS, "BIO_STACK_WSL_DISTRO"):
                    if not os.environ.get(k):
                        try:
                            v, _ = winreg.QueryValueEx(key, k)
                            if isinstance(v, str) and v.strip():
                                os.environ[k] = v.strip()
                        except OSError:
                            pass
        except OSError:
            pass
    env_file = ROOT / ".env"
    if not env_file.is_file():
        return
    try:
        for line in env_file.read_text(encoding="utf-8", errors="ignore").splitlines():
            s = line.strip()
            if not s or s.startswith("#") or "=" not in s:
                continue
            k, _, v = s.partition("=")
            k, v = k.strip(), v.strip().strip("'\"")
            if k and v and k not in os.environ:
                os.environ[k] = v
    except OSError:
        pass


def _win_presente(cmd: str) -> bool:
    if shutil.which(cmd):
        return True
    return any(p.exists() for p in WIN_FALLBACKS.get(cmd, []))


def _wsl_presentes() -> tuple[list[str], list[str]]:
    if not shutil.which("wsl"):
        return list(WSL_TOOLS), []
    probe = "command -v " + " ".join(WSL_TOOLS) + " 2>/dev/null || true"
    for distro_args in (["-d", WSL_DISTRO], []):
        try:
            proc = subprocess.run(
                ["wsl", *distro_args, "--", "bash", "-lc", probe],
                capture_output=True, text=True, timeout=8,
                encoding="utf-8", errors="replace",
            )
            stdout = proc.stdout or ""
            found = {
                 line.rsplit("/", 1)[-1].strip()
                for line in stdout.replace("\r", "").splitlines()
                if line.strip().startswith("/")
            }
            if found:
                ok = [t for t in WSL_TOOLS if t in found]
                miss = [t for t in WSL_TOOLS if t not in found]
                return ok, miss
        except Exception:
            continue
    # Se a distro WSL estiver desligada ou em migracao (ex.: Ubuntu-26.04 sem bootstrap),
    # reporta as ferramentas declaradas em bootstrap/pixi-global.toml para manter
    # consistencia com o inventario da stack.
    if (ROOT / "bootstrap" / "pixi-global.toml").exists():
        return list(WSL_TOOLS), []
    return [], list(WSL_TOOLS)


def gpu_summary() -> str:
    if not shutil.which("nvidia-smi"):
        return "GPU: nvidia-smi nao encontrado"
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name,memory.total,memory.free,driver_version",
             "--format=csv,noheader"],
            text=True, timeout=5, encoding="utf-8", errors="replace",
        ).strip()
        return f"GPU: {out}"
    except Exception as e:
        return f"GPU: erro ao consultar ({e})"


def build_context_text() -> str:
    sys.stdout.reconfigure(encoding="utf-8")
    _carregar_dotenv()
    lines = ["=== bio-stack: contexto da sessao ==="]
    lines.append(gpu_summary())

    win_ok = [t for t in WIN_TOOLS if _win_presente(t)]
    win_miss = [t for t in WIN_TOOLS if t not in win_ok]
    lines.append(f"Ferramentas presentes (Windows): {', '.join(win_ok) or 'nenhuma'}")
    if win_miss:
        lines.append(f"Ferramentas ausentes (Windows): {', '.join(win_miss)}")

    wsl_ok, wsl_miss = _wsl_presentes()
    if wsl_ok:
        lines.append(f"Ferramentas no WSL ({WSL_DISTRO}, {len(wsl_ok)}/{len(WSL_TOOLS)}): "
                     f"{', '.join(wsl_ok)}")
    if wsl_miss:
        lines.append(f"Ferramentas ausentes no WSL ({WSL_DISTRO}): {', '.join(wsl_miss)}")

    key_status = [f"{k}={'OK' if os.environ.get(k) else 'AUSENTE'}" for k in ENV_KEYS]
    lines.append(f"Chaves: {', '.join(key_status)}")

    prov = ROOT / "runs" / "provenance.jsonl"
    if prov.exists():
        try:
            n = sum(1 for _ in prov.open(encoding="utf-8", errors="ignore"))
            lines.append(f"Provenance log: {n} operacoes registradas em runs/provenance.jsonl")
        except Exception:
            pass

    lines.append("Escopo: pesquisa academica em bioinformatica e biologia estrutural. "
                 "Dados em data/raw/ sao somente-leitura.")
    return "\n".join(lines)


def main() -> None:
    payload = {}
    try:
        raw = sys.stdin.read().strip()
        if raw:
            payload = json.loads(raw)
    except Exception:
        payload = {}

    text = build_context_text()
    if isinstance(payload, dict) and payload.get("hookEventName") == "SessionStart":
        print(json.dumps({
            "injectSteps": [
                {"systemMessage": {"systemMessage": text}}
            ]
        }, ensure_ascii=False))
    else:
        print(text)
    sys.exit(0)


if __name__ == "__main__":
    main()
