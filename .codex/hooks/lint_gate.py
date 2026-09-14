# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
PostToolUse hook - ponte entre os language servers e o Claude.

O Claude Code nao fala LSP nativamente. Mas os mesmos analisadores que
alimentam o seu editor tambem rodam em modo batch, e o resultado pode voltar
para o modelo pelo canal de feedback do hook. Efeito pratico: quando o Claude
escreve um script de analise com um erro de tipo ou uma variavel nao definida,
ele ve o diagnostico na hora e corrige, em vez de voce descobrir depois de
esperar 40 minutos de pipeline.

Analisadores por extensao:
  .py            ruff check  +  basedpyright   (mesmo motor do Pylance)
  .nf / .config  nextflow lint                  (mesmo parser do Nextflow LSP)
  .smk, Snakefile snakemake --lint
  .R             Rscript -e lintr::lint
  .json/.yaml    validacao de sintaxe embutida

Cada analisador ausente e simplesmente pulado - o hook nunca quebra por falta
de ferramenta.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

# Diagnosticos de linter costumam trazer caracteres nao-ASCII; sem forcar UTF-8
# eles saem corrompidos no console do Windows.
try:
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except (AttributeError, OSError):
    pass

TIMEOUT = 45
MAX_DIAGNOSTICS = 25


# Os linters sao instalados dentro de .venv, que NAO esta no PATH do processo
# do hook - resolver por shutil.which fazia a ponte inteira ficar inerte numa
# instalacao que seguiu a propria documentacao. `uvx` resolve e cacheia o
# pacote sob demanda, independente de venv ativa.
UVX = shutil.which("uv")


def have(tool: str) -> bool:
    """Disponivel diretamente no PATH, ou alcancavel via uvx."""
    return bool(shutil.which(tool) or UVX)


def tool_cmd(tool: str, *args: str) -> list[str]:
    if shutil.which(tool):
        return [tool, *args]
    return [UVX, "tool", "run", "--quiet", tool, *args]


def run(cmd: list[str], cwd: Path | None = None) -> tuple[int, str, str]:
    try:
        p = subprocess.run(
            cmd, capture_output=True, text=True, timeout=TIMEOUT,
            cwd=str(cwd) if cwd else None,
            encoding="utf-8", errors="replace",
        )
        return p.returncode, p.stdout, p.stderr
    except (subprocess.TimeoutExpired, OSError):
        return -1, "", ""


def check_python(path: Path) -> list[str]:
    out: list[str] = []

    if have("ruff"):
        code, stdout, _ = run(tool_cmd("ruff", "check", "--output-format", "concise", str(path)))
        if code == 1 and stdout.strip():
            out += [line for line in stdout.strip().splitlines() if line.strip()]

    # basedpyright e pyright emitem o MESMO formato de --outputjson; o primeiro
    # que existir serve.
    checker = "basedpyright" if shutil.which("basedpyright") else (
        "pyright" if shutil.which("pyright") else ("basedpyright" if UVX else None)
    )
    if checker:
        _, stdout, _ = run(tool_cmd(checker, "--outputjson", str(path)))
        if stdout.strip():
            try:
                data = json.loads(stdout)
                for d in data.get("generalDiagnostics", []):
                    if d.get("severity") != "error":
                        continue
                    line = (d.get("range") or {}).get("start", {}).get("line", 0) + 1
                    msg = (d.get("message") or "").splitlines()
                    out.append(f"{path.name}:{line}: {msg[0] if msg else '?'}")
            except (json.JSONDecodeError, KeyError, IndexError, TypeError):
                pass
    return out


def check_nextflow(path: Path) -> list[str]:
    if not have("nextflow"):
        return []
    code, stdout, stderr = run(["nextflow", "lint", str(path)], cwd=path.parent)
    if code != 0:
        blob = (stderr or stdout).strip()
        return [line for line in blob.splitlines() if line.strip()][:MAX_DIAGNOSTICS]
    return []


def check_snakemake(path: Path) -> list[str]:
    if not have("snakemake"):
        return []
    code, _, stderr = run(tool_cmd("snakemake", "--lint", "-s", str(path)), cwd=path.parent)
    if code != 0 and stderr.strip():
        return [line for line in stderr.strip().splitlines() if line.strip()][:MAX_DIAGNOSTICS]
    return []


def check_r(path: Path) -> list[str]:
    if not have("Rscript"):
        return []
    expr = f'suppressMessages(print(lintr::lint("{path.as_posix()}")))'
    code, stdout, _ = run(["Rscript", "-e", expr])
    if stdout.strip():
        return [line for line in stdout.strip().splitlines() if line.strip()][:MAX_DIAGNOSTICS]
    return []


def check_structured(path: Path) -> list[str]:
    """Sintaxe de JSON/YAML - pega virgula sobrando em schema de nf-core."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix == ".json":
        try:
            json.loads(text)
        except json.JSONDecodeError as exc:
            return [f"{path.name}:{exc.lineno}: JSON invalido - {exc.msg}"]
    elif path.suffix in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore
            yaml.safe_load(text)
        except ImportError:
            return []
        except Exception as exc:  # yaml.YAMLError
            return [f"{path.name}: YAML invalido - {str(exc).splitlines()[0]}"]
    return []


DISPATCH = {
    ".py": check_python,
    ".nf": check_nextflow,
    ".smk": check_snakemake,
    ".R": check_r,
    ".r": check_r,
    ".json": check_structured,
    ".yaml": check_structured,
    ".yml": check_structured,
}


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    if payload.get("tool_name") not in ("Write", "Edit", "NotebookEdit"):
        return 0

    ti = payload.get("tool_input") or {}
    # NotebookEdit usa notebook_path; sem isto o ramo de notebook nunca disparava.
    raw = ti.get("file_path") or ti.get("notebook_path")
    if not raw:
        return 0

    path = Path(raw)
    if not path.is_file():
        return 0

    checker = DISPATCH.get(path.suffix)
    if checker is None and path.name in ("Snakefile", "nextflow.config"):
        checker = check_snakemake if path.name == "Snakefile" else check_nextflow
    if checker is None:
        return 0

    try:
        diagnostics = checker(path)
    except Exception:
        return 0  # um analisador quebrado nunca deve travar a edicao

    if not diagnostics:
        return 0

    shown = diagnostics[:MAX_DIAGNOSTICS]
    extra = len(diagnostics) - len(shown)

    print(
        f"[lint_gate] {len(diagnostics)} diagnostico(s) em {path.name}:\n"
        + "\n".join(f"  {d}" for d in shown)
        + (f"\n  ... e mais {extra}" if extra else "")
        + "\nCorrija antes de rodar a analise.",
        file=sys.stderr,
    )
    # exit 2 = devolve o stderr ao modelo para ele agir; nao bloqueia a edicao,
    # que ja aconteceu (este e um PostToolUse).
    return 2


if __name__ == "__main__":
    sys.exit(main())
