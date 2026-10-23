# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""
Suite da camada de CLIs do WSL (bootstrap/pixi-global.toml).

    uv run --script tests/check_wsl_ferramentas.py

Duas perguntas, e as duas ja deram resposta errada em silencio nesta stack:

1. O manifesto versionado e o que esta instalado concordam? Cada ambiente
   declarado existe no ~/.pixi/envs, e o snapshot (bootstrap/pixi-global.
   snapshot.tsv) registra a mesma versao que o manifesto fixa — snapshot
   velho faz a secao de Metodos citar uma versao que nao rodou.
2. Cada ferramenta que a documentacao promete responde? Executa o binario
   de verdade e confere a versao impressa, em vez de so `command -v`: binario
   presente que aborta por biblioteca faltando passa no `command -v`.

Tudo numa chamada so ao WSL (abrir o WSL custa ~2,5-3 s com a distro parada,
porque ele espera o systemd subir; com ela de pe, ~0,1-0,2 s), com o
script em arquivo e `bash -l` — as armadilhas 1 e 3 do CLAUDE.md. Sem WSL a
suite e pulada, com aviso, e sai com 0.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

for _s in (sys.stdout, sys.stderr):
    _reconf = getattr(_s, "reconfigure", None)
    if _reconf is not None:
        _reconf(encoding="utf-8", errors="replace")

RAIZ = Path(__file__).resolve().parents[1]
MANIFESTO = RAIZ / "bootstrap" / "pixi-global.toml"
SNAPSHOT = RAIZ / "bootstrap" / "pixi-global.snapshot.tsv"
DISTRO = "Ubuntu-24.04"

#: binario -> (comando de versao, trecho que tem de aparecer na saida).
#: O trecho e a versao fixada no manifesto; quando o binario nao imprime a
#: versao (gotree), basta ele rodar um comando real.
SONDAS: dict[str, tuple[str, str]] = {
    "iqtree3": ("iqtree3 --version", "3.1.4"),
    "mafft": ("mafft --version", "7.525"),
    "famsa": ("famsa -help", "2.4.1"),
    "muscle": ("muscle -version", "5.3"),
    "clipkit": ("clipkit --version", "2.14.0"),
    "trimal": ("trimal --version", "1.5"),
    "raxml-ng": ("raxml-ng --version", "2.0.3"),
    "FastTree": ("FastTree -expert", "2.2.0"),
    "wastral": ("wastral -h", "1.25"),
    "goalign": ("goalign version", "0.4.1"),
    "gotree": ("printf '((a,b),(c,d));' | gotree stats", "tips"),
    "hyphy": ("hyphy --version", "2.5.103"),
    "codeml": ("echo | codeml /dev/null", "4.10.10"),
    "macse": ("macse -help", "2.07"),
    # O pacote bioconda e 251117, mas o binario se apresenta como "v.250331" —
    # registrado assim no snapshot e no Metodos; aqui basta ele carregar.
    "prank": ("prank -version", "PRANK v."),
    "blastn": ("blastn -version", "2.17.0"),
    "diamond": ("diamond version", "2.2.8"),
    "hmmsearch": ("hmmsearch -h", "3.4"),
    "mmseqs": ("mmseqs version", "18"),
    "miniprot": ("miniprot --version", "0.18"),
    # O -v pinta "1.5" em negrito mesmo com NO_COLOR; sem tirar o ANSI o
    # trecho "3.1.5" nunca casa.
    "orthofinder": (r"orthofinder -v | sed 's/\x1b\[[0-9;]*m//g'", "v3.1.5"),
    "minimap2": ("minimap2 --version", "2.31"),
    "featureCounts": ("featureCounts -v", "2.1.1"),
    "stringtie": ("stringtie --version", "3.0.3"),
    "samtools": ("samtools --version", "1.24"),
    "bcftools": ("bcftools --version", "1.24"),
    "seqkit": ("seqkit version", "2.14.0"),
    "csvtk": ("csvtk version", "0.38.0"),
    "taxonkit": ("taxonkit version", "0.20.0"),
    # A base (taxdump) e baixada pelo bootstrap; sem ela o taxonkit nao responde.
    "taxonkit-base": ("echo 8492 | taxonkit lineage", "Archosauria"),
    "datasets": ("datasets --version", "18.38.0"),
    "efetch": ("efetch -version", "26.0"),
    "gffread": ("gffread --version", "0.12.9"),
    "jq": ("jq --version", "1.8.2"),
    "duckdb": ("duckdb --version", "1.5.6"),
    "parallel": ("parallel --version", "20260922"),
    "rclone": ("rclone version", "v1.75.2"),
    "pigz": ("pigz --version", "2.8"),
    "nf-test": ("nf-test version", "0.9.5"),
    # kent v482 nao imprime versao; a linha de uso prova que o binario carrega.
    "twoBitToFa": ("twoBitToFa", "Convert all or part of .2bit"),
    "liftOver": ("liftOver", "Move annotations"),
    "mafsInRegion": ("mafsInRegion", "Extract MAFS"),
    "agat": ("agat --version", "1.7.0"),
}

falhas: list[str] = []


def ok(cond: bool, msg: str) -> None:
    print(f"  {'OK  ' if cond else 'FALHA'} {msg}")
    if not cond:
        falhas.append(msg)


def _norm(versao: str) -> tuple[str, ...]:
    """Segmentos numericos sem zero a esquerda: '2.07' e '2.7' empatam, como no conda."""
    return tuple(str(int(s)) if s.isdigit() else s for s in versao.split("."))


def para_wsl(p: Path) -> str:
    s = str(p.resolve()).replace("\\", "/")
    return f"/mnt/{s[0].lower()}{s[2:]}" if len(s) > 1 and s[1] == ":" else s


def main() -> int:
    print("== camada de CLIs do WSL (pixi global) ==")
    if not shutil.which("wsl"):
        print("  PULADO: wsl.exe fora do PATH")
        return 0

    manifesto = tomllib.loads(MANIFESTO.read_text(encoding="utf-8"))
    envs = manifesto.get("envs", {})
    fixadas = {
        (env, pacote): spec.removeprefix("==")
        for env, cfg in envs.items()
        for pacote, spec in cfg.get("dependencies", {}).items()
    }

    # Um script so, executado de uma vez. Cada sonda imprime um marcador e a
    # saida; o Python separa depois. Sem `$var` na linha de comando do wsl.exe.
    # Roda num diretorio temporario: o codeml grava rst/rst1/rub e o nf-test
    # grava .nf-test.log no diretorio corrente — que, herdado do pytest, seria
    # a raiz do repositorio (aconteceu na primeira versao desta suite).
    linhas = [
        "export NO_COLOR=1",
        'cd "$(mktemp -d)"',
        "echo '@@ENVS'",
        "ls -1 ~/.pixi/envs 2>/dev/null",
    ]
    # Cada sonda vira funcao exportada e roda sob `timeout 60`, com stdin
    # fechado: binario que espera entrada travaria a suite inteira, e o que
    # trava tem de aparecer como FALHA dele, nao como suite pulada.
    for i, (nome, (cmd, _)) in enumerate(SONDAS.items()):
        linhas += [
            f"sonda_{i}() {{ {cmd} ; }}",
            f"export -f sonda_{i}",
            f"echo '@@{nome}'",
            f"timeout 60 bash -c sonda_{i} </dev/null 2>&1 | head -n 40",
            "[ \"${PIPESTATUS[0]}\" = 124 ] && echo 'ESTOUROU 60 s'",
        ]
    with tempfile.NamedTemporaryFile(
        "w", suffix=".sh", delete=False, encoding="utf-8", newline="\n"
    ) as fh:
        fh.write("\n".join(linhas) + "\n")
        script = Path(fh.name)
    try:
        proc = subprocess.run(
            ["wsl", "-d", DISTRO, "--", "bash", "-l", para_wsl(script)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=900,
            check=False,
            env={**os.environ, "MSYS_NO_PATHCONV": "1"},
        )
    except subprocess.TimeoutExpired:
        # O WSL respondeu e algo travou alem dos 60 s por sonda: isso e defeito
        # da camada, nao ausencia dela — pular aqui esconderia o problema.
        ok(False, "as sondas terminam em 900 s (alguma travou fora do timeout de 60 s)")
        return 1
    except OSError as exc:
        print(f"  PULADO: WSL nao respondeu ({exc})")
        return 0
    finally:
        script.unlink(missing_ok=True)

    blocos: dict[str, str] = {}
    atual = None
    for linha in proc.stdout.splitlines():
        if linha.startswith("@@"):
            atual = linha[2:].strip()
            blocos[atual] = ""
        elif atual is not None:
            blocos[atual] += linha + "\n"

    print("\n-- manifesto x instalado --")
    instalados = set(blocos.get("ENVS", "").split())
    faltando = sorted(set(envs) - instalados)
    ok(
        not faltando,
        f"{len(envs)} ambientes do manifesto instalados"
        + (f" (faltam: {faltando})" if faltando else ""),
    )

    with SNAPSHOT.open(encoding="utf-8", newline="") as fh:
        snap = {
            (r["ambiente"], r["pacote"]): r["versao"] for r in csv.DictReader(fh, delimiter="\t")
        }
    # Comparacao no sentido do conda, nao de texto: o `pixi global` regrava o
    # manifesto normalizando "2.07" para "2.7" (medido com o macse), e as duas
    # grafias sao a mesma versao.
    divergentes = [
        f"{e}/{p}: manifesto {v}, snapshot {snap.get((e, p))}"
        for (e, p), v in fixadas.items()
        if _norm(snap.get((e, p)) or "") != _norm(v)
    ]
    ok(
        not divergentes,
        "snapshot concorda com o manifesto" + (f" — {divergentes}" if divergentes else ""),
    )

    print("\n-- cada binario roda (versao fixada conferida quando ele a imprime) --")
    for nome, (_, esperado) in SONDAS.items():
        saida = blocos.get(nome, "")
        ok(
            esperado in saida,
            f"{nome}: espera '{esperado}'"
            + ("" if esperado in saida else f" — veio: {saida.strip()[:120]!r}"),
        )

    print()
    if falhas:
        print(
            f"FALHOU: {len(falhas)} verificacao(oes). Reaplique: bash bootstrap/05_wsl_ferramentas.sh (no WSL)"
        )
        return 1
    print("Camada de CLIs do WSL OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
