# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""
Suite do pipeline de filogenia (pipelines/filogenia) — roda o perfil de teste
de verdade, no WSL com Docker, e confere o resultado contra a arvore conhecida.

    uv run --script tests/check_filogenia.py

Os locos de teste foram simulados com o AliSim sobre
pipelines/filogenia/assets/test/arvore_verdadeira.nwk e depois desalinhados.
O pipeline tem de realinhar, limpar, inferir e chegar na mesma topologia: a
distancia de Robinson-Foulds entre a arvore concatenada e a verdadeira tem de
ser zero. Um pipeline que roda ate o fim e devolve a arvore errada passaria em
qualquer teste de "terminou com exit 0" — este nao deixa.

Na primeira execucao o Docker baixa ~0,3 GB de imagens, que ocupam ~1 GB no
disco; depois roda em poucos
minutos. Sem WSL ou sem Docker a suite e pulada, com aviso, e sai com 0.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

for _s in (sys.stdout, sys.stderr):
    _reconf = getattr(_s, "reconfigure", None)
    if _reconf is not None:
        _reconf(encoding="utf-8", errors="replace")

RAIZ = Path(__file__).resolve().parents[1]
PIPELINE = RAIZ / "pipelines" / "filogenia"
VERDADEIRA = PIPELINE / "assets" / "test" / "arvore_verdadeira.nwk"
DISTRO = "Ubuntu-24.04"

falhas: list[str] = []


def ok(cond: bool, msg: str) -> None:
    print(f"  {'OK  ' if cond else 'FALHA'} {msg}")
    if not cond:
        falhas.append(msg)


def para_wsl(p: Path) -> str:
    """C:\\Users\\X -> /mnt/c/Users/X."""
    s = str(p.resolve()).replace("\\", "/")
    return f"/mnt/{s[0].lower()}{s[2:]}" if len(s) > 1 and s[1] == ":" else s


# --------------------------------------------------------------------------
# Newick e Robinson-Foulds, so com biblioteca padrao
# --------------------------------------------------------------------------

No = tuple[str, list["No"]]


def ler_newick(texto: str) -> No:
    s = texto.strip().rstrip(";").strip()
    pos = 0

    def no() -> No:
        nonlocal pos
        filhos: list[No] = []
        if s[pos] == "(":
            pos += 1
            while True:
                filhos.append(no())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"newick malformado perto da posicao {pos}")
        inicio = pos
        while pos < len(s) and s[pos] not in ",();":
            pos += 1
        return (s[inicio:pos].split(":")[0].strip(), filhos)

    return no()


def folhas(n: No) -> set[str]:
    nome, filhos = n
    if not filhos:
        return {nome}
    return set().union(*(folhas(f) for f in filhos))


def biparticoes(arv: No, taxa: set[str]) -> set[frozenset[str]]:
    """Biparticoes nao triviais, restritas a `taxa`, normalizadas pelo lado sem o taxon de referencia."""
    ref = min(taxa)
    out: set[frozenset[str]] = set()

    def rec(n: No) -> set[str]:
        nome, filhos = n
        if not filhos:
            return {nome} & taxa
        abaixo: set[str] = set()
        for f in filhos:
            abaixo |= rec(f)
        if 1 < len(abaixo) < len(taxa) - 1:
            out.add(frozenset(taxa - abaixo if ref in abaixo else abaixo))
        return abaixo

    rec(arv)
    return out


def rf(a: No, b: No) -> tuple[int, int]:
    """(distancia RF, maximo possivel) sobre os taxons em comum."""
    comuns = folhas(a) & folhas(b)
    return len(biparticoes(a, comuns) ^ biparticoes(b, comuns)), 2 * (len(comuns) - 3)


# --------------------------------------------------------------------------


def wsl_pronto() -> str | None:
    """None se WSL + docker + nextflow respondem; senao o motivo."""
    if not shutil.which("wsl"):
        return "wsl.exe fora do PATH"
    try:
        p = subprocess.run(
            [
                "wsl",
                "-d",
                DISTRO,
                "--",
                "bash",
                "-lc",
                "command -v nextflow docker && docker info -f ok",
            ],
            capture_output=True,
            text=True,
            timeout=90,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"WSL nao respondeu ({exc})"
    if p.returncode != 0:
        return f"nextflow ou docker indisponivel no WSL: {(p.stderr or p.stdout).strip()[:200]}"
    return None


def main() -> int:
    print("== pipeline de filogenia: perfil de teste ponta a ponta ==")
    motivo = wsl_pronto()
    if motivo:
        print(f"  PULADO: {motivo}")
        return 0

    saida = Path(tempfile.mkdtemp(prefix="filogenia_check_"))
    script = saida / "rodar.sh"
    # Script em arquivo e bash -l: o wsl.exe engole `$var` passado pela linha
    # de comando, e o PATH do nextflow so existe em shell de login.
    script.write_text(
        "set -uo pipefail\n"
        "mkdir -p ~/nf-work/filogenia-check && cd ~/nf-work/filogenia-check\n"
        f"nextflow run '{para_wsl(PIPELINE)}' -profile test,docker "
        f"--outdir '{para_wsl(saida / 'out')}' -w ~/nf-work/filogenia-check/work -resume\n",
        encoding="utf-8",
        newline="\n",
    )
    inicio = time.monotonic()
    proc = subprocess.run(
        ["wsl", "-d", DISTRO, "--", "bash", "-l", para_wsl(script)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=1800,
        check=False,
        env={**os.environ, "MSYS_NO_PATHCONV": "1"},
    )
    dur = time.monotonic() - inicio
    ok(proc.returncode == 0, f"nextflow terminou com exit 0 ({dur:.0f} s)")
    if proc.returncode != 0:
        print("\n".join(proc.stdout.splitlines()[-30:]))
        print("\n".join(proc.stderr.splitlines()[-30:]))
        return 1

    out = saida / "out"
    verdadeira = ler_newick(VERDADEIRA.read_text(encoding="utf-8"))

    resumo = json.loads((out / "resumo.json").read_text(encoding="utf-8"))
    ok(
        resumo.get("n_locos_entrada") == 5,
        f"5 locos na entrada (veio {resumo.get('n_locos_entrada')})",
    )
    ok(
        resumo.get("n_locos_inferidos") == 4,
        f"4 locos inferidos (veio {resumo.get('n_locos_inferidos')})",
    )
    ok(
        resumo.get("locos_fora") == {"loco_E": "poucos_taxa"},
        f"loco_E descartado por poucos taxons ({resumo.get('locos_fora')})",
    )

    concat = ler_newick(
        (out / "05_concatenado" / "concatenado.treefile").read_text(encoding="utf-8")
    )
    d, m = rf(concat, verdadeira)
    ok(d == 0, f"arvore concatenada = arvore verdadeira (RF {d}/{m})")

    especies = ler_newick((out / "07_arvore_especies" / "especies.tre").read_text(encoding="utf-8"))
    d, m = rf(especies, verdadeira)
    ok(d == 0, f"arvore de especies (wASTRAL) = arvore verdadeira (RF {d}/{m})")

    for arq in sorted((out / "04_arvores_genes").glob("*.treefile")):
        d, m = rf(ler_newick(arq.read_text(encoding="utf-8")), verdadeira)
        # Arvore de um loco so pode errar ramo curto; o criterio forte e o da
        # concatenada. Aqui so se exige que nenhuma esteja grosseiramente errada.
        ok(d <= 4, f"{arq.stem}: RF {d}/{m} contra a verdadeira (tolerancia 4)")

    cf = resumo.get("concordancia", {})
    ok(float(cf.get("gCF_medio", 0)) > 50, f"gCF medio > 50 (veio {cf.get('gCF_medio')})")
    ok("sCF_medio" in cf, f"sCFL calculado (veio {cf.get('sCF_medio')})")

    metodos = (out / "metodos.md").read_text(encoding="utf-8")
    ok("IQ-TREE v3.1.4" in metodos, "metodos.md cita IQ-TREE v3.1.4")
    # "A CONFIRMAR", sem o til: o marcador e "[VERSÃO A CONFIRMAR]", e a grafia
    # sem acento que estava aqui nunca casava — a verificacao passava sempre.
    ok("A CONFIRMAR" not in metodos, "nenhuma versao pendente em metodos.md")
    # Perfil de teste: versoes registradas, supermatriz em -T 1 e nenhum loco com
    # limpeza agressiva (todos retem > 80%) — as pendencias tem de vir vazias.
    pendencias = metodos.split("## Pendências", 1)[-1].strip()
    ok(pendencias == "- Nenhuma.", f"metodos.md sem pendencias (veio: {pendencias[:160]!r})")
    versoes = (out / "pipeline_info" / "versoes.yml").read_text(encoding="utf-8")
    ok("iqtree: 3.1.4" in versoes, "versoes.yml registra iqtree 3.1.4")

    print()
    if falhas:
        print(f"FALHOU: {len(falhas)} verificacao(oes). Saida mantida em {out}")
        return 1
    print(f"Pipeline de filogenia OK. Saida em {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
