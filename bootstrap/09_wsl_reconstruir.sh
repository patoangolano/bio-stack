#!/usr/bin/env bash
# bio-stack — reconstroi a distro WSL inteira, da instalacao ao daemon do
# Claude Science, sem nenhum passo que espere teclado.
#
# Roda no Git Bash do WINDOWS, nao dentro do WSL: o primeiro passo e criar a
# distro. Sao horas de trabalho (41+ ambientes pixi, GROMACS, PyMOL, OpenMM,
# OpenSim, os envs de biomecanica), entao vai por tarefa agendada:
#
#   .\bootstrap\tarefa_longa.ps1 -Nome reconstruir-wsl -Script bootstrap\09_wsl_reconstruir.sh -EsperarPor reconstruir-windows-uv
#
# Nasceu do reconstruir_wsl_pixi.sh do scratchpad de 2026-10-22, que partia de
# uma camada base ja de pe. Este parte do zero, como depois do
# `wsl --unregister Ubuntu-24.04` de 2026-10-23, e segue a ordem do CLAUDE.md:
# instalar, 04b como root, terminate, 05, 06, 07. Se a distro ja existe, pula a
# instalacao e reaplica as camadas (todas sao idempotentes).
#
# As tres armadilhas do CLAUDE.md, todas presentes e todas necessarias:
#   1. wsl.exe engole $var na traducao da linha de comando — por isso cada
#      passo chama um ARQUIVO de script, nunca um comando com cifrao.
#   2. Git Bash reescreve /mnt/c/... para caminho Windows e cola a raiz do
#      Git na frente — por isso MSYS_NO_PATHCONV=1, e cygpath -m onde o
#      argumento e para programa do Windows.
#   3. bash -l, nunca bash: o PATH de ~/.local/bin, dos shims do pipx e do
#      ~/.pixi/bin vive em ~/.profile, lido so em shell de LOGIN. Sem o -l a
#      falha e parcial e enganosa: nextflow e as CLIs do pixi dao command not
#      found enquanto bedtools, Java e docker seguem de pe em /usr/bin.
#
# Sobre sudo: o 04b cria `manec` SEM senha e ja instala como root o que o 06
# pediria por `sudo -n` (linux-libc-dev). Nenhum passo aqui espera senha — o
# travamento de 19 min de 2026-10-16 foi um sudo interativo.
#
# Docker Desktop: com a WSL integration ligada para esta distro, ele poe o
# docker e o socket dele dentro dela e disputa /var/run/docker.sock com o
# docker-ce do 04b. Desligue a integracao (ou desinstale o Docker Desktop)
# antes de rodar isto.

set -u
export MSYS_NO_PATHCONV=1

DISTRO=Ubuntu-24.04
UV=/c/Users/Manec/.local/bin/uv.exe
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(dirname "$AQUI")"
RAIZ_WIN="$(cygpath -m "$RAIZ")"
BOOT=/mnt/c/Users/Manec/bio-stack/bootstrap
LOG="${1:-$RAIZ/runs/bootstrap/09_wsl_reconstruir.log}"
mkdir -p "$(dirname "$LOG")"
: > "$LOG"

registrar() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }

passo() {
    local rotulo="$1"; shift
    registrar "=== $rotulo ==="
    if "$@" >>"$LOG" 2>&1; then
        registrar "OK   $rotulo"
        return 0
    else
        local codigo=$?
        registrar "FALHA $rotulo (codigo $codigo)"
        return "$codigo"
    fi
}

# wsl -l escreve em UTF-16; sem tirar os NUL o grep nunca casa.
distro_existe() { wsl.exe -l -q 2>/dev/null | tr -d '\0\r' | grep -qx "$DISTRO"; }

registrar "inicio — distro $DISTRO do zero"

if distro_existe; then
    registrar "OK   $DISTRO ja registrada — pula a instalacao"
else
    passo "instalar $DISTRO (web-download, sem abrir)" \
        wsl.exe --install -d "$DISTRO" --no-launch --web-download
fi
if ! distro_existe; then
    registrar "FALHA $DISTRO nao aparece em wsl -l depois da instalacao — parando"
    exit 1
fi

# Sem a camada base nao ha usuario manec, pixi, Java nem docker: o resto so
# empilharia falhas. Aqui para.
passo "04b — camada base (root)" \
    wsl.exe -d "$DISTRO" -u root -- bash "$BOOT/04b_wsl_base.sh" || exit 1
# O usuario padrao que o 04b grava em /etc/wsl.conf so vale na proxima partida.
passo "terminate (usuario padrao do 04b passa a valer)" wsl.exe --terminate "$DISTRO"

# Daqui em diante uma falha nao para a fila: cada camada e independente o
# bastante para valer a pena tentar, e o log diz o que faltou.
passo "05 — pixi: CLIs + GROMACS, PyMOL, OpenMM, OpenSim" \
    wsl.exe -d "$DISTRO" -- bash -l "$BOOT/05_wsl_ferramentas.sh"
passo "06 — envs de biomecanica: MyoConverter, opencap-processing" \
    wsl.exe -d "$DISTRO" -- bash -l "$BOOT/06_wsl_envs_biomec.sh"
passo "07 — daemon do Claude Science na 8765" \
    wsl.exe -d "$DISTRO" -- bash -l "$BOOT/07_wsl_claude_science.sh"
# Um `wsl` sem -d caia na docker-desktop, que virou a padrao quando o Docker
# Desktop entrou em 2026-10-21.
passo "distro padrao do WSL = $DISTRO" wsl.exe --set-default "$DISTRO"

passo "conferencia: tests/check_wsl_ferramentas.py" \
    "$UV" run --script "$RAIZ_WIN/tests/check_wsl_ferramentas.py"

registrar "=== falta a mao ==="
registrar "  senha do manec (o sudo so vale depois): wsl -d $DISTRO -u root passwd manec"
registrar "  login do Claude Science (link de uso unico): wsl -d $DISTRO -- bash -lc \"claude-science url\""
registrar "fim"
