#!/usr/bin/env bash
# bio-stack — camadas do uv no Windows: o .venv do projeto, os ambientes
# isolados de envs/ e a camada `uv tool` (anatomia, pose, a ponte do ChimeraX
# e o CLI do Virtual Cell Challenge).
#
# Roda no Git Bash. E trabalho de horas (o torch cu128 sozinho levou 3 h 19 min
# numa tentativa de 2026-10-21), entao vai por tarefa agendada, que sobrevive
# ao fim da sessao do Claude Code:
#
#   .\bootstrap\tarefa_longa.ps1 -Nome reconstruir-windows-uv -Script bootstrap\08_windows_uv.sh
#
# Idempotente: o que ja esta inteiro e conferido e mantido; so o que falta ou
# quebrou e refeito. O log vai para $1 (padrao runs/bootstrap/08_windows_uv.log).
#
# Nasceu do reinstalar_uv_tools_v3.sh, que viveu no scratchpad da sessao de
# 2026-10-22 e reconstruiu as quatro ferramentas com kernel sm_120 OK. O que
# cada passada ensinou:
#   21/10 16:57 — as cinco instalaram OK. As 19:16 um
#     `winget install astral-sh.uv` recriou AppData/Roaming/uv vazio e
#     levou tudo junto, inclusive o cpython-3.11 do .venv do projeto.
#   Primeira passada — TotalSegmentator/DeepLabCut/SLEAP cairam por rede
#     (torch cu128 esgotou 4 retries em 7674 s); nnInteractive e
#     chimerax-mcp cairam em "Executable already exists".
#   23/10 — o repositorio foi re-clonado e o .venv e os envs/*/.venv foram
#     junto; o vcc-cli entrou na lista. O 3.13 estava estripado porque um
#     colab-mcp orfao segurava o python.exe dele (ver passo 0).
#
# A ordem poupa download: o .venv vem primeiro e deixa o torch 2.11.0+cu128 no
# cache do uv, e BioEmu, TotalSegmentator, DeepLabCut e SLEAP reaproveitam o
# mesmo wheel. So o nnInteractive puxa outro (2.8.0).
#
# Defesas fixas:
#   --force em toda `uv tool install`, porque os trampolins orfaos em
#     ~/.local/bin (TotalSegmentator.exe, dlc.exe, sleap-*.exe) bloqueiam a
#     instalacao.
#   --python 3.13 nas uv tools: e o interpretador em que as seis estavam, e
#     sem ele o uv escolhe sozinho (o C:\Python314 do sistema esta no PATH).
#   uv por caminho absoluto, porque tarefa agendada nao herda o PATH do
#     shell interativo.
#
# Armadilhas de docs/STACK.md, todas aplicadas:
#   1. --torch-backend cu128, nunca o --index do cu128: aquele indice
#      tambem serve build de CPU e o CPU ganha a resolucao (medido:
#      instalou torch==2.14.1+cpu).
#   2. O extra nn-cuda128 do SLEAP so marca intencao, nao carrega indice.
#      Precisa vir acompanhado do --torch-backend.
#   3. totalseg_mcp nao declara uvicorn nem mcp e morre em
#      ModuleNotFoundError; entram por --with.

set -u

UV=/c/Users/Manec/.local/bin/uv.exe
PYTHONS=/c/Users/Manec/AppData/Roaming/uv/python
TOOLS=/c/Users/Manec/AppData/Roaming/uv/tools
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(dirname "$AQUI")"
LOG="${1:-$RAIZ/runs/bootstrap/08_windows_uv.log}"
mkdir -p "$(dirname "$LOG")"
: > "$LOG"

registrar() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }

passo() {
    local rotulo="$1"; shift
    registrar "=== $rotulo ==="
    if "$@" >>"$LOG" 2>&1; then
        registrar "OK   $rotulo"
    else
        registrar "FALHA $rotulo (codigo $?)"
    fi
}

# Kernel real, nao cuda.is_available(). Um wheel errado instala calado e so
# falha no primeiro kernel — e o `is_available()` mente nesse caso.
kernel() {
    local rotulo="$1" py="$2" saida
    if [ ! -x "$py" ]; then
        registrar "  $rotulo: sem interpretador (instalacao nao completou)"
        return
    fi
    saida=$("$py" -c "
import torch
a = torch.randn(256, 256, device='cuda')
(a @ a).sum().item()
torch.cuda.synchronize()
cc = torch.cuda.get_device_capability(0)
print('torch %s | cuda %s | kernel sm_%d%d OK' % (torch.__version__, torch.version.cuda, cc[0], cc[1]))
" 2>&1 | tail -1)
    registrar "  $rotulo: $saida"
}

# Interpretador gerenciado inteiro tem Lib/; o estripado fica com 6 entradas
# (DLLs/, python.exe, python3*.dll, vcruntime*). Sao exatamente os arquivos que
# um processo vivo trava: o uv apaga o resto e para no primeiro arquivo preso.
python_inteiro() {
    local d
    for d in "$PYTHONS"/cpython-"$1".*-windows-x86_64-none; do
        [ -d "$d/Lib" ] && return 0
    done
    return 1
}

registrar "inicio — camadas do uv no Windows"
registrar "uv: $("$UV" --version 2>&1)"

# --- 0. interpretadores -------------------------------------------------------
# 3.11: .venv e envs/. 3.13: todas as uv tools. Se a reinstalacao falhar com
# "Acesso negado (os error 5)", ha processo rodando daquele interpretador. Em
# 2026-10-23 era um colab-mcp orfao de sessao anterior; ache com
# Get-CimInstance Win32_Process | ? CommandLine -match 'cpython-3\.13'.
for v in 3.11 3.13; do
    if python_inteiro "$v"; then
        registrar "OK   Python $v (ja inteiro)"
    else
        passo "Python $v (ausente ou estripado)" "$UV" python install "$v" --reinstall
    fi
done

# --- 1. .venv do projeto ------------------------------------------------------
# O mesmo sync que bootstrap/02_python_env.ps1 faz sem flags: --locked exige
# que o uv.lock ja cubra o pyproject, e dev e gpu entram por serem
# default-groups (sem a camada gpu, o sync DESINSTALA o torch cu128).
passo ".venv (uv sync --locked, Python 3.11)" \
    "$UV" sync --project "$RAIZ" --python 3.11 --locked
kernel ".venv" "$RAIZ/.venv/Scripts/python.exe"

# --- 2. ambientes isolados ----------------------------------------------------
# A receita de cada um esta no cabecalho do proprio requirements.in. Venv que
# existe mas nao sobe e relatado, nao apagado: decidir isso e de quem le o log.
env_isolado() {
    local nome="$1"; shift
    local dir="$RAIZ/envs/$nome/.venv"
    if [ -x "$dir/Scripts/python.exe" ] && "$dir/Scripts/python.exe" -c "import sys" >/dev/null 2>&1; then
        registrar "OK   envs/$nome (ja existia)"
    elif [ -e "$dir" ]; then
        registrar "FALHA envs/$nome: $dir existe mas o python nao sobe — apague e rode de novo"
        return
    else
        passo "envs/$nome: venv" "$UV" venv "$dir" --python 3.11
    fi
    passo "envs/$nome: pacotes" \
        "$UV" pip install --python "$dir" "$@" -r "$RAIZ/envs/$nome/requirements.in"
}
env_isolado bioemu --torch-backend cu128
env_isolado bonej-headless
kernel "envs/bioemu" "$RAIZ/envs/bioemu/.venv/Scripts/python.exe"

# --- 3. uv tools --------------------------------------------------------------
ferramenta() {
    local rotulo="$1"; shift
    passo "$rotulo" "$UV" tool install --force --python 3.13 "$@"
}
ferramenta "TotalSegmentator" totalsegmentator --torch-backend cu128 --with uvicorn --with mcp
ferramenta "DeepLabCut" deeplabcut --torch-backend cu128
ferramenta "SLEAP" "sleap[nn-cuda128]" --torch-backend cu128
ferramenta "nnInteractive (servidor)" nninteractive-slicer-server --torch-backend cu128
ferramenta "chimerax-mcp" chimerax-mcp
# O nome do pacote saiu do trampolim orfao ~/.local/bin/vcc.exe, que aponta
# para uv/tools/vcc-cli. A chave da conta fica no chaveiro do Windows e
# sobrevive a reinstalacao.
ferramenta "vcc-cli (Virtual Cell Challenge)" vcc-cli

registrar "=== inventario ==="
"$UV" tool list >>"$LOG" 2>&1

registrar "=== CUDA: kernel real em sm_120 ==="
for t in totalsegmentator deeplabcut sleap nninteractive-slicer-server; do
    kernel "$t" "$TOOLS/$t/Scripts/python.exe"
done
# stdin fechado e teto de 60 s: sem conta no chaveiro, um prompt de login
# travaria a tarefa inteira em silencio.
registrar "  vcc: $(timeout 60 "$HOME/.local/bin/vcc.exe" whoami </dev/null 2>&1 | tail -1)"

registrar "fim"
