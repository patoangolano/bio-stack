#!/usr/bin/env bash
# scripts/vcc/montar_v7.sh — emite, valida e empacota a submissao v7.
#
# v7 = `k_desvio = 4` com `alpha = 1`, o encolhimento POR GENE medido contra
# verdade retida (docs/VCC-2026-EXECUCAO.md secao 20). Vale +0,0266 em
# `avg_score` contra a v6, e nenhum dos seis membros perde.
#
# Rodar pelo agendador, nao pelo Bash tool: a emissao sao ~11 min e job de
# duracao iniciado pelo Bash morre de SIGHUP quando a sessao fecha.
#
#   .\bootstrap\tarefa_longa.ps1 -Nome montar-v7 -Script scripts\vcc\montar_v7.sh
#
# NAO submete. O envio e voltado para fora, gasta um dos dois slots do dia e
# tem de ser decidido DEPOIS de ler os laudos — entao ele fica fora deste
# script, de proposito. O que sai daqui e um `.vcc` pronto e tres laudos.
set -euo pipefail

RAIZ=/c/Users/Manec/bio-stack
PY="$RAIZ/.venv/Scripts/python.exe"
SUB="$RAIZ/runs/vcc/submissao"
REF="$RAIZ/data/raw/vcc/controls-final"
DELTA="$RAIZ/runs/vcc/transfer/varredura_k/delta_kd4_a1.npz"
H5AD="$SUB/submissao_v7.h5ad"
VCC="$SUB/submissao_v7.vcc"

marco() { printf '\n========== %s  (%s)\n' "$1" "$(date +%H:%M:%S)"; }

[[ -r "$DELTA" ]] || { echo "erro: $DELTA nao existe — rode varredura_k.py --k 4 --alphas 1" >&2; exit 1; }
mkdir -p "$SUB"

# O gerador tem gate de disco proprio (25 GB), mas falhar DEPOIS de 11 min de
# emissao e pior que falhar antes.
LIVRE=$(df -BG --output=avail /c | tail -1 | tr -dcs '0-9' '\n' | head -1)
marco "disco livre: ${LIVRE} GB"
[[ "$LIVRE" -ge 25 ]] || { echo "erro: ${LIVRE} GB livres, preciso de 25" >&2; exit 1; }

marco "1/4 emissao dos tres contextos"
"$PY" -I "$RAIZ/scripts/vcc/gerar_submissao.py" "$DELTA" \
    --saida "$H5AD" --relatorio "$SUB/laudo_v7.json"

marco "2/4 validador de CONTRATO (28 itens)"
"$PY" -I "$RAIZ/scripts/vcc/validar_submissao.py" "$H5AD" \
    --ref "$REF" --json "$SUB/contrato_v7.json"

marco "3/4 validador NUMERICO (inclinacao contra o delta pretendido)"
"$PY" -I "$RAIZ/scripts/vcc/validar_numerico.py" "$H5AD" \
    --delta "$DELTA" --controles "$REF" --json "$SUB/numerico_v7.json"

marco "4/4 vcc prep"
# O `prep` e do cpython-3.13 das uv tools e nao do venv do projeto; `vcc` esta
# no PATH do shell de LOGIN, que e o que a tarefa_longa usa.
vcc prep "$H5AD" --output "$VCC"

marco "pronto"
ls -la "$H5AD" "$VCC" | awk '{print $5, $9}'
echo
echo "laudos: laudo_v7.json contrato_v7.json numerico_v7.json"
# Aspas SIMPLES: entre aspas duplas o par de crases viraria substituicao de
# comando e o script EXECUTARIA `vcc submit` ao terminar. `bash -n` nao pega,
# porque e sintaxe valida.
echo 'o ENVIO nao acontece aqui: leia os laudos e rode `vcc submit` a mao'
