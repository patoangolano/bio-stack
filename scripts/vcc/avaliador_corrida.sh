#!/usr/bin/env bash
# scripts/vcc/avaliador_corrida.sh — a sequencia de pontuacao do cell-eval2, ponta a ponta.
#
# Etapa 5 de docs/VCC-2026-ROTA.md. Quatro passos, nessa ordem, e os quatro com o
# MESMO preset, o MESMO --pert-col e os MESMOS --set. Divergencia em qualquer um
# dos tres faz o `score --real-bundle` abortar, e nao existe waiver:
# `--allow-config-mismatch` ao lado de `--real-bundle` e erro de uso.
#
#   baseline          constroi o braco de resposta generica (a ponta 0 da escala)
#   prep-real-bundle  junta a ponta 0 e a ancora de replicata (a ponta 1) num bundle
#   run               pontua a predicao contra o real
#   score             compara os oito campos de identidade e normaliza as metricas
#
# Rodar do Windows, com `bash -l` (o PATH do venv e montado aqui, mas o `uv` e os
# shims vem de ~/.profile, lido so em shell de login):
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/vcc/avaliador_corrida.sh
#
# Variaveis de ambiente (todas com default):
#   CONTEXTO=D  N_PERTS=32  N_CTRL=18400  SEMENTE=0  TRABALHO=$HOME/vcc-avaliador
#
# O diretorio de trabalho fica em $HOME dentro da distro, NAO em /mnt/c: o h5ad
# intermediario e lido quatro vezes (uma por passo) e a ancora o le cinco vezes
# mais, uma por split.
set -euo pipefail

PROJETO=/mnt/c/Users/Manec/bio-stack
VENV="$HOME/envs/avaliador-proxy/.venv"
CONTEXTO="${CONTEXTO:-D}"
N_PERTS="${N_PERTS:-32}"
N_CTRL="${N_CTRL:-18400}"
SEMENTE="${SEMENTE:-0}"
TRABALHO="${TRABALHO:-$HOME/vcc-avaliador}"

MARCA="${CONTEXTO}-k${N_PERTS}-c${N_CTRL}-s${SEMENTE}"
PROXY="$TRABALHO/proxy-$MARCA"
OUT="$TRABALHO/out-$MARCA"
RELOGIO="$OUT/tempos.tsv"

[[ -x "$VENV/bin/cell-eval2" ]] || {
    echo "erro: $VENV/bin/cell-eval2 nao existe — rode envs/avaliador-proxy/montar.sh" >&2
    exit 1
}

rm -rf "$OUT"
mkdir -p "$OUT" "$PROXY"
printf 'passo\trss_kb\twall_s\n' > "$RELOGIO"

# `%M` e o pico de RSS em KiB do processo e dos filhos — a unica medida de
# memoria que vale aqui, porque o avaliador materializa a matriz uma vez por
# passo e o pico e o que decide se cabe.
cronometrar() {
    local passo="$1"; shift
    /usr/bin/time -f "${passo}\t%M\t%e" -a -o "$RELOGIO" "$@"
}

# ---------------------------------------------------------------- passo 0: dado
# Nao faz parte da sequencia do cell-eval2; monta o par (real, predito) que ela
# consome. Ver scripts/vcc/avaliador_proxy.py para o que e real e o que e proxy.
#
# `--pred-pool-unico` NAO e opcional aqui, e a razao esta medida no docstring de
# `_particionar`: sem ele o lado predito sorteia da mesma particao de ~460
# celulas-fonte que o lado real da mesma perturbacao, e o `pds_cosine` le esse
# rastro de identidade de celula como acerto. Medido com 40 perturbacoes e
# predicao de controle puro: 0,99936 de pds bruto com o rastro contra 0,41859
# sem, e +0,188 de avg_score contra -0,033. Uma predicao sintetica deste script
# e um braco de REFERENCIA, e referencia com rastro nao serve de referencia.
cronometrar proxy "$VENV/bin/python" -I "$PROJETO/scripts/vcc/avaliador_proxy.py" \
    --contexto "$CONTEXTO" --n-perts "$N_PERTS" --n-ctrl "$N_CTRL" \
    --semente "$SEMENTE" --pred-pool-unico --saida "$PROXY"

# Os tres argumentos que os quatro passos compartilham. `de.backend=pdex` fixado
# e nao `auto`: com `auto` o `cell-eval2` LEVANTA excecao num host com
# dispositivo CUDA visivel e sem motor de DE de GPU instalado, e num host sem
# CUDA cai para pdex com um aviso. Fixar transforma `resolved_de_backend` numa
# constante — e esse campo e um dos oito que o `score` compara.
COMUM=(--preset vcc2026 --pert-col target_gene --set de.backend=pdex)

# ------------------------------------------------------- passo 1: a ponta 0
# `--save-pred` e obrigatorio: o `prep-real-bundle` recebe o braco de baseline
# como PREDICAO (`--baseline`), nao o agregado dele.
cronometrar baseline "$VENV/bin/cell-eval2" baseline \
    -ar "$PROXY/real.h5ad" \
    "${COMUM[@]}" \
    --save-pred "$OUT/baseline_pred.h5ad" \
    -o "$OUT/baseline"

# ------------------------------------- passo 2: as duas pontas e a identidade
# Custa cinco corridas completas de metrica (uma por split da ancora) mais a do
# braco de baseline. E o passo caro, e o unico que se pode reaproveitar entre
# submissoes: enquanto o lado REAL nao muda, o bundle vale.
cronometrar prep_real_bundle "$VENV/bin/cell-eval2" prep-real-bundle \
    --real "$PROXY/real.h5ad" --baseline "$OUT/baseline_pred.h5ad" \
    "${COMUM[@]}" \
    --id "proxy-$MARCA" -o "$OUT/bundle" --force

# ------------------------------------------------------- passo 3: a submissao
cronometrar run "$VENV/bin/cell-eval2" run \
    -ap "$PROXY/pred.h5ad" -ar "$PROXY/real.h5ad" \
    "${COMUM[@]}" \
    -o "$OUT/run"

# ------------------------------------------------------- passo 4: a pontuacao
# `--real-bundle` e afirmacao: "isto e uma submissao de competicao". Qualquer um
# dos oito campos divergente aborta aqui.
cronometrar score "$VENV/bin/cell-eval2" score \
    --user-agg "$OUT/run/agg_results.csv" --real-bundle "$OUT/bundle" \
    -o "$OUT/scored.csv"

echo
echo "=== $MARCA ==="
cat "$RELOGIO"
echo
cat "$OUT/scored.csv"
echo
echo "artefatos: $OUT"
