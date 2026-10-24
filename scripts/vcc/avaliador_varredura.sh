#!/usr/bin/env bash
# scripts/vcc/avaliador_varredura.sh — pontua VARIOS bracos contra UM bundle real.
#
# O `avaliador_corrida.sh` faz a sequencia inteira para um braco sintetico. Este
# faz a parte que interessa para escolher `k_desvio`: monta o lado real uma vez,
# com rotulos do painel VCC (`--so-alvos-vcc`), constroi o bundle uma vez, e
# entao, para cada submissao passada, recorta e pontua.
#
# A economia nao e detalhe. Medido na corrida de 2026-10-24 (32 alvos, 18.400
# controles): proxy 9 s, baseline 94 s, prep-real-bundle **444 s**, run 99 s,
# score 1,4 s. O bundle e 76% do custo e depende SO do lado real, que nao muda
# entre bracos — reusa-lo leva o custo de um braco de ~11 min para ~2 min.
#
# Rodar do Windows:
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/vcc/avaliador_varredura.sh \
#       /mnt/c/Users/Manec/bio-stack/runs/vcc/submissao/sub_kd1.h5ad:k1 \
#       /mnt/c/Users/Manec/bio-stack/runs/vcc/submissao/sub_kd4.h5ad:k4
#
# Cada argumento e `caminho.h5ad:rotulo`. O rotulo nomeia o `scored_<rotulo>.csv`.
#
# ⚠️ O que isto mede, e o que nao mede: a predicao e nossa de verdade, mas a
# "verdade" e o efeito medido pela Arc em ARC_H1 (celula H1 embrionaria) aplicado
# a celulas de D. A pontuacao responde "quanto a predicao concorda com uma
# resposta de formato H1", e serve para COMPARAR bracos no mesmo proxy. Nao
# estima o placar. Ver o docstring de avaliador_recorte.py.
set -euo pipefail

PROJETO=/mnt/c/Users/Manec/bio-stack
VENV="$HOME/envs/avaliador-proxy/.venv"
CONTEXTO="${CONTEXTO:-D}"
N_PERTS="${N_PERTS:-40}"
N_CTRL="${N_CTRL:-18400}"
SEMENTE="${SEMENTE:-0}"
# orion/hct116 e nao arc_vcc_25: com `--so-alvos-vcc` o arc tem 12 alvos do
# painel e so 11 passam o gate de 200 genes moventes, de 32 pedidos. O orion
# tem 148. A profundidade do orion e rasa (mediana de 2 genes significativos
# por perturbacao contra 2.016 do arc), mas aqui ele e GERADOR do lado real, e
# nao doador de efeito para a predicao: o que o avaliador precisa e de
# estrutura de DE com os rotulos certos, e e o numero de alvos que manda.
FONTE_EFEITOS="${FONTE_EFEITOS:-orion}"
CONTEXTO_PIE="${CONTEXTO_PIE:-hct116}"
TRABALHO="${TRABALHO:-$HOME/vcc-avaliador}"

MARCA="vcc-${CONTEXTO}-${FONTE_EFEITOS}-k${N_PERTS}-s${SEMENTE}"
PROXY="$TRABALHO/proxy-$MARCA"
OUT="$TRABALHO/out-$MARCA"
RELOGIO="$OUT/tempos.tsv"

[[ $# -ge 1 ]] || { echo "uso: $0 <submissao.h5ad:rotulo> [...]" >&2; exit 2; }
[[ -x "$VENV/bin/cell-eval2" ]] || {
    echo "erro: $VENV/bin/cell-eval2 nao existe — rode envs/avaliador-proxy/montar.sh" >&2
    exit 1
}

mkdir -p "$OUT" "$PROXY"
[[ -r "$RELOGIO" ]] || printf 'passo\trss_kb\twall_s\n' > "$RELOGIO"
cronometrar() {
    local passo="$1"; shift
    /usr/bin/time -f "${passo}\t%M\t%e" -a -o "$RELOGIO" "$@"
}

COMUM=(--preset vcc2026 --pert-col target_gene --set de.backend=pdex)

# --------------------------------------------------- o lado real, UMA vez
# `--so-alvos-vcc` e o que torna o recorte possivel: sem ele as perturbacoes do
# proxy levam nomes de genes que o painel VCC nao tem, e a nossa submissao nao
# cobre nenhuma delas.
if [[ ! -r "$PROXY/real.h5ad" ]]; then
    cronometrar proxy "$VENV/bin/python" -I "$PROJETO/scripts/vcc/avaliador_proxy.py" \
        --contexto "$CONTEXTO" --fonte-efeitos "$FONTE_EFEITOS" \
        --contexto-pie "$CONTEXTO_PIE" \
        --n-perts "$N_PERTS" --n-ctrl "$N_CTRL" --semente "$SEMENTE" \
        --so-alvos-vcc --pred-pool-unico --saida "$PROXY"
else
    echo "== lado real reusado: $PROXY/real.h5ad"
fi

if [[ ! -r "$OUT/bundle/bundle.json" && ! -d "$OUT/bundle/splits" ]]; then
    cronometrar baseline "$VENV/bin/cell-eval2" baseline \
        -ar "$PROXY/real.h5ad" "${COMUM[@]}" \
        --save-pred "$OUT/baseline_pred.h5ad" -o "$OUT/baseline"
    cronometrar prep_real_bundle "$VENV/bin/cell-eval2" prep-real-bundle \
        --real "$PROXY/real.h5ad" --baseline "$OUT/baseline_pred.h5ad" \
        "${COMUM[@]}" --id "proxy-$MARCA" -o "$OUT/bundle" --force
else
    echo "== bundle reusado: $OUT/bundle"
fi

# ------------------------------------------------- um braco por submissao
for ARGUMENTO in "$@"; do
    CAMINHO="${ARGUMENTO%:*}"
    ROTULO="${ARGUMENTO##*:}"
    [[ -r "$CAMINHO" ]] || { echo "erro: $CAMINHO nao existe" >&2; exit 1; }
    echo
    echo "== braco $ROTULO  <-  $CAMINHO"
    cronometrar "recorte_$ROTULO" "$VENV/bin/python" -I \
        "$PROJETO/scripts/vcc/avaliador_recorte.py" \
        --submissao "$CAMINHO" --proxy-real "$PROXY/real.h5ad" \
        --contexto "$CONTEXTO" --saida "$OUT/pred_$ROTULO.h5ad" \
        --resumo "$OUT/recorte_$ROTULO.json"
    cronometrar "run_$ROTULO" "$VENV/bin/cell-eval2" run \
        -ap "$OUT/pred_$ROTULO.h5ad" -ar "$PROXY/real.h5ad" \
        "${COMUM[@]}" -o "$OUT/run_$ROTULO"
    cronometrar "score_$ROTULO" "$VENV/bin/cell-eval2" score \
        --user-agg "$OUT/run_$ROTULO/agg_results.csv" --real-bundle "$OUT/bundle" \
        -o "$OUT/scored_$ROTULO.csv"
    # O pred recortado tem 31.200 celulas e nao e reusavel; o h5ad fica, mas a
    # pasta de `run` carrega a tabela de DE inteira e e ela que enche o disco.
    rm -rf "$OUT/run_$ROTULO/de_"*
done

echo
echo "=== $MARCA ==="
cat "$RELOGIO"
for CSV in "$OUT"/scored_*.csv; do
    [[ -e "$CSV" ]] || continue
    echo
    echo "--- $(basename "$CSV")"
    cat "$CSV"
done
echo
echo "artefatos: $OUT"
