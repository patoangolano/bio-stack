#!/usr/bin/env bash
set -eu
W=$HOME/vcc-avaliador/tres
for f in tempos_varredura.tsv tempos_varredura2.tsv tempos_varredura3.tsv; do
    [ -f "$W/$f" ] && tail -n +2 "$W/$f" >> "$W/tempos.tsv"
done
V=$HOME/envs/avaliador-proxy/.venv
"$V/bin/python" -I /mntC:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/b9d649df-1ccb-4bb9-949b-094924f36cea/scratchpad/recuperado/scripts/vcc/avaliador_verificacao.py \
    --corrida "k32=$HOME/vcc-avaliador/out-D-k32-c18400-s0" \
    --proxy "k32=$HOME/vcc-avaliador/proxy-D-k32-c18400-s0" \
    --peers "$HOME/vcc-avaliador/peers-k4.json" \
    --python-venv "$V/bin/python" \
    --comparacao "$W" \
    --bracos "oraculo_unico,controle_unico,nossa_v4,kd1.5,kd2,kd3,kd4,kd5,kd6,kd8,nossa"
