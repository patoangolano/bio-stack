#!/usr/bin/env bash
# Interrompe so a corrida do LZ07 deste braco (as outras analises da maquina ficam).
pkill -f "f03_iqtree.sh LZ07" || true
pkill -f "matrizes/LZ07/LZ07.phy" || true
sleep 1
ps -eo pid,args | grep "LZ07" | grep -v grep | cut -c1-120 || true
