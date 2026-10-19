#!/usr/bin/env bash
# Sonda do ambiente WSL para o braco paleontologico (fosseis) dos rapinantes.
set -uo pipefail
command -v iqtree3
iqtree3 --version | head -2
nproc
uptime
ps -eo pid,pcpu,pmem,etime,comm --sort=-pcpu | head -8
df -h ~ | tail -1
