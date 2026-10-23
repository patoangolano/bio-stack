#!/usr/bin/env bash
# bio-stack — cria o remoto `gdrive` do rclone, que o backup_drive.sh usa.
#
# Roda DENTRO do WSL e num terminal seu: o rclone imprime uma URL local e
# espera o consentimento voltar do navegador em 127.0.0.1:53682. E o unico
# passo do backup que pede mao humana.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/rclone_gdrive.sh
#
# Escopo drive.file: o rclone so enxerga o que ele mesmo subiu, nao o resto do
# Drive. Consequencia medida em 2026-10-19: o que subiu por outro caminho (o
# nf-work, pelo G:) fica invisivel a ele. Para restaurar, o restaura_drive.sh
# le do G:, que ve tudo.
#
# O rclone vem do manifesto (envs.rclone em bootstrap/pixi-global.toml), nao de
# `pixi global install` avulso: o 05 copia o manifesto versionado por cima do
# global e roda `pixi global sync`, que apagaria um ambiente instalado a mao.
# Restaurado do instala_rclone.sh + autoriza_rclone.sh do scratchpad de
# 2026-10-19.

set -euo pipefail

command -v rclone >/dev/null || { echo "rclone ausente — rode bootstrap/05_wsl_ferramentas.sh" >&2; exit 1; }
rclone version | head -1

if rclone listremotes | grep -qx 'gdrive:'; then
    echo "remoto gdrive ja existe"
else
    rclone config create gdrive drive scope=drive.file
fi
rclone about gdrive:
