#!/usr/bin/env bash
# bio-stack — ambientes pixi de biomecanica no WSL (Ubuntu-24.04).
#
#   bash bootstrap/06_wsl_envs_biomec.sh      (dentro do WSL)
#
# Idempotente. Complementa 05_wsl_ferramentas.sh, que cuida da camada
# `pixi global` (uma CLI por ambiente). Aqui ficam os ambientes que NAO cabem
# naquele modelo, por dois motivos:
#
#   1. sao bibliotecas Python importadas por script de terceiros, nao CLIs, e
#      precisam do modulo no proprio interpretador;
#   2. travam versao de 2023 que briga com a camada global (o myoconverter pede
#      opensim 4.4.1 + numpy 1.21; o global tem opensim 4.6 + numpy 2.4).
#
# O manifesto e versionado em envs/<nome>/pixi.toml e e a fonte da verdade. O
# AMBIENTE, porem, e construido em ~/envs/<nome> dentro do WSL, nao em
# /mnt/c: um ambiente pixi sao dezenas de milhares de arquivos pequenos, e
# /mnt/c atravessa a fronteira de sistema de arquivos — ordens de grandeza mais
# lento, exatamente no padrao de acesso que mais sofre.
#
# O pixi.lock de cada um fica ao lado do ambiente, em ~/envs/<nome>, e e copiado
# de volta para envs/<nome>/ ao fim: e ele que fixa o commit do git do
# myoconverter e o fecho transitivo que a secao de Metodos cita.
set -euo pipefail

# Dependencia de sistema do ambiente opencap: `maskpass` -> `evdev`, que compila
# contra o header de input. O WSL nao tem linux-headers para o proprio kernel,
# mas evdev so precisa do header de API de usuario, em linux-libc-dev, que
# independe da versao do kernel. Idempotente.
#
# Quem instala e o 04b_wsl_base.sh, como root. Aqui so `sudo -n`: em
# 2026-10-16 um `sudo` interativo ficou 19 min esperando a senha de um usuario
# que nao tinha senha, sem erro e sem saida. Melhor falhar na hora.
if [ ! -f /usr/include/linux/input-event-codes.h ]; then
  echo "== instalando linux-libc-dev (header de input, exigido por evdev) =="
  sudo -n apt-get install -y -q linux-libc-dev || {
    echo "sem sudo sem senha; rode como root:" >&2
    echo "  wsl -d Ubuntu-24.04 -u root -- bash /mnt/c/Users/Manec/bio-stack/bootstrap/04b_wsl_base.sh" >&2
    exit 1
  }
fi

projeto="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
destino="$HOME/envs"
mkdir -p "$destino"

for nome in myoconverter opencap; do
  manifesto="$projeto/envs/$nome/pixi.toml"
  if [ ! -f "$manifesto" ]; then
    echo "SEM MANIFESTO: $manifesto" >&2
    exit 1
  fi

  mkdir -p "$destino/$nome"
  cp "$manifesto" "$destino/$nome/pixi.toml"
  # Reaproveita o lock anterior, se houver, para nao re-resolver a cada execucao.
  if [ -f "$projeto/envs/$nome/pixi.lock" ] && [ ! -f "$destino/$nome/pixi.lock" ]; then
    cp "$projeto/envs/$nome/pixi.lock" "$destino/$nome/pixi.lock"
  fi

  echo "== $nome =="
  ( cd "$destino/$nome" && pixi install )
  cp "$destino/$nome/pixi.lock" "$projeto/envs/$nome/pixi.lock"
  echo "   ambiente: $destino/$nome    lock regravado em envs/$nome/pixi.lock"
done

echo
echo "Como usar (shell de login do WSL):"
echo "  pixi run --manifest-path ~/envs/myoconverter/pixi.toml python -c 'import myoconverter'"
echo "  pixi run --manifest-path ~/envs/opencap/pixi.toml      python -c 'import opensim'"
