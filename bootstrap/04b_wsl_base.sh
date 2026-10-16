#!/usr/bin/env bash
# bio-stack — camada base do WSL (Ubuntu-24.04): usuario, Java, bedtools,
# Docker, Nextflow, nf-core, uv e as dependencias de sistema da sandbox do
# Claude Science.
#
# Roda como root, a partir do Windows. O `-u root` dispensa senha e pula o
# assistente de primeiro uso da distro:
#
#   wsl -d Ubuntu-24.04 -u root -- bash /mnt/c/Users/Manec/bio-stack/bootstrap/04b_wsl_base.sh
#
# Depois dele, como o usuario: 05_wsl_ferramentas.sh (pixi + CLIs, inclusive
# a biofisica) e 06_wsl_envs_biomec.sh (ambientes de biomecanica).
#
# Era a unica camada sem script. Em 2026-10-10 foi montada a mao; quando a
# distro sumiu, em 2026-10-16, a reconstrucao dependeu de garimpar comandos no
# runs/provenance.jsonl. Idempotente: cada passo confere antes de agir.
#
# O que NAO instala, de proposito: samtools, bcftools e tabix do apt. O htslib
# 1.19 do Ubuntu 24.04 tem CVEs abertas no decodificador CRAM (corrigidas na
# 1.23.1); quem responde por eles e o ambiente `hts` do pixi (1.24).
set -euo pipefail

USUARIO=manec
NXF_VER=26.04.6
NFCORE_VER=4.1.0

if [ "$(id -u)" -ne 0 ]; then
    echo "rode como root: wsl -d Ubuntu-24.04 -u root -- bash $0" >&2
    exit 1
fi

# runuser -l abre shell de login: e o que le o ~/.profile e poe ~/.local/bin
# no PATH (armadilha 3 do CLAUDE.md).
como_usuario() { runuser -l "$USUARIO" -c "$1"; }

# --- 1. usuario -----------------------------------------------------------------
# uid 1000 e o que o /etc/wsl-distribution.conf declara como usuario padrao.
if id "$USUARIO" >/dev/null 2>&1; then
    echo "usuario $USUARIO ja existe (uid $(id -u "$USUARIO"))"
else
    useradd -m -u 1000 -s /bin/bash -G sudo "$USUARIO"
    echo "usuario $USUARIO criado SEM senha; defina com:  wsl -d Ubuntu-24.04 -u root passwd $USUARIO"
fi
if grep -q '^\[user\]' /etc/wsl.conf 2>/dev/null; then
    sed -i "s/^default=.*/default=$USUARIO/" /etc/wsl.conf
else
    printf '\n[user]\ndefault=%s\n' "$USUARIO" >> /etc/wsl.conf
fi
echo "/etc/wsl.conf: usuario padrao = $USUARIO (vale depois de wsl --terminate Ubuntu-24.04)"

# --- 2. apt ---------------------------------------------------------------------
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
# linux-libc-dev: header de input que o evdev do ambiente opencap (06) compila.
# bubblewrap e socat: a sandbox do Claude Science (07) nao sobe sem eles.
apt-get install -y -qq openjdk-21-jre-headless bedtools pipx curl ca-certificates git unzip \
    linux-libc-dev bubblewrap socat >/dev/null
echo "apt: $(java -version 2>&1 | head -1) | $(bedtools --version)"

# --- 2b. uv em caminho de sistema ---------------------------------------------------
# Os conectores locais do Claude Science (servers/*) rodam por `uv run --script`
# dentro da sandbox, que nao enxerga ~/.local/bin. Versao fixa no URL do instalador.
UV_VER=0.13.1
if /usr/local/bin/uv --version 2>/dev/null | grep -q " $UV_VER "; then
    echo "uv $UV_VER ja instalado em /usr/local/bin"
else
    uv_inst="$(mktemp)"
    curl -fsSL --retry 3 -o "$uv_inst" "https://astral.sh/uv/$UV_VER/install.sh"
    echo "instalador do uv $UV_VER: sha256 $(sha256sum "$uv_inst" | cut -d' ' -f1)"
    env UV_INSTALL_DIR=/usr/local/bin UV_NO_MODIFY_PATH=1 sh "$uv_inst" >/dev/null 2>&1
    rm -f "$uv_inst"
    echo "uv: $(/usr/local/bin/uv --version)"
fi

# --- 3. Docker nativo: docker-ce do repositorio oficial, sem Docker Desktop -------
if ! command -v docker >/dev/null 2>&1; then
    install -m 0755 -d /etc/apt/keyrings
    curl -fsSL --retry 3 https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
    chmod a+r /etc/apt/keyrings/docker.asc
    # shellcheck disable=SC1091
    . /etc/os-release
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu ${UBUNTU_CODENAME:-$VERSION_CODENAME} stable" \
        > /etc/apt/sources.list.d/docker.list
    apt-get update -qq
    apt-get install -y -qq docker-ce docker-ce-cli containerd.io docker-buildx-plugin >/dev/null
fi
usermod -aG docker "$USUARIO"
if [ "$(ps -p 1 -o comm=)" = systemd ]; then
    systemctl enable --now docker >/dev/null 2>&1
    echo "docker: $(docker --version) | daemon $(systemctl is-active docker)"
else
    systemctl enable docker >/dev/null 2>&1
    echo "AVISO: systemd nao e o PID 1; o daemon do Docker so sobe depois de wsl --terminate Ubuntu-24.04" >&2
fi

# --- 4. Nextflow, versao fixa, em ~/.local/bin -----------------------------------
# O get.nextflow.io E o proprio launcher, com a versao corrente embutida; ele so
# se instala quando chega por pipe (`$0 = bash`), e o guard de dados barra
# script encanado para o shell. Aqui baixa-se direto o launcher da versao
# fixada — o mesmo arquivo que o modo pipe buscaria — e `-download` traz o jar
# para ~/.nextflow/framework. O sha256 do launcher fica no log.
if como_usuario "command -v nextflow >/dev/null && nextflow -version 2>/dev/null | grep -q 'version $NXF_VER'"; then
    echo "nextflow $NXF_VER ja instalado"
else
    launcher="$(mktemp)"
    curl -fsSL --retry 3 -o "$launcher" "https://www.nextflow.io/releases/v$NXF_VER/nextflow"
    grep -q "^NXF_VER=\${NXF_VER:-'$NXF_VER'}" "$launcher" || { echo "launcher baixado nao e o $NXF_VER" >&2; exit 1; }
    echo "launcher do Nextflow $NXF_VER: sha256 $(sha256sum "$launcher" | cut -d' ' -f1)"
    como_usuario "mkdir -p ~/.local/bin"
    install -m 0755 -o "$USUARIO" -g "$USUARIO" "$launcher" "/home/$USUARIO/.local/bin/nextflow"
    rm -f "$launcher"
    como_usuario "nextflow -download >/dev/null"
fi

# --- 5. nf-core por pipx ---------------------------------------------------------
if como_usuario "pipx list --short 2>/dev/null | grep -q '^nf-core $NFCORE_VER\$'"; then
    echo "nf-core $NFCORE_VER ja instalado"
else
    como_usuario "pipx install --force nf-core==$NFCORE_VER >/dev/null"
fi

# --- 6. resumo -------------------------------------------------------------------
como_usuario "echo '== camada base =='; nextflow -version 2>/dev/null | grep version; nf-core --version 2>/dev/null | tail -1; docker --version"
echo
echo "proximo: wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/05_wsl_ferramentas.sh"
