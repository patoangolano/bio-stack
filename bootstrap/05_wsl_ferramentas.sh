#!/usr/bin/env bash
# bio-stack — camada de ferramentas de linha de comando no WSL (Ubuntu-24.04).
#
# Rodar de dentro do WSL, ou do Windows com as tres defesas do CLAUDE.md:
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/bootstrap/05_wsl_ferramentas.sh
#
# O que faz, e e idempotente:
#   1. instala o pixi (binario oficial, versao e sha256 fixados — nada de curl | sh);
#   2. poe ~/.pixi/bin no PATH dentro do bloco bio-stack do ~/.profile;
#   3. aplica o manifesto versionado bootstrap/pixi-global.toml (pixi global sync);
#   4. grava bootstrap/pixi-global.snapshot.tsv com versao, build, canal e
#      sha256 de cada pacote declarado no manifesto (sem as dependencias
#      transitivas) — e o que a secao de Metodos cita;
#   5. baixa a base do TaxonKit (taxdump do NCBI, md5 conferido) se faltar, e
#      grava a data do download em ~/.taxonkit/BAIXADO_EM.
#
# Por que pixi global e nao apt nem micromamba: cada ferramenta (ou grupo
# coeso, como hts e kent) ganha um ambiente proprio e nao divide dependencias
# com as outras. Os 43 pacotes ate resolvem num ambiente so (pixi 0.81, ~2 s),
# mas ai o solver troca o que e compartilhado: o perl desce de 5.44 para
# 5.32.1, clipkit e orthofinder ficam em python 3.13 e o bcftools 1.24 muda
# de build — os builds exatos do snapshot nao cabem juntos. Alem disso,
# o manifesto e declarativo e fixa versao exata, e os binarios saem num
# diretorio so. As do apt (samtools, bcftools, bedtools, tabix) nao sao
# removidas — analises antigas registraram o 1.19 —, mas o ambiente `hts` do
# manifesto (samtools/bcftools/htslib 1.24) faz sombra a elas em shell de
# login, porque ~/.pixi/bin entra antes de /usr/bin. Motivo: o htslib 1.19 do
# Ubuntu 24.04 tem CVEs abertas no decodificador CRAM (corrigidas na 1.23.1).
set -euo pipefail

PIXI_VERSAO="0.81.0"
PIXI_SHA256="7aa3ec39aecceff9062fa2ed4d42cbaa0bdc25ddea727d048e061cf188d434f6"
AQUI="$(cd "$(dirname "$0")" && pwd)"
MANIFESTO="$AQUI/pixi-global.toml"
SNAPSHOT="$AQUI/pixi-global.snapshot.tsv"

# --- 1. pixi -----------------------------------------------------------------
if [ -x "$HOME/.pixi/bin/pixi" ] && "$HOME/.pixi/bin/pixi" --version | grep -q " $PIXI_VERSAO\$"; then
    echo "pixi $PIXI_VERSAO ja instalado"
else
    tmp="$(mktemp -d)"
    trap 'rm -rf "$tmp"' EXIT
    url="https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSAO}/pixi-x86_64-unknown-linux-musl.tar.gz"
    echo "baixando $url"
    curl -fsSL --retry 3 -o "$tmp/pixi.tar.gz" "$url"
    echo "$PIXI_SHA256  $tmp/pixi.tar.gz" | sha256sum -c -
    mkdir -p "$HOME/.pixi/bin"
    tar -xzf "$tmp/pixi.tar.gz" -C "$tmp"
    install -m 0755 "$tmp/pixi" "$HOME/.pixi/bin/pixi"
    "$HOME/.pixi/bin/pixi" --version
fi

# --- 2. PATH no bloco bio-stack do ~/.profile ---------------------------------
# So o ~/.profile: e ele que o `bash -l` le. O ~/.bashrc do Ubuntu desiste em
# shell nao interativa (armadilha 3 do CLAUDE.md).
LINHA='export PATH="$HOME/.pixi/bin:$PATH"'
if ! grep -qF "$LINHA" "$HOME/.profile"; then
    if grep -q '^# <<< bio-stack PATH <<<' "$HOME/.profile"; then
        sed -i "/^# <<< bio-stack PATH <<</i $LINHA" "$HOME/.profile"
    else
        printf '\n# >>> bio-stack PATH >>>\n%s\n# <<< bio-stack PATH <<<\n' "$LINHA" >> "$HOME/.profile"
    fi
    echo "PATH do pixi acrescentado ao ~/.profile"
fi
export PATH="$HOME/.pixi/bin:$PATH"

# --- 3. manifesto ---------------------------------------------------------------
if [ -f "$MANIFESTO" ]; then
    mkdir -p "$HOME/.pixi/manifests"
    cp "$MANIFESTO" "$HOME/.pixi/manifests/pixi-global.toml"
    pixi global sync
else
    echo "manifesto $MANIFESTO ausente; so o pixi foi instalado"
    exit 0
fi

# --- 4. snapshot para Metodos --------------------------------------------------
# conda-meta/*.json de cada ambiente tem nome, versao, build, canal e sha256
# exatos do que foi instalado. Lista so o pacote principal de cada ambiente.
python3 - "$HOME/.pixi" "$MANIFESTO" "$SNAPSHOT" <<'PY'
import json, sys, tomllib
from pathlib import Path

pixi_home, manifesto, saida = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
cfg = tomllib.loads(manifesto.read_text())
linhas = ["ambiente\tpacote\tversao\tbuild\tcanal\tsha256"]
for env, spec in sorted(cfg.get("envs", {}).items()):
    meta = pixi_home / "envs" / env / "conda-meta"
    for pacote in sorted(spec.get("dependencies", {})):
        achados = sorted(meta.glob(f"{pacote}-*.json"))
        rec = next((json.loads(p.read_text()) for p in achados if json.loads(p.read_text()).get("name") == pacote), None)
        if rec is None:
            linhas.append(f"{env}\t{pacote}\tAUSENTE\t\t\t")
            continue
        canal = str(rec.get("channel", "")).rstrip("/").rsplit("/", 1)[-1]
        linhas.append(f"{env}\t{pacote}\t{rec['version']}\t{rec.get('build', '')}\t{canal}\t{rec.get('sha256', '')}")
saida.write_text("\n".join(linhas) + "\n")
print(f"snapshot: {saida} ({len(linhas) - 1} pacotes)")
PY

# --- 5. base do TaxonKit -------------------------------------------------------
# O taxonkit nao responde nada sem o taxdump do NCBI. A taxonomia muda toda
# semana, entao a base nao se atualiza sozinha: so baixa se faltar, e a data
# do download fica em ~/.taxonkit/BAIXADO_EM — e ela que o Metodos cita. Para
# atualizar de proposito, apague ~/.taxonkit/names.dmp e rode de novo.
TAXDUMP="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
if [ -s "$HOME/.taxonkit/names.dmp" ]; then
    echo "base do TaxonKit ja presente (baixada em $(cat "$HOME/.taxonkit/BAIXADO_EM" 2>/dev/null || echo '?'))"
else
    tmp_tax="$(mktemp -d)"
    curl -fsSL --retry 3 -o "$tmp_tax/taxdump.tar.gz" "$TAXDUMP"
    curl -fsSL --retry 3 -o "$tmp_tax/taxdump.tar.gz.md5" "$TAXDUMP.md5"
    (cd "$tmp_tax" && md5sum -c taxdump.tar.gz.md5)
    mkdir -p "$HOME/.taxonkit"
    tar -xzf "$tmp_tax/taxdump.tar.gz" -C "$HOME/.taxonkit" names.dmp nodes.dmp delnodes.dmp merged.dmp
    date -u +%Y-%m-%dT%H:%M:%SZ > "$HOME/.taxonkit/BAIXADO_EM"
    rm -rf "$tmp_tax"
    echo "base do TaxonKit baixada em $(cat "$HOME/.taxonkit/BAIXADO_EM")"
fi
