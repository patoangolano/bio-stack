#!/usr/bin/env bash
# bio-stack — restaura o backup do Google Drive: runs/ e data/ do projeto no
# Windows, e ~/rapinantes e ~/nf-work no WSL.
#
# Roda no Git Bash do Windows e pega os arquivos no G: do Google Drive para
# desktop, sem rclone. Os tres .tar.gz estao la, inclusive o nf-work, que subiu
# pelo G: em 2026-10-19 e que o rclone (escopo drive.file) nem enxerga. Cada um
# e copiado para runs/bootstrap/staging/ antes de ser lido (ver restaurar()).
#
#   bash bootstrap/restaura_drive.sh                     # log em runs/bootstrap/
#   .\bootstrap\tarefa_longa.ps1 -Nome restaura-drive -Script bootstrap\restaura_drive.sh
#
# Tres garantias, nesta ordem:
#   1. md5 de cada .tar.gz conferido com o LEIAME.txt do backup ANTES de
#      extrair. Arquivo que nao bate nao e extraido.
#   2. tar --skip-old-files: arquivo que ja existe nunca e sobrescrito. E o que
#      protege o runs/provenance.jsonl e tudo que o git trouxe mais novo que o
#      backup. Um tar sem a flag trocaria o caderno de laboratorio de hoje pelo
#      de 19/10.
#   3. Item ja restaurado com o mesmo md5 fica marcado em runs/bootstrap/ e e
#      pulado na proxima vez — reler 8 GB do Drive para nada leva tempo.
#
# Os itens do WSL so entram se a distro e o usuario manec existirem (depois do
# 04b); antes disso sao pulados, e uma segunda rodada os pega.

set -u

DISTRO=Ubuntu-24.04
ORIGEM="${BIO_STACK_BACKUP:-/g/Meu Drive/bio-stack-backup/2026-10-19}"
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(dirname "$AQUI")"
MARCAS="$RAIZ/runs/bootstrap"
STAGING="$MARCAS/staging"
LOG="${1:-$MARCAS/restaura_drive.log}"
mkdir -p "$MARCAS" "$(dirname "$LOG")"

registrar() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }

md5_esperado() {
    awk -v a="$1" '$1 == a && length($2) == 32 && $2 ~ /^[0-9a-f]+$/ { print $2 }' "$ORIGEM/LEIAME.txt"
}

wsl_pronto() {
    wsl.exe -l -q 2>/dev/null | tr -d '\0\r' | grep -qx "$DISTRO" &&
        MSYS_NO_PATHCONV=1 wsl.exe -d "$DISTRO" -u manec -- test -d /home/manec 2>/dev/null
}

# nome | destino (windows: diretorio do Git Bash; wsl: diretorio na distro)
restaurar() {
    local nome="$1" lado="$2" destino="$3"
    local arq="$ORIGEM/$nome.tar.gz" marca="$MARCAS/restaurado_$nome.md5"
    local esperado medido
    registrar "=== $nome -> $lado:$destino ==="
    esperado=$(md5_esperado "$nome.tar.gz")
    if [ -z "$esperado" ]; then
        registrar "FALHA $nome: sem md5 no LEIAME.txt"
        return 1
    fi
    if [ -f "$marca" ] && [ "$(cat "$marca")" = "$esperado" ]; then
        registrar "OK   $nome ja restaurado (md5 $esperado)"
        return 0
    fi
    if [ ! -f "$arq" ]; then
        registrar "FALHA $nome: $arq nao existe (G: montado?)"
        return 1
    fi
    # Copia para disco local pelo leitor nativo do Windows, e so depois le.
    # Ler o G: direto pelo Git Bash falhou no meio do arquivo com "Invalid
    # request code" (EBADRQC do sistema de arquivos virtual do Drive para
    # desktop, medido em 2026-10-23); o Copy-Item do .NET atravessa inteiro.
    local copia="$STAGING/$nome.tar.gz"
    mkdir -p "$STAGING"
    registrar "  copia local ($(du -h "$arq" | cut -f1)) -> $copia"
    if ! powershell.exe -NoProfile -Command \
        "Copy-Item -LiteralPath '$(cygpath -w "$arq")' -Destination '$(cygpath -w "$copia")' -Force" \
        </dev/null >>"$LOG" 2>&1; then
        registrar "FALHA $nome: a copia do G: nao completou"
        return 1
    fi
    # O arquivo entra pelo stdin, e nao como argumento: sob o conhost
    # --headless da tarefa agendada o stdin herdado e invalido, e o md5sum do
    # Git Bash morre tentando poe-lo em modo binario ("failed to set file
    # descriptor text/binary mode", medido em 2026-10-23) — com md5 vazio.
    medido=$(md5sum <"$copia" | cut -d' ' -f1)
    if [ "$medido" != "$esperado" ]; then
        registrar "FALHA $nome: md5 $medido, o LEIAME diz $esperado — nao extraido"
        registrar "  baixe o .tar.gz pelo navegador (drive.google.com) para uma pasta, copie o"
        registrar "  LEIAME.txt para ela e rode de novo com BIO_STACK_BACKUP=<essa pasta>"
        rm -f "$copia"
        return 1
    fi
    registrar "  md5 confere: $medido"
    if [ "$lado" = windows ]; then
        tar -xzf "$copia" --skip-old-files -C "$destino" </dev/null >>"$LOG" 2>&1
    else
        # Fluxo pelo stdin do wsl.exe: os bytes passam crus, o dono dos
        # arquivos e o manec, e nao e preciso montar o G: dentro da distro
        # (o que pediria sudo, e o manec nao tem senha).
        MSYS_NO_PATHCONV=1 wsl.exe -d "$DISTRO" -u manec -- \
            tar -xz --skip-old-files -C "$destino" <"$copia" >>"$LOG" 2>&1
    fi
    local codigo=$?
    if [ "$codigo" -eq 0 ]; then
        echo "$esperado" >"$marca"
        rm -f "$copia"
        registrar "OK   $nome"
    else
        registrar "FALHA $nome: tar saiu com $codigo (ver acima no log; copia mantida em $copia)"
    fi
}

registrar "inicio — origem $ORIGEM"
if [ ! -f "$ORIGEM/LEIAME.txt" ]; then
    registrar "FALHA $ORIGEM/LEIAME.txt nao encontrado — o Google Drive para desktop esta aberto?"
    exit 1
fi

restaurar bio-stack_runs_data windows "$RAIZ"

if wsl_pronto; then
    restaurar rapinantes_wsl wsl /home/manec
    restaurar nf-work_wsl wsl /home/manec
else
    registrar "pula rapinantes_wsl e nf-work_wsl: $DISTRO sem o usuario manec (rode o 09 antes)"
fi

# ~/genomas (9,7 GB) e ~/nf-smoke nao foram para o backup: tudo publico, refeito
# por scripts/atavismo/s35_genomas.sh, s36_multiz.sh e pelo nf-core/demo.
registrar "fim"
