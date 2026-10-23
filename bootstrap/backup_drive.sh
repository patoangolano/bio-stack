#!/usr/bin/env bash
# bio-stack — manda para o Google Drive o que o git nao guarda: runs/ e data/
# do projeto, e os diretorios de trabalho do WSL. Sem escrever nada local: cada
# item vai como tar | pigz | rclone rcat, com o md5 do fluxo calculado no
# caminho e comparado com o md5 que o Drive devolve.
#
# Roda DENTRO do WSL, com o remoto gdrive ja criado (bootstrap/rclone_gdrive.sh).
# Upload de horas, entao vai por tarefa agendada:
#
#   .\bootstrap\tarefa_longa.ps1 -Nome backup-drive -Script bootstrap\backup_drive.sh -Wsl
#
# Restaurado do backup_drive.sh e do backup_nfwork_retry.sh, que viviam no
# scratchpad da sessao de 2026-10-19. Aquele backup salvou o projeto: em
# 2026-10-23 a distro foi desregistrada e o repositorio re-clonado, e o que
# voltou, voltou dali (bootstrap/restaura_drive.sh).
#
# O que 2026-10-19 ensinou, ja aplicado:
#   - a cota do client_id compartilhado do rclone estourou (403
#     rateLimitExceeded) no terceiro item. Cada item agora tem ate 4 tentativas,
#     com --tpslimit 4 e --low-level-retries 40, e 90 s entre elas;
#   - tar sai com 1 quando um arquivo muda durante a leitura (o
#     provenance.jsonl cresce enquanto a sessao roda): e aviso, nao falha. O
#     criterio e o md5 do Drive bater com o do fluxo e o rclone sair com 0;
#   - o LEIAME leva o md5 de cada .tar.gz, no formato que o restaura_drive.sh
#     confere antes de extrair.
#
# Itens: os de 2026-10-19 por padrao. Para outros, BIO_STACK_BACKUP_ITENS com
# "nome|diretorio base|membros" separados por ';'. Destino:
# gdrive:bio-stack-backup/<data de hoje>, ou BIO_STACK_BACKUP_DEST.

set -uo pipefail

DEST="${BIO_STACK_BACKUP_DEST:-gdrive:bio-stack-backup/$(date +%F)}"
LOG="${1:-$HOME/backup_drive.log}"
TMP="$HOME/.backup_drive_md5"
BIO=/mnt/c/Users/Manec/bio-stack
mkdir -p "$TMP" "$(dirname "$LOG")"

if [ -n "${BIO_STACK_BACKUP_ITENS:-}" ]; then
    IFS=';' read -r -a ITENS <<< "$BIO_STACK_BACKUP_ITENS"
else
    ITENS=(
        "rapinantes_wsl|$HOME|rapinantes"
        "nf-work_wsl|$HOME|nf-work"
        "bio-stack_runs_data|$BIO|runs data"
    )
fi

log() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

for exe in rclone pigz; do
    command -v "$exe" >/dev/null || { log "FALHA: $exe ausente — rode bootstrap/05_wsl_ferramentas.sh"; exit 1; }
done
rclone listremotes 2>/dev/null | grep -qx 'gdrive:' ||
    { log "FALHA: remoto gdrive ausente — rode bootstrap/rclone_gdrive.sh"; exit 1; }

log "inicio -> $DEST"
resumo=()

for item in "${ITENS[@]}"; do
    IFS='|' read -r nome base membros <<< "$item"
    arq="$nome.tar.gz"
    md5f="$TMP/$nome.md5"
    log "$nome: lista de arquivos"
    # shellcheck disable=SC2086
    (cd "$base" && find $membros -type f -printf '%s\t%TY-%Tm-%Td %TH:%TM\t%p\n') |
        rclone rcat "$DEST/$nome.arquivos.tsv" 2>>"$LOG"
    ok=FALHA
    for tentativa in 1 2 3 4; do
        rm -f "$md5f"
        log "$nome: envio de $arq, tentativa $tentativa"
        # shellcheck disable=SC2086
        tar -C "$base" -cf - $membros 2>>"$LOG" | pigz -p 16 -6 |
            tee >(md5sum | cut -d' ' -f1 >"$md5f") |
            rclone rcat "$DEST/$arq" --drive-chunk-size 64M --tpslimit 4 --low-level-retries 40 2>>"$LOG"
        st=("${PIPESTATUS[@]}")
        until [ -s "$md5f" ]; do sleep 1; done
        local_md5=$(cat "$md5f")
        remoto_md5=$(rclone md5sum "$DEST/$arq" 2>>"$LOG" | cut -d' ' -f1)
        tam=$(rclone size --json "$DEST/$arq" 2>>"$LOG" | grep -oE '"bytes":[0-9]+' | cut -d: -f2)
        if [ "$local_md5" = "$remoto_md5" ] && [ "${st[3]}" = 0 ]; then ok=ok; fi
        log "$nome: $ok tar=${st[0]} pigz=${st[1]} rclone=${st[3]} bytes=$tam md5_local=$local_md5 md5_drive=$remoto_md5"
        [ "$ok" = ok ] && break
        sleep 90
    done
    resumo+=("$(printf '  %-28s%s' "$arq" "$([ "$ok" = ok ] && echo "$local_md5" || echo "FALHOU")")")
done

{
    echo "Backup do bio-stack feito em $(date +%F) pelo bootstrap/backup_drive.sh."
    echo
    echo "Cada <nome>.tar.gz tem ao lado <nome>.arquivos.tsv: tamanho em bytes,"
    echo "data de modificacao e caminho de cada arquivo."
    echo
    echo "md5 de cada .tar.gz (o mesmo do fluxo gerado no WSL):"
    printf '%s\n' "${resumo[@]}"
    echo
    echo "Restaurar: bash bootstrap/restaura_drive.sh, com BIO_STACK_BACKUP"
    echo "apontando para esta pasta no G: (Google Drive para desktop)."
} | rclone rcat "$DEST/LEIAME.txt" 2>>"$LOG"

log "fim"
