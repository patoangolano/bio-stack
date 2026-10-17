#!/usr/bin/env bash
# s32 — ODAM e ODAPH. Alinhamento no WSL; download e leitura no s32_odam_odaph.py.
#
#   export MSYS_NO_PATHCONV=1
#   wsl.exe -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s32_odam_odaph.sh etapa1
#   wsl.exe -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s32_odam_odaph.sh etapa2
#
# etapa1  so o jacare: consultas de mamifero/anfibio contra o intervalo SCPP
#         (ODAM) e a janela CDKL2..RCHY1 (ODAPH); blastp das consultas contra o
#         proteoma do jacare (primeiro lado do reciproco).
# etapa2  modelo do jacare (escrito pelo `s32_odam_odaph.py modelo`) contra o
#         Swiss-Prot humano (segundo lado do reciproco) e contra todas as
#         especies, junto com os controles de ancora no MESMO comando.
#
# tblastn: evalue 10, saida completa; o limiar de decisao e aplicado na leitura
# (s32_odam_odaph.py tabela), igual para alvo e controle. Duas passadas, SEG
# ligado (comparavel ao s26/s27) e desligado (ODAPH e curto; o motivo SSEE e
# baixa complexidade e o SEG pode mascara-lo).
# miniprot: as opcoes de reliquia do s29 (-k 4 -M 0 -n 1 -l 4), que o s29 mediu
# nao criarem lesao falsa no ENAM do Anolis.
set -euo pipefail

P=/mnt/c/Users/Manec/bio-stack/runs/atavismo
O="$P/odam_odaph"
DB="$O/db"
A="$O/aln"
mkdir -p "$DB" "$A"
FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'
MP=(-j 2 -k 4 -M 0 -n 1 -l 4 --outc=0.01 --outs=0.5 -t 8)
echo "tblastn $(tblastn -version | head -1 | awk '{print $2}')  miniprot $(miniprot --version)"

db_nucl() {  # $1=fasta  $2=nome
    [ -s "$DB/$2.nsq" ] || [ -s "$DB/$2.00.nsq" ] || makeblastdb -in "$1" -dbtype nucl -out "$DB/$2" >/dev/null
}

busca() {  # $1=faa  $2=fasta alvo  $3=rotulo de saida
    local nome
    nome=$(basename "$2" .fna)
    db_nucl "$2" "$nome"
    tblastn -query "$1" -db "$DB/$nome" -evalue 10 -seg yes -max_target_seqs 50 -num_threads 8 \
        -outfmt "$FMT" > "$A/tbn_seg_$3.tsv"
    tblastn -query "$1" -db "$DB/$nome" -evalue 10 -seg no -max_target_seqs 50 -num_threads 8 \
        -outfmt "$FMT" > "$A/tbn_noseg_$3.tsv"
    miniprot "${MP[@]}" --trans "$2" "$1" 2>/dev/null > "$A/mp_$3.paf" || true
    printf '  %-28s tblastn %4s/%4s HSPs (seg/noseg)  miniprot %3s\n' "$3" \
        "$(wc -l < "$A/tbn_seg_$3.tsv")" "$(wc -l < "$A/tbn_noseg_$3.tsv")" \
        "$(grep -vc '^##' "$A/mp_$3.paf" || true)"
}

etapa1() {
    # consultas de fora + controles de ancora num arquivo so: mesmo comando, mesmo limiar
    cat "$O/consultas_odam.faa" "$O/controles_odam.faa" > "$A/e1_q_odam.faa"
    cat "$O/consultas_odaph.faa" "$O/controles_odaph.faa" > "$A/e1_q_odaph.faa"
    echo "== ODAM: consultas de fora + ENAM/AMBN/AMTN contra o intervalo SCPP de cada especie"
    for f in "$P"/lncrna/intervalo_*.fna; do
        rot=$(basename "$f" .fna); rot=${rot#intervalo_}
        busca "$A/e1_q_odam.faa" "$f" "e1_odam_$rot"
    done
    echo "== ODAM: controle positivo de sensibilidade (regiao do ODAM de Xenopus, 1,3 Mb)"
    busca "$O/consultas_odam.faa" "$O/controle_odam_Xtro.fna" "e1_odam_ctrlXtro"
    echo "== ODAM: perfil Pfam PF15424 contra os ORFs de 6 quadros (>= 15 aa) de cada intervalo"
    for f in "$P"/lncrna/intervalo_*.fna "$O/controle_odam_Xtro.fna"; do
        rot=$(basename "$f" .fna); rot=${rot#intervalo_}; rot=${rot#controle_odam_}
        esl-translate -l 15 "$f" > "$A/orf_$rot.faa"
        # -Z fixo: E-value comparavel entre especies, independente do tamanho do intervalo
        hmmsearch --max -E 10 --domE 10 -Z 50000 --domtblout "$A/hmm_odam_$rot.dom" \
            "$O/pf15424_odam.hmm" "$A/orf_$rot.faa" > /dev/null
        printf '  %-6s %6s ORFs  %3s dominios (E<=10)\n' "$rot" "$(grep -c '^>' "$A/orf_$rot.faa")" \
            "$(grep -vc '^#' "$A/hmm_odam_$rot.dom" || true)"
    done
    echo "== ODAPH: consultas de fora + CDKL2/THAP6/RCHY1/G3BP2 contra a janela de cada especie"
    for f in "$O"/janela_odaph_*.fna; do
        rot=$(basename "$f" .fna); rot=${rot#janela_odaph_}
        busca "$A/e1_q_odaph.faa" "$f" "e1_odaph_$rot"
    done
    echo "== blastp: consultas contra o proteoma do jacare (reciproco, lado 1)"
    cat "$O/consultas_odam.faa" "$O/consultas_odaph.faa" > "$A/e1_consultas.faa"
    blastp -query "$A/e1_consultas.faa" -db "$P/painel/proteoma_Amis/db" -evalue 10 -seg no \
        -max_target_seqs 10 -num_threads 8 -outfmt "$FMT" > "$A/e1_bp_consultas_vs_Amis.tsv"
    [ -s "$DB/sprot_humano.psq" ] || makeblastdb -in "$O/sprot_humano.faa" -dbtype prot -out "$DB/sprot_humano" >/dev/null
    echo "  ok: $(wc -l < "$A/e1_bp_consultas_vs_Amis.tsv") HSPs"
}

etapa2() {
    for m in odam odaph; do
        [ -s "$O/modelo_${m}_Amis.faa" ] || { echo "  modelo_${m}_Amis.faa ausente (rode o 'modelo')"; continue; }
        echo "== $m: modelo do jacare contra Swiss-Prot humano e proteoma do jacare (reciproco, lado 2)"
        blastp -query "$O/modelo_${m}_Amis.faa" -db "$DB/sprot_humano" -evalue 10 -seg no \
            -max_target_seqs 10 -num_threads 8 -outfmt "$FMT" > "$A/e2_bp_${m}_vs_sprot.tsv"
        blastp -query "$O/modelo_${m}_Amis.faa" -db "$P/painel/proteoma_Amis/db" -evalue 10 -seg no \
            -max_target_seqs 10 -num_threads 8 -outfmt "$FMT" > "$A/e2_bp_${m}_vs_Amis.tsv"
        # consulta + controles num arquivo so: mesmo comando, mesmo limiar
        cat "$O/modelo_${m}_Amis.faa" "$O/controles_${m}.faa" > "$A/e2_q_${m}.faa"
    done
    echo "== ODAM: modelo do jacare + ENAM/AMBN/AMTN contra o intervalo SCPP de cada especie"
    for f in "$P"/lncrna/intervalo_*.fna; do
        rot=$(basename "$f" .fna); rot=${rot#intervalo_}
        [ -s "$A/e2_q_odam.faa" ] && busca "$A/e2_q_odam.faa" "$f" "e2_odam_$rot"
    done
    echo "== ODAPH: modelo do jacare + CDKL2/THAP6/RCHY1/G3BP2 contra a janela de cada especie"
    for f in "$O"/janela_odaph_*.fna; do
        rot=$(basename "$f" .fna); rot=${rot#janela_odaph_}
        [ -s "$A/e2_q_odaph.faa" ] && busca "$A/e2_q_odaph.faa" "$f" "e2_odaph_$rot"
    done
    [ -s "$O/modelo_odaph_Amis.faa" ] || return 0
    echo "== ODAPH: so o vao CDKL2..vizinho de fora (proteina e gene inteiro do jacare)"
    for f in "$O"/vao_odaph_*.fna; do
        rot=$(basename "$f" .fna); rot=${rot#vao_odaph_}
        busca "$O/modelo_odaph_Amis.faa" "$f" "e2_vao_$rot"
        # nucleotideo: exons + introns do jacare contra o vao (dc-megablast, sensivel a ~250 Ma)
        blastn -task dc-megablast -query "$O/gene_odaph_Amis.fna" -subject "$f" -evalue 1e-3 \
            -outfmt "$FMT" > "$A/bn_vao_$rot.tsv"
        printf '  %-28s blastn  %4s HSPs (E<=1e-3)\n' "e2_vao_$rot" "$(wc -l < "$A/bn_vao_$rot.tsv")"
    done
    if [ -s "$O/transcritos_vao.fna" ]; then
        echo "== ODAPH: ORFs dos transcritos LOC (lncRNA) no vao contra o modelo do jacare"
        # cabecalho do ORF vira fonte|orf|coords|tamanho para a leitura no Python
        esl-translate -l 30 "$O/transcritos_vao.fna" | awk '/^>/{split($2,s,"=");split($3,c,"=");split($4,l,"=");print ">" s[2] "|" substr($1,2) "|" c[2] "|" l[2]; next}{print}' > "$A/orf_transcritos_vao.faa"
        blastp -query "$O/modelo_odaph_Amis.faa" -subject "$A/orf_transcritos_vao.faa" -evalue 10 -seg no \
            -outfmt "$FMT" > "$A/bp_odaph_vs_orf_transcritos.tsv"
        echo "  $(grep -c '^>' "$A/orf_transcritos_vao.faa") ORFs >= 30 aa; $(wc -l < "$A/bp_odaph_vs_orf_transcritos.tsv") HSPs"
        # Segunda consulta, do outro lado da arvore: o ORF do Podarcis (escamado,
        # com dente) que melhor alinha ao modelo do jacare, se E <= 1e-3. Jacare
        # e Podarcis cercam tartarugas e aves pelos dois lados de Sauria.
        orf=$(sort -t$'\t' -k9,9g "$A/bp_odaph_vs_orf_transcritos.tsv" | awk -F'\t' '$2 ~ /^Pmur_/ && $9 <= 1e-3 {print $2; exit}')
        if [ -n "$orf" ]; then
            awk -v id=">$orf" 'BEGIN{p=0} /^>/{p=($0==id)} p' "$A/orf_transcritos_vao.faa" \
                | sed "1s/.*/>ODAPH_Pmur $orf/" > "$O/modelo_odaph_Pmur.faa"
            echo "== ODAPH: ORF do Podarcis ($orf) contra os vaos"
            for f in "$O"/vao_odaph_*.fna; do
                rot=$(basename "$f" .fna); rot=${rot#vao_odaph_}
                busca "$O/modelo_odaph_Pmur.faa" "$f" "e2_vaoPmur_$rot"
            done
        fi
    fi
}

case "${1:-}" in
    etapa1) etapa1 ;;
    etapa2) etapa2 ;;
    *) echo "uso: $0 {etapa1|etapa2}" >&2; exit 2 ;;
esac
