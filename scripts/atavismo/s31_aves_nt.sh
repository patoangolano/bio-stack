#!/usr/bin/env bash
# Nivel 4 (lesao compartilhada), passo 1: achar a reliquia em TODAS as aves por
# nucleotideo, usando a reliquia de uma ave como consulta.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s31_aves_nt.sh
#
# Por que: a calibracao com o Psi-AMEL (s27_calibra_amelx.sh) mostrou que a
# proteina do jacare nao separa reliquia de ave do ruido, enquanto uma consulta
# de ave em nucleotideo acha na hora (100% de identidade no Psi-AMEL). E o s29
# mostrou que o miniprot so alinha a reliquia de ENAM no kakapo e a de AMBN na
# ema. Entao essas duas viram consulta para as outras aves.
#
#   ENAM  consulta = reliquia do kakapo  (Shab, s29: 581874-584251 do intervalo)
#   AMBN  consulta = reliquia da ema     (Dnov, s29: 1215533-1226473)
#   coordenadas lidas de runs/atavismo/lesoes/s29_lesoes.tsv, nao digitadas aqui
#
# Saida: runs/atavismo/nivel4/<GENE>_<rot>.tsv (blastn) e o trecho homologo de
# cada ave, orientado na fita da consulta, em <GENE>_trechos.fna — a entrada do
# alinhamento multiplo do s31_nivel4.py.
set -euo pipefail

A=/mnt/c/Users/Manec/bio-stack/runs/atavismo
L="$A/lncrna"
N="$A/nivel4"
mkdir -p "$N"
AVES=(Ggal Shab Apla Tgut Dnov Scam Arow Tgua)
# referencias com dente para polarizar a lesao (gap so nas aves = lesao de ave)
REF=(Amis Ggan Cpor)
FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'

coord() {  # gene especie -> "ini fim" do melhor alinhamento do miniprot (1-based)
    awk -F'\t' -v g="$1" -v e="$2" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
        $c["gene"]==g && $c["especie"]==e && $c["alinhou"]=="True" {print $c["s_ini"], $c["s_fim"]}' \
        "$A/lesoes/s29_lesoes.tsv"
}

extrai() {  # arquivo ini fim nome -> fasta com o nome dado
    local id
    id=$(head -1 "$1" | cut -c2- | awk '{print $1}')
    [ -s "$1.fai" ] || samtools faidx "$1"
    samtools faidx "$1" "${id}:${2}-${3}" | sed "1s/.*/>$4/"
}

for par in "ENAM Shab" "AMBN Dnov"; do
    read -r gene fonte <<< "$par"
    read -r ini fim <<< "$(coord "$gene" "$fonte")"
    q="$N/${gene}_consulta_${fonte}.fna"
    extrai "$L/intervalo_${fonte}.fna" "$ini" "$fim" "${gene}_${fonte}" > "$q"
    echo "== $gene: consulta = reliquia de $fonte, $((fim - ini + 1)) pb =="
    : > "$N/${gene}_trechos.fna"
    for rot in "${AVES[@]}" "${REF[@]}"; do
        f="$L/intervalo_${rot}.fna"
        [ -s "$f" ] || continue
        tsv="$N/${gene}_${rot}.tsv"
        blastn -task dc-megablast -query "$q" -subject "$f" -evalue 1e-5 -outfmt "$FMT" > "$tsv"
        n=$(wc -l < "$tsv")
        cob=$(awk -F'\t' '{for(i=($5<$6?$5:$6);i<=($5>$6?$5:$6);i++) c[i]=1} END{print length(c)}' "$tsv")
        id=$(awk -F'\t' '{s+=$3*$4; l+=$4} END{if(l) printf "%.0f", s/l; else print "-"}' "$tsv")
        # trecho homologo: do menor ao maior sstart/send dos hits na fita dominante
        if [ "$n" -gt 0 ]; then
            read -r sl sh fita <<< "$(awk -F'\t' '{if($7<$8){p+=$10; lo=(lo==""||$7<lo)?$7:lo; hi=($8>hi)?$8:hi}
                else {m+=$10; lo2=(lo2==""||$8<lo2)?$8:lo2; hi2=($7>hi2)?$7:hi2}}
                END{if(p>=m) print lo, hi, "+"; else print lo2, hi2, "-"}' "$tsv")"
            if [ "$fita" = "+" ]; then
                extrai "$f" "$sl" "$sh" "${gene}_${rot}" >> "$N/${gene}_trechos.fna"
            else
                extrai "$f" "$sl" "$sh" "${gene}_${rot}" | seqkit seq -r -p -t dna 2>/dev/null >> "$N/${gene}_trechos.fna"
            fi
            printf '   %-5s %3s HSPs  consulta coberta %5s/%-5s pb  identidade %3s%%  trecho %s-%s (%s)\n' \
                "$rot" "$n" "$cob" "$((fim - ini + 1))" "$id" "$sl" "$sh" "$fita"
        else
            printf '   %-5s   0 HSPs\n' "$rot"
        fi
    done
done

# --- passo 2: alinhamento multiplo do ENAM na fita codificante --------------------
# Referencias com dente (gene funcional) dao o quadro de leitura e a polaridade:
# o bloco do miniprot no jacare cobre os aa 1-849 contiguos (cg 849M), entao a
# posicao r do jacare cai no codon r/3. Fita de cada um = coluna 5 do melhor PAF.
fita_paf() {  # gene especie -> fita do alinhamento de maior AS
    awk -F'\t' '{for(i=13;i<=NF;i++) if($i~/^AS:i:/) as=substr($i,6); if(as>m){m=as; f=$5}} END{print f}' \
        "$A/lesoes/mp_${1}_${2}.paf"
}
G=ENAM
fs_aves=$(fita_paf $G Shab)
: > "$N/${G}_msa_in.fna"
for rot in "${REF[@]}"; do
    read -r ini fim <<< "$(coord $G "$rot")"
    if [ "$(fita_paf $G "$rot")" = "-" ]; then
        extrai "$L/intervalo_${rot}.fna" "$ini" "$fim" "${G}_${rot}" | seqkit seq -r -p -t dna 2>/dev/null >> "$N/${G}_msa_in.fna"
    else
        extrai "$L/intervalo_${rot}.fna" "$ini" "$fim" "${G}_${rot}" >> "$N/${G}_msa_in.fna"
    fi
done
if [ "$fs_aves" = "-" ]; then
    seqkit seq -r -p -t dna "$N/${G}_trechos.fna" 2>/dev/null >> "$N/${G}_msa_in.fna"
else
    cat "$N/${G}_trechos.fna" >> "$N/${G}_msa_in.fna"
fi
echo
echo "== $G: alinhamento multiplo (fita codificante; aves orientadas pela fita do kakapo: $fs_aves) =="
seqkit fx2tab -n -l "$N/${G}_msa_in.fna" 2>/dev/null | awk '{printf "   %-12s %6s pb\n", $1, $2}'
mafft --localpair --maxiterate 1000 --thread 8 --quiet "$N/${G}_msa_in.fna" > "$N/${G}_mafft.fna"
muscle -align "$N/${G}_msa_in.fna" -output "$N/${G}_muscle.fna" -threads 8 >/dev/null 2>&1
echo "   mafft L-INS-i: $(seqkit fx2tab -n -l "$N/${G}_mafft.fna" 2>/dev/null | awk 'NR==1{print $2}') colunas"
echo "   muscle 5:      $(seqkit fx2tab -n -l "$N/${G}_muscle.fna" 2>/dev/null | awk 'NR==1{print $2}') colunas  ($(muscle -version 2>&1 | head -1))"
