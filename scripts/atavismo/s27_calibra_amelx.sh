#!/usr/bin/env bash
# Calibracao do s27 com um pseudogene conhecido: o Psi-AMEL da galinha.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s27_calibra_amelx.sh
#
# Sire et al. 2008 (BMC Evol Biol 8:246, PMID 18775069) depositaram o Psi-AMEL
# de Gallus no GenBank: EU340348, 513 pb, exons 2/3/5/6, ~38 kb a montante do
# exon 2 de ARHGAP6 (assembly WASHUC2). O s27 montou o intervalo MSL3..MID1 no
# chr1 de GRCg7b — que deveria conte-lo — e o tblastn com o AMELX do jacare deu
# ZERO bits ali.
#
# Duas perguntas:
#   A. O Psi-AMEL esta dentro do intervalo_AMELX_Ggal.fna? (blastn, nucleotideo
#      contra nucleotideo, ave contra ave: sensibilidade alta)
#   B. Se esta, o tblastn do s27 tinha HSP nele, ainda que abaixo do ruido?
# Se A=sim e B=nao, o s27 perde um pseudogene que existe: o zero do AMELX nas
# aves e falta de sensibilidade, nao ausencia. E o que o controle do Anolis ja
# sugeria.
set -euo pipefail
D=/mnt/c/Users/Manec/bio-stack/runs/atavismo/painel
cd "$D"

[ -s psiamel_EU340348.fna ] || efetch -db nuccore -id EU340348 -format fasta > psiamel_EU340348.fna
echo "consulta: $(head -1 psiamel_EU340348.fna | cut -c1-90)"
echo "          $(grep -v '^>' psiamel_EU340348.fna | tr -d '\n' | wc -c) pb"

echo
echo "=== A. blastn do Psi-AMEL contra o intervalo AMELX da galinha (GRCg7b chr1) ==="
cab=$(head -1 intervalo_AMELX_Ggal.fna | cut -c2-)
echo "intervalo: $cab"
blastn -query psiamel_EU340348.fna -subject intervalo_AMELX_Ggal.fna -evalue 1e-5 \
       -outfmt '6 qstart qend sstart send pident length evalue bitscore' > calibra_psiamel.tsv
if [ -s calibra_psiamel.tsv ]; then
    printf '%-14s %-18s %7s %6s %9s %7s\n' QUERY SUJEITO %ID COMPR EVALUE BITS
    awk -F'\t' '{printf "%-14s %-18s %7s %6s %9s %7s\n", $1"-"$2, $3"-"$4, $5, $6, $7, $8}' calibra_psiamel.tsv
    ini=$(awk -F'\t' 'NR==1{print $3}' calibra_psiamel.tsv)
    # coordenada absoluta em GRCg7b = inicio da janela + posicao - 1
    janela=$(echo "$cab" | grep -oE ':[0-9]+-' | tr -d ':-')
    echo "coordenada aproximada em GRCg7b (NC_052572 nao; e o chr1): $(( janela + ini - 1 ))"
else
    echo "nenhum hit: o Psi-AMEL NAO esta no intervalo MSL3..MID1"
fi

echo
echo "=== B. HSPs do tblastn do s27 (AMELX do jacare) na galinha, qualquer bitscore ==="
if [ -s calibra_psiamel.tsv ]; then
    lo=$(awk -F'\t' '{m=($3<$4?$3:$4); if (NR==1||m<a) a=m} END{print a-2000}' calibra_psiamel.tsv)
    hi=$(awk -F'\t' '{m=($3>$4?$3:$4); if (m>b) b=m} END{print b+2000}' calibra_psiamel.tsv)
    echo "janela do Psi-AMEL (+-2 kb): $lo-$hi"
    awk -F'\t' -v lo="$lo" -v hi="$hi" '$1=="AMELX_Amis" {s=($7<$8?$7:$8); if (s>=lo && s<=hi) print}' \
        s27_AMELX_Ggal.tsv | sort -t$'\t' -k10,10gr | head -5
    n=$(awk -F'\t' -v lo="$lo" -v hi="$hi" '$1=="AMELX_Amis" {s=($7<$8?$7:$8); if (s>=lo && s<=hi) c++} END{print c+0}' s27_AMELX_Ggal.tsv)
    echo "HSPs do AMELX do jacare sobre o Psi-AMEL: $n"
fi

echo
echo "=== C. O mesmo tblastn do s27, com SEG DESLIGADO, so o AMELX do jacare ==="
# Hipotese: o AMELX e rico em Pro/Gln de baixa complexidade e o SEG mascara a
# consulta quase inteira. Se com -seg no o Psi-AMEL aparecer, a perda do s27 e
# do filtro, nao da distancia.
awk '/^>/{p=($1==">AMELX_Amis")} p' painel_AMELX_proteinas.faa > amelx_jacare.faa
tblastn -query amelx_jacare.faa -subject intervalo_AMELX_Ggal.fna -seg no -evalue 10 \
        -outfmt '6 qstart qend sstart send pident length evalue bitscore' > calibra_amelx_segno.tsv
if [ -s calibra_amelx_segno.tsv ]; then
    lo=$(awk -F'\t' '{m=($3<$4?$3:$4); if (NR==1||m<a) a=m} END{print a-2000}' calibra_psiamel.tsv)
    hi=$(awk -F'\t' '{m=($3>$4?$3:$4); if (m>b) b=m} END{print b+2000}' calibra_psiamel.tsv)
    printf '%-10s %-18s %6s %6s %9s %7s %s\n' QUERY SUJEITO %ID COMPR EVALUE BITS "SOBRE O Psi-AMEL?"
    sort -t$'\t' -k8,8gr calibra_amelx_segno.tsv | head -8 |
    awk -F'\t' -v lo="$lo" -v hi="$hi" '{s=($3<$4?$3:$4); d=(s>=lo && s<=hi)?"SIM":"nao";
        printf "%-10s %-18s %6s %6s %9s %7s %s\n", $1"-"$2, $3"-"$4, $5, $6, $7, $8, d}'
else
    echo "nenhum HSP nem com SEG desligado"
fi
