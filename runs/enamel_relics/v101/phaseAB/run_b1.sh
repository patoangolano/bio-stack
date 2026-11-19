#!/usr/bin/env bash
set -euo pipefail
D="/mnt/c/Users/user/bio-stack/runs/enamel_relics"
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/phaseAB"
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe stitle"
# 1. Standard reciprocal BLASTx (same params as Table S2b: -evalue 10 -seg yes -max_target_seqs 25)
blastx -query "$O/b1_ggal_enam_ex6_queries.fna" -db "$D/blastdb/amis_proteome" -evalue 10 -seg yes -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$O/b1_blastx_std_segyes.tsv"
# 2. Standard reciprocal BLASTx with -seg no
blastx -query "$O/b1_ggal_enam_ex6_queries.fna" -db "$D/blastdb/amis_proteome" -evalue 10 -seg no -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$O/b1_blastx_std_segno.tsv"
# 3. Ultra-sensitive reciprocal BLASTx for short exon (-evalue 1000 -word_size 2 -seg no -comp_based_stats 0)
blastx -query "$O/b1_ggal_enam_ex6_queries.fna" -db "$D/blastdb/amis_proteome" -evalue 1000 -word_size 2 -seg no -comp_based_stats 0 -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$O/b1_blastx_sens_w2.tsv"
# 4. Re-run BLASTn of Arow & Amis Exon 6 against Ggal NRG1-LPL interval
blastn -task blastn -word_size 7 -query "$D/level4/ENAM_exon6_splice_aves/queries.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-4 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "$O/b1_blastn_arow_vs_ggal.tsv"
# 5. Direct BLASTn of Amis Exon 6 (+-60bp) against Ggal NRG1-LPL interval
printf ">Amis_ex6_pad168\nCGCACCCCCAAGCAATGGTCTTTCATTAGAATGATGTGTCTATTTTCCTTCTGTCTGCAGATTTACAACCAGGGAGGCTTCTTACCACAGCTGAACCCACAACACAGAGTAAGATGCTTTATTTAAGGGCATTTCTAGTACTGCTAAACTGTCTTTTATTTTATTTTT\n>Amis_ex6_core48\nATTTACAACCAGGGAGGCTTCTTACCACAGCTGAACCCACAACACAGA\n" > "$O/amis_ex6_q.fna"
blastn -task blastn -word_size 7 -query "$O/amis_ex6_q.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "$O/b1_blastn_amis_vs_ggal.tsv"
# 6. Single-copy check of Ggal Exon 6 (149 bp) inside Ggal NRG1-LPL interval (both dc-megablast and blastn)
printf ">Ggal_ENAM_ex6_149\nGTTTGCCCAACAATGATCTTTCATGCATATAATTTGTCCATTTGTGTGTGTGTGTGCACAGATTTACAACTAAGAGGCTTCTGGCTACAGCTGAACATGTTCCAAGGTGTAAGAGGCTTTCTCTTAAAGGCAATTCTATTTGTGATGAA\n" > "$O/ggal_ex6_149.fna"
blastn -task dc-megablast -query "$O/ggal_ex6_149.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-3 -dust yes -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" > "$O/b1_singlecopy_dcmegablast.tsv"
blastn -task blastn -query "$O/ggal_ex6_149.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-3 -dust yes -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" > "$O/b1_singlecopy_blastn.tsv"
