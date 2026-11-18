#!/usr/bin/env bash
# Verifica acessos de genoma citados na pesquisa de microrganismos energéticos.
# Saída: runs/microbios_energia/genomas_ncbi.tsv
set -euo pipefail
cd /mnt/c/Users/user/bio-stack
mkdir -p runs/microbios_energia
OUT=runs/microbios_energia/genomas_ncbi.jsonl
: > "$OUT"
ACCS="GCF_000195775.1 GCF_000009285.2 GCF_000007985.2 GCF_000146165.2 GCF_000007345.1 GCF_000247605.1 GCF_000143685.1 GCF_000013105.1 GCF_000021485.1 GCF_011066545.1 GCF_000018425.1 GCF_004347865.1 GCF_000018465.1 GCF_000008625.1 GCF_000012865.1 GCF_000018865.1 GCF_000017945.1 GCF_000009725.1"
for a in $ACCS; do
  r=$(datasets summary genome accession "$a" --as-json-lines 2>/dev/null || true)
  if [ -z "$r" ]; then
    echo "{\"query\":\"$a\",\"status\":\"NAO_ENCONTRADO\"}" >> "$OUT"
  else
    echo "$r" | jq -c --arg q "$a" '{query:$q, status:"OK", accession, organism:.organism.organism_name, strain:(.organism.infraspecific_names.strain // ""), level:.assembly_info.assembly_level, refseq_cat:(.assembly_info.refseq_category // ""), size_bp:(.assembly_stats.total_sequence_length|tonumber), gc:(.assembly_stats.gc_percent // null), n_seq:(.assembly_stats.number_of_component_sequences // null), cds:(.annotation_info.stats.gene_counts.protein_coding // null), status_asm:(.assembly_info.assembly_status // "")}' >> "$OUT"
  fi
done
jq -r '[.query,.status,(.accession//""),(.organism//""),(.strain//""),(.level//""),((.size_bp//0)/1e6|tostring),((.gc//"")|tostring),((.cds//"")|tostring),(.status_asm//"")]|@tsv' "$OUT" > runs/microbios_energia/genomas_ncbi.tsv
cat runs/microbios_energia/genomas_ncbi.tsv
