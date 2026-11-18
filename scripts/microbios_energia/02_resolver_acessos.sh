#!/usr/bin/env bash
# Resolve acessos ausentes/obsoletos: C. necator H16, Electrothrix aarhusiensis,
# R. palustris CGA009 (versão nova), M. thermoacetica ATCC 39073 (suprimido).
set -uo pipefail
cd /mnt/c/Users/user/bio-stack
OUT=runs/microbios_energia/genomas_resolvidos.tsv
: > "$OUT"
for t in "Cupriavidus necator H16" "Candidatus Electrothrix aarhusiensis" "Rhodopseudomonas palustris CGA009" "Moorella thermoacetica"; do
  datasets summary genome taxon "$t" --as-json-lines --assembly-source all 2>/dev/null \
   | jq -r --arg q "$t" '[$q, .accession, .organism.organism_name, (.organism.infraspecific_names.strain // ""), .assembly_info.assembly_level, (.assembly_info.assembly_status // ""), (.assembly_info.refseq_category // ""), ((.assembly_stats.total_sequence_length|tonumber)/1e6|tostring), ((.assembly_stats.gc_percent // "")|tostring), ((.annotation_info.stats.gene_counts.protein_coding // "")|tostring)]|@tsv' \
   | grep -E "H16|aarhus|CGA009|39073|Complete" | head -n 8 >> "$OUT"
done
cat "$OUT"
