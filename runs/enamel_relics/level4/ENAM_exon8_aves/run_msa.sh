#!/usr/bin/env bash
set -euo pipefail
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/level4/ENAM_exon8_aves"
mafft --localpair --maxiterate 1000 --thread 8 --quiet "$O/msa_in.fna" > "$O/mafft.fna"
muscle -align "$O/msa_in.fna" -output "$O/muscle.fna" -threads 8 >/dev/null 2>&1
