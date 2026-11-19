#!/usr/bin/env bash
set -euo pipefail
L4="/mnt/c/Users/user/bio-stack/runs/enamel_relics/level4"
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/phaseD/iqtree_asr"
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -te "$O/tree_T1.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_mafft_T1" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -te "$O/tree_T2.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_mafft_T2" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -te "$O/tree_T3.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_mafft_T3" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -te "$O/tree_T1.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_muscle_T1" -redo -nt 4 --quiet
