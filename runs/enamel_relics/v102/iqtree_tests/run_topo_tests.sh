#!/usr/bin/env bash
set -euo pipefail
L4="/mnt/c/Users/user/bio-stack/runs/enamel_relics/level4"
D101="/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/phaseD/iqtree_asr"
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/v102/iqtree_tests"
# Formal topology tests (KH, SH, wKH, wSH, ELW, AU with 10,000 RELL bootstrap replicates)
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -z "$O/trees_T1_T2_T3.nwk" -te "$D101/tree_T1.nwk" -m GTR+G -zb 10000 -zw -au -o Amis,Ggan,Cpor -pre "$O/topo_mafft" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -z "$O/trees_T1_T2_T3.nwk" -te "$D101/tree_T1.nwk" -m GTR+G -zb 10000 -zw -au -o Amis,Ggan,Cpor -pre "$O/topo_muscle" -redo -nt 4 --quiet
# Also run MUSCLE T2 and T3 ASR for completeness
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -te "$D101/tree_T2.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_muscle_T2" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -te "$D101/tree_T3.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_muscle_T3" -redo -nt 4 --quiet
