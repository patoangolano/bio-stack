#!/usr/bin/env bash
# status das tarefas longas: log do IQTREE_CONCAT (AA) e contagem do DNA12
d=$(ls -d ~/nf-work/filogenia/work/81/e3fdfb* 2>/dev/null | head -1)
echo "concat AA: $d"
[ -n "$d" ] && { grep -E "NOTE: .*RAM|memory|Merging|BEST SCORE|Iteration|ModelFinder|Wall-clock|Estimate model|Subset" "$d/concatenado.log" 2>/dev/null | tail -4; tail -2 "$d/concatenado.log"; }
free -g | sed -n 2p
