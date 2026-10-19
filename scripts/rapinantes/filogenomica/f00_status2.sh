d=$(ls -d ~/nf-work/filogenia/work/81/e3fdfb* | head -1)
grep -E "RAM|BEST SCORE|Iteration [0-9]+ /|Wall-clock time for ModelFinder|NOTE" "$d/concatenado.log" | tail -4 | cut -c1-150; tail -1 "$d/concatenado.log" | cut -c1-150
