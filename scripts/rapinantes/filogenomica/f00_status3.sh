d=$(ls -d ~/nf-work/filogenia/work/81/e3fdfb* | head -1)
wc -l "$d/concatenado.log"; tail -12 "$d/concatenado.log" | cut -c1-160
ls -la "$d" | head -20
top -b -n1 | head -15 | tail -9
