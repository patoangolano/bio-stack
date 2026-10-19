d=$(ls -d ~/nf-work/filogenia/work/81/e3fdfb* | head -1)
ls -la --time-style=+%H:%M "$d/concatenado.model.gz" "$d/concatenado.log"
zcat "$d/concatenado.model.gz" 2>/dev/null | grep -c "" 
zcat "$d/concatenado.model.gz" 2>/dev/null | tail -2 | cut -c1-150
