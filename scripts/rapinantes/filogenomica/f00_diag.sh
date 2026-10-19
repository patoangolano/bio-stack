uptime; free -g | head -2
ls -la ~/rapinantes/raw/GCA_026930845.1/ ~/rapinantes/miniprot/ | head -30
tail -5 ~/rapinantes/miniprot/Pandion_haliaetus.miniprot.log 2>/dev/null
dmesg 2>/dev/null | tail -5
ps -eo pid,pcpu,pmem,etime,args --sort=-pmem | head -8 | cut -c1-180
