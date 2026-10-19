uptime; free -g | head -2
echo "iqtree rodando: $(pgrep -fc 'iqtree -m MFP')"
ps -eo pcpu,pmem,etime,args --sort=-pcpu | grep -v "iqtree -m MFP" | head -6 | cut -c1-140
