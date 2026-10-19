L=$(ls -td ~/nf-work/filogenia/*/ | head -1)
echo "$L"; grep -cE "Task completed.*MAFFT_ALIGN" $L/.nextflow.log; grep -cE "Task completed.*CLIPKIT" $L/.nextflow.log; grep -cE "Task completed.*IQTREE_GENES" $L/.nextflow.log
grep -E "ERROR|exit: [1-9]" $L/.nextflow.log | head -5
uptime; free -g | head -2
ps -eo pcpu,pmem,etime,args --sort=-pcpu | head -8 | cut -c1-150
