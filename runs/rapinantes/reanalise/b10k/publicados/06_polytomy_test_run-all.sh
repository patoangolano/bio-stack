#!/bin/bash
#SBATCH --job-name="ASTRAL"
#SBATCH --output="astral.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH -t 24:00:00
#SBATCH -A uot138
#SBATCH --ntasks-per-node=128
#SBATCH --mem=248G


date

export HM=/expanse/lustre/scratch/smirarab/temp_project/avian/expanse_analyses
cd /expanse/lustre/scratch/smirarab/temp_project/avian/expanse_analyses/polytomytest

p=0
for GTz in `ls -S  random/*gz all_loci/*gz|tac`; do 
	for ST in 63430_filter_mRNA_random_region.gene.trees.aLRT0.05.species.tre astral-128233_genome-wide.05aLRT.gene.trees.tre astral-159205_genome-wide.05aLRT.gene.trees.tre; do
		if [ -s "$ST-${GTz/\//-}-polytomy.tre" ]; then 
			echo done $ST-${GTz/\//-}-polytomy.tre
		else
	GT=`mktemp`
	gunzip -c $GTz > $GT
			#java -Xmx50g -Djava.library.path=~/workspace/ASTRAL/Astral/lib/ -jar ~/workspace/ASTRAL/Astral/astral.5.7.4.jar  -i $GT   -o $ST-${GTz/\//-}-polytomy.tre -q $ST -p 10 &  pids[$p]=$!
			java -Xmx50g -Djava.library.path=$HM/Astral/lib -jar $HM/Astral/astral.5.15.1.jar -T 31  -C -i $GT   -o $ST-${GTz/\//-}-polytomy.tre -q $ST -p 10 &  pids[$p]=$!
		        p=$((p+1))
			if [ $(( p % 4 )) == 0 ]; then
				echo PIDs $pids
				for pid in ${pids[*]}; do
					wait $pid
				done
				p=0
				pids=()
				rm /tmp/tmp.*
			fi
		fi
	done
done	

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done

date
