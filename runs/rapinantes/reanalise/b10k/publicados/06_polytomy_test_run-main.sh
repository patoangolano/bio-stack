#!/bin/bash
#SBATCH --job-name="ASTRAL"
#SBATCH --output="astral.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH -t 6:00:00
#SBATCH -A uot138
#SBATCH --ntasks-per-node=16


export SLURM_EXPORT_ENV=ALL

date

cd /expanse/lustre/scratch/smirarab/temp_project/avian/expanse_analyses/polytomytest
# cd /home/smirarab/projects/avian/avian-data/ASTRAL/polytomytests

GTz=128233_genome-wide.05aLRT.gene.trees.gz
GT=`mktemp`
gunzip -c $GTz > $GT
ST=astral-128233_genome-wide.05aLRT.gene.trees.tre
		java -Xmx32g  -jar ../Astral/astral.5.15.1.jar  -i $GT   -o $ST-${GTz/\//-}-polytomy.tre -q $ST -t 10 
GTz=159205_genome-wide.05aLRT.gene.trees.gz
GT=`mktemp`
gunzip -c $GTz > $GT
ST=astral-159205_genome-wide.05aLRT.gene.trees.tre
		java -Xmx32g  -jar ../Astral/astral.5.15.1.jar  -i $GT   -o $ST-${GTz/\//-}-polytomy.tre -q $ST -t 10 
date

exit

GTz=63430_filter_mRNA_random_region.gene.trees.aLRT0.05.gene.trees.gz
GT=`mktemp`
gunzip -c $GTz > $GT
ST=63430_filter_mRNA_random_region.gene.trees.aLRT0.05.species.tre
		java -Xmx32g  -jar ../Astral/astral.5.15.1.jar  -i $GT   -o $ST-${GTz/\//-}-polytomy.tre -q $ST -t 10 


