#! /usr/bin/bash

for FILE in encode10/*; do
{
	FILENAME=`basename $FILE`	
	VAR = `python main.py -p Stringtie --n_trials 1 --n_init 3 --scallop_bam $FILENAME --ard -a thompson --scallop_path /data008/users/zyan/software/stringtie-2.2.1.Linux_x86_64/ --cawarmup 2 --ref_file /biodb/human/gencode/v35/gene_annotations.gtf --max_iters 5 2> autopar_str.err 1> autopar_str.out &`
	echo $VAR
}&
done