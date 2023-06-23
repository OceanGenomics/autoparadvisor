### bash command line to generate .msh file from fastq files 
cat SRR1583599.1_1.fastq.gz SRR1583599.1_2.fastq.gz | mash sketch -r -m 2 -o SRR1583599 -