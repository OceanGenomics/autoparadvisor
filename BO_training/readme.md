Generalized version by Zhiwen.

All parameter bounds and options for software running (including main software and evaluation software) are moved into YAML config files.

YAML config files: stringtie.yml and scallop.yml should be placed under autoparadvisor directory. (test version, for now)

Additional packages required:
-YAML library for python
    pip install pyyaml
-YAML module for perl
    install YAML
    install YAML::XS
    tips: could use cpan for perl module installation

example usage:
python main.py -p stringtie --n_trials 1 --n_init 10 --max_iters 150 --cawarmup 20 --bam_file /data008/users/zyan/sample_data/SRR5882130/star/SRR5882130.sort.bam --ard -a thompson --software_path /data008/users/zyan/software/stringtie-2.2.1.Linux_x86_64/ --ref_file /biodb/human/gencode/v35/gene_annotations.gtf 2> autopar_str.err 1> autopar_str.out &
