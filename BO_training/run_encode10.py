import os
encode10_list = ["SRR307903","SRR307911","SRR315323","SRR315334","SRR387661","SRR534291","SRR534307","SRR534319","SRR545695","SRR545723"]
for sample_name in encode10_list:
    cmd = "python main.py -p stringtie --max_iters 210 --n_trials 1 --save_path output_stringtie_"+ \
        sample_name+"_warmup_new --scallop_bam /data008/users/zyan/encode10/"+ \
        sample_name+"/hisat.sort.bam --ard -a thompson --n_init 10 --cawarmup 60 --scallop_path /data008/users/zyan/software/stringtie-2.2.1.Linux_x86_64/ --bamid "+ \
        sample_name+" --ref_file /data008/users/zyan/encode10/GRCh38.gtf  1> stringtie_"+sample_name+".out 2>stringtie_"+sample_name+".err &"
    print(cmd)
    os.system(cmd) 