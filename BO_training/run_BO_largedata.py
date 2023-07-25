import os,csv,subprocess,multiprocessing,glob
def run_BO_on_stringtie(SRR_list,SRR_path_list):
    assert len(SRR_list)==len(SRR_path_list)
    for i in range(len(SRR_list)):
        #cmd = 'python main.py -p stringtie --max_iters 210 --n_trials 1 --save_path output_stringtie_'+SRR_list[i]+'_warmup_new --bam '+SRR_path_list[i]+
        # 'Aligned.sortedByCoord.out.bam --ard -a thompson --ref GRCh38 --n_init 10 --cawarmup 60 --bamid '+SRR_list[i]+" 1> stringtie_"+SRR_list[i]+".out 2>stringtie_"+SRR_list[i]+".err"
        cmd = "python main.py -p stringtie --max_iters 210 --n_trials 1 --save_path output/output_stringtie_"+ \
        SRR_list[i]+"_warmup_new --scallop_bam " + \
        SRR_path_list[i]+"/Aligned.sortedByCoord.out.bam --ard -a thompson --n_init 10 --cawarmup 60 --scallop_path /data008/users/zyan/software/stringtie-2.2.1.Linux_x86_64/ --bamid "+ \
        SRR_list[i]+" --ref_file /data008/users/zyan/encode10/GRCh38.gtf  1> output/stringtie_"+SRR_list[i]+".out 2>output/stringtie_"+SRR_list[i]+".err"
        print(cmd)
        os.system(cmd)

def main_process(chunk_size,chunk_start,chunk_end):
    SRR_list_all = []
    SRR_path_list_all = []
    process_list = []
    with open("/data008/users/zyan/data_large/set_200_info.txt","rt")as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            SRR_list_all.append(line[0])
            SRR_path_list_all.append(line[1])
    for i in range(chunk_start,chunk_end):
        #print(i)
        #print('SRR_list_all:',SRR_list_all[i*chunk_size:(i+1)*chunk_size])
        #print('SRR_path_list_all:',SRR_path_list_all[i*chunk_size:(i+1)*chunk_size])
        tmp_process = multiprocessing.Process(target=run_BO_on_stringtie,args=(SRR_list_all[i*chunk_size:(i+1)*chunk_size],SRR_path_list_all[i*chunk_size:(i+1)*chunk_size],))
        process_list.append(tmp_process)
    for process in process_list:
        process.start()
    for process in process_list:
        process.join()

main_process(20,0,10)