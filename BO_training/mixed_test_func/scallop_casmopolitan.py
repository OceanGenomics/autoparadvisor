import numpy as np
from test_funcs.base import TestFunction
import pdb,sys,os,csv,glob
from multiprocessing import Manager
import multiprocessing
import subprocess

scallop_bounds = [{"name":"uniquely_mapped_only","type":"int","min":0,"max":1,"default":0,"hard_min":0,"hard_max":1},
                {"name":"use_second_alignment","type":"int","min":0,"max":1,"default":0,"hard_min":0,"hard_max":1},
                {"name":"max_dp_table_size","type":"int","min":0,"max":100000,"default":10000,"hard_min":0,"hard_max":float('inf')},
               {"name":"max_edit_distance","type":"int","min":0,"max":100,"default":10,"hard_min":0,"hard_max":float('inf')},
               {"name":"max_num_exons","type":"int","min":0,"max":10000,"default":1000,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_bundle_gap","type":"int","min":0,"max":500,"default":50,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_exon_length","type":"int","min":0,"max":200,"default":20,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_flank_length","type":"int","min":0,"max":30,"default":3,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_mapping_quality","type":"int","min":0,"max":10,"default":1,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_num_hits_in_bundle","type":"int","min":0,"max":200,"default":20,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_router_count","type":"int","min":0,"max":10,"default":1,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_splice_boundary_hits","type":"int","min":0,"max":10,"default":1,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_subregion_gap","type":"int","min":0,"max":30,"default":3,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_subregion_length","type":"int","min":0,"max":150,"default":15,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_transcript_length_base","type":"int","min":0,"max":1500,"default":150,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_transcript_length_increase","type":"int","min":0,"max":500,"default":50,"hard_min":0,"hard_max":float('inf')},
               {"name":"max_intron_contamination_coverage","type":"float","min":0.0,"max":20.0,"default":2.0,"hard_min":0,"hard_max":float('inf')},
               {"name":"min_subregion_overlap","type":"float","min":0.0,"max":15.0,"default":1.5,"hard_min":0,"hard_max":float('inf')}           
]


'''
Scallop_para_to_index = {"max_dp_table_size":0,"max_edit_distance":1,"max_intron_contamination_coverage":2,
"max_num_exons":3,"min_bundle_gap":4,"min_exon_length":5,"min_flank_length":6,"min_mapping_quality":7,
"min_num_hits_in_bundle":8,"min_router_count":9,"min_splice_boundary_hits":10,"min_subregion_gap":11,
"min_subregion_length":12,"min_subregion_overlap":13,"min_transcript_length_base":14,"min_transcript_length_increase":15,
"uniquely_mapped_only":16,"use_second_alignment":17
}


def Scallop(X,ref_file,num_transcripts,bam_file):
    #pdb.set_trace()
    with Manager() as manager:
        Y = manager.list(range(X.size(0)))
        process_list = []
        for i in range(X.size()[0]):
            tmp_process = multiprocessing.Process(target=Scallop_base,args=(i,X[i],Y,ref_file,num_transcripts,bam_file,))
            process_list.append(tmp_process)
        for process in process_list:
            process.start()
        for process in process_list:
            process.join()
        Y = list(Y)
    return torch.tensor(Y,dtype=dtype, device=device)

'''


def Scallop_base(index,x,Result,ref_file,num_transcripts,bam_file,library_type,scallop_path='',subsamp=1):
    #pdb.set_trace()
    pid = os.getpid()
    if(scallop_path==''):
        cmd = "scallop -i " + bam_file + " --verbose 0 --min_transcript_coverage 0"
    else:
        cmd = scallop_path + "scallop -i " + bam_file + " --verbose 0 --min_transcript_coverage 0"
    for i in range(x.shape[0]):
        if(scallop_bounds[i]["name"]=="uniquely_mapped_only" or scallop_bounds[i]["name"]=="use_second_alignment"):
            if(int(x[i])==0):
                cmd += " --" + scallop_bounds[i]["name"] + " false"
            else:
                cmd += " --" + scallop_bounds[i]["name"] + " true"
            continue
        if(scallop_bounds[i]["type"]=="int"):
            cmd += " --" + scallop_bounds[i]["name"] + " " + str(int(x[i]))
        else:
            cmd += " --" + scallop_bounds[i]["name"] + " " + str(x[i])
        '''
        if(var_types[i]=='float'):
            cmd += " --" + Scallop_index_to_para[i] + " " + str(x[i].item())
        else:
            cmd += " --" + Scallop_index_to_para[i] + " " + str(x[i].int().item())
        '''
    if(subsamp<1):
        cmd += " --subsampling " + str(subsamp)
    cmd +=" --library_type " + library_type
    cmd +=" -o ./" + str(pid) + ".gtf" + " > /dev/null 2>&1"
    print("Run scallop with the following command: \n")
    print(cmd)
    os.system(cmd)
    #pdb.set_trace()
    #if the chromosome head starts with 'chr', remove it
    #cmd = "head -c 3 " + str(pid) + ".gtf"
    cmd = "grep -c '^chr' "+ str(pid) + ".gtf"
    chr_header = int(subprocess.getoutput(cmd))
    print(chr_header)
    if(chr_header>0):
        #cmd = "sed -i 's/^chr//' " + str(pid) + ".gtf"
        #cmd = "cut -c4- " + str(pid) + ".gtf > " + str(pid) + "_new.gtf"
        #print(cmd)
        #os.system(cmd)
        #cmd = "mv " + str(pid) + "_new.gtf " + str(pid) + ".gtf"
        #os.system(cmd)
    cmd = "/data008/users/zyan/software/gffcompare/gffcompare -r " + ref_file + ' ./' + str(pid) + ".gtf"
    print("Run gffcompare: \n")
    print(cmd)
    os.system(cmd)
    cmd = "/data008/users/zyan/software/rnaseqtools-1.0.3/gtfcuff/gtfcuff auc ./" + 'gffcmp.' + str(pid) + ".gtf" + ".vim " + str(num_transcripts)
    print("Run gtfcuff: \n")
    print(cmd)
    #pdb.set_trace()
    result = subprocess.getoutput(cmd)
    print(result)
    if(len(str(result))<10):
        auc_val = 0.0
    else:
        auc_val = float(str(result).split("auc")[-1].split("=")[-1].split("\\")[0])
    print(auc_val)
    Result[index] = auc_val
    cmd = 'rm '+ str(pid) + ".gtf"
    os.system(cmd)
    cmd = 'rm *.' + str(pid) + '.*'
    os.system(cmd)



class Scallop(TestFunction):
    problem_type = 'mixed'
    def __init__(self, bam_file, normalize=False,boundary_fold = 0,ref_file='',library_type = 'empty'):
        super(Scallop,self).__init__(normalize)
        assert boundary_fold>=0
        self.para_to_index = {}
        self.boundary_fold = boundary_fold
        self.categorical_dims = np.array([0, 1])
        self.continuous_dims = np.array([i for i in range(2,18)])
        self.dim = len(self.categorical_dims) + len(self.continuous_dims)
        #self.default = []
        self.n_vertices = np.array([2, 2])
        self.config = self.n_vertices
        #pdb.set_trace()
        for i in range(2):
            self.default.append(scallop_bounds[i]['default'])
            self.para_to_index[scallop_bounds[i]['name']] = i
        self.int_constrained_dims = np.array([i for i in range(2,16)])
        #pdb.set_trace()
        hard_lb = []
        hard_ub = []
        for i in range(2,18):
            hard_lb.append(scallop_bounds[i]['hard_min'])
            hard_ub.append(scallop_bounds[i]['hard_max'])
            self.para_to_index[scallop_bounds[i]['name']] = i

        self.hard_lb = np.array(hard_lb)
        self.hard_ub = np.array(hard_ub)
        lb = []
        ub = []
        #pdb.set_trace()
        if(boundary_fold==0): #specify the domain boundary
            for i in range(2,18):
                lb.append(scallop_bounds[i]['min'])
                ub.append(scallop_bounds[i]['max'])
                self.default.append(scallop_bounds[i]['default'])
        else:
            for i in range(2,18):
                lb.append(max(scallop_bounds[i]['hard_min'],(1-boundary_fold)*scallop_bounds[i]['default']))
                ub.append(min(scallop_bounds[i]['hard_max'],(1+boundary_fold)*scallop_bounds[i]['default']))
                self.default.append(scallop_bounds[i]['default'])
        self.lb = np.array(lb)
        self.ub = np.array(ub)
        #self.lamda = lamda
        #no normalize implementation
        self.mean, self.std = None, None
        self.bam_file = bam_file
        self.ref_file = ref_file
        self.library_type = library_type
        #pdb.set_trace()
        cmd = 'cat ' + self.ref_file + ' | awk \'{print $3}\' | grep -c transcript'
        #self.num_transcripts = int(subprocess.check_output(cmd,shell=True))
        self.num_transcripts = 197649
        self.default = np.array(self.default)
    def compute(self,X,normalize=False,scallop_path='',subsamp=1):
        #pdb.set_trace()
        if X.ndim == 1:
            X = X.reshape(1, -1)
        N = X.shape[0]
        #assert np.array_equal(X[:,self.int_constrained_dims], np.round(X[:,self.int_constrained_dims]))
        X[:,self.int_constrained_dims] = np.round(X[:,self.int_constrained_dims])
        with Manager() as manager:
            Y = manager.list(range(N))
            process_list = []
            for i in range(N):
                tmp_process = multiprocessing.Process(target=Scallop_base,args=(i,X[i],Y,self.ref_file,self.num_transcripts,self.bam_file,self.library_type,scallop_path,subsamp,))
                process_list.append(tmp_process)
            for process in process_list:
                process.start()
            for process in process_list:
                process.join()
            Y = list(Y)
        # return the negative AUC score as accuracy
        return -np.array(Y)
    def read_para_from_file(self,file_name):
        para_x = np.zeros(self.dim)
        counter = 0
        with open(file_name,"rt")as f:
            freader = csv.reader(f,delimiter='\t')
            for line in freader:
                if(line[0] in self.para_to_index.keys()):
                    para_x[self.para_to_index[line[0]]] = float(line[1])
                    counter+=1
        assert counter==18
        return para_x

def read_warmup_info(path_name):
    #pdb.set_trace()
    paraname_to_index = {}
    #paraname_list = []
    for i in range(len(scallop_bounds)):
        paraname_to_index[scallop_bounds[i]['name']] = i
    paraname_list = sorted(list(paraname_to_index.keys()))
    auc_files = glob.glob(path_name + "*.auc")
    X = np.zeros((0,len(paraname_list)))
    Y = []
    #pdb.set_trace()
    for auc_file in auc_files:
        result = subprocess.getoutput("cat " + auc_file)
        if(len(str(result))<10):
            Y.append(0.0)
        else:
            Y.append(-float(str(result).split("auc")[-1].split("=")[-1].split("\\")[0]))
        para_instance = auc_file.split("/")[-1]
        para_instance = para_instance[:len(para_instance)-4].split("_")[1:]
        X_instance = np.zeros((len(para_instance),))
        for i in range(len(para_instance)):
            X_instance[paraname_to_index[paraname_list[i]]] = float(para_instance[i])
        X = np.vstack((X,X_instance))
    return X,Y

def coordinate_ascent_warmup(Function,fold_step=1,iterations=20):
    #pdb.set_trace()
    max_iter_for_each_dim = int(iterations/Function.dim) + 1
    dim_shuffle = np.random.permutation(Function.dim)
    shuffle_index = 0
    ca_index = dim_shuffle[shuffle_index]
    ca_direction = 0
    old_fold_step = fold_step
    temp_fold_step = fold_step
    y_old = float(Function.compute(Function.default))
    iter_num = 0
    x_old = np.copy(Function.default)
    ca_X = np.zeros((0,Function.dim))
    ca_Y = []
    ca_X = np.vstack((ca_X,x_old))
    ca_Y.append(y_old)
    iter_num+=1
    while iter_num < iterations:
        if ca_index in Function.categorical_dims:
            if(x_old[ca_index]+1<=Function.config[ca_index]-1 and ca_direction>=0 and ca_direction<max_iter_for_each_dim):
                x_new = np.copy(x_old)
                x_new[ca_index] = x_new[ca_index] + 1
                y_new = float(Function.compute(x_new))
                iter_num+=1
                ca_X = np.vstack((ca_X,x_new))
                ca_Y.append(y_new)
                if(y_new < y_old):
                    ca_direction += 1
                    y_old = y_new
                    x_old = np.copy(x_new)
                    continue
            if(x_old[ca_index]-1>=0 and ca_direction<=0 and abs(ca_direction)<max_iter_for_each_dim):
                x_new = np.copy(x_old)
                x_new[ca_index] = x_new[ca_index] - 1
                y_new = float(Function.compute(x_new))
                iter_num+=1
                ca_X = np.vstack((ca_X,x_new))
                ca_Y.append(y_new)
                if(y_new < y_old):
                    y_old = y_new
                    x_old = np.copy(x_new)
                    ca_direction -= 1
                    continue
        else:
            if(x_old[ca_index] + temp_fold_step*Function.default[ca_index] <=Function.hard_ub[ca_index-len(Function.categorical_dims)] and ca_direction>=0 and ca_direction<max_iter_for_each_dim):
                x_new = np.copy(x_old)
                x_new[ca_index] = x_new[ca_index] + temp_fold_step*Function.default[ca_index]
                #x_new[ca_index] = x_new[ca_index] * (1+temp_fold_step)
                if(ca_index in Function.int_constrained_dims):
                    x_new[ca_index] = int(x_new[ca_index])
                y_new = float(Function.compute(x_new))
                iter_num+=1
                ca_X = np.vstack((ca_X,x_new))
                ca_Y.append(y_new)
                if(y_new < y_old):
                    temp_fold_step = temp_fold_step/2.0
                    ca_direction += 1
                    y_old = y_new
                    x_old = np.copy(x_new)
                    continue
            if(x_old[ca_index] - temp_fold_step*Function.default[ca_index]>=Function.hard_lb[ca_index-len(Function.categorical_dims)] and ca_direction<=0 and abs(ca_direction)<max_iter_for_each_dim):
                x_new = np.copy(x_old)
                x_new[ca_index] = x_new[ca_index] - temp_fold_step*Function.default[ca_index]
                #x_new[ca_index] = x_new[ca_index] * (1-temp_fold_step)
                if(ca_index in Function.int_constrained_dims):
                    x_new[ca_index] = int(x_new[ca_index])
                y_new = float(Function.compute(x_new))
                iter_num+=1
                ca_X = np.vstack((ca_X,x_new))
                ca_Y.append(y_new)
                if(y_new < y_old):
                    temp_fold_step = temp_fold_step/2.0
                    ca_direction -= 1
                    y_old = y_new
                    x_old = np.copy(x_new)
                    continue
        shuffle_index = (shuffle_index+1) % Function.dim
        #ca_index = (ca_index+1) % Function.dim
        ca_index = dim_shuffle[shuffle_index]
        ca_direction = 0
        if(shuffle_index==0):
            old_fold_step = old_fold_step/2.0
        temp_fold_step = old_fold_step
    return (ca_X,ca_Y)
