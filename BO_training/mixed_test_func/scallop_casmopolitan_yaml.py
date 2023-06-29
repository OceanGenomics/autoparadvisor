import numpy as np
from test_funcs.base import TestFunction
import pdb,sys,os,csv,glob
from multiprocessing import Manager
import multiprocessing
import subprocess
import yaml


def Scallop_base(index,x,Result,ref_file,num_transcripts,bam_file,library_type,software_path,docs,problem):
    pid = os.getpid()
    #NOTE: original scallop_bounds list is moved into YAML file
    #initial command of choosen assemble software
    software = docs['initial_option']
    #initial parameters for choosen software
    parameter_bounds = docs['parameter_bounds']
    cmd = software_path + software['name'] + software['input_option'] + \
         bam_file + software['additional_option']

    #for optional parameter choosing
    for i in range(x.shape[0]):
        #get parameter instance from yaml dict
        parameter = parameter_bounds[i]
        parameter_name = list(parameter.keys())[0]
        parameter_type = parameter[parameter_name]['type']

        #for category parameters
        #NOTE: for now binary parameters are divided into two cases
        #TF stands for True/False usage (in scallop)
        #turn_on stands for use/not use usage (in stringtie)
        if parameter_type=='cag': 
            if parameter[parameter_name]['usage']=='TF':
                if int(x[i])==0:
                    cmd += parameter[parameter_name]['prefix'] + parameter_name + " false"
                else:
                    cmd += parameter[parameter_name]['prefix']  + parameter_name + " true"
            elif parameter[parameter_name]['usage']=='turn_on':
                if int(x[i])==1:
                    cmd += parameter[parameter_name]['prefix'] + parameter_name
        #for continous type parameter (int)
        elif parameter_type=='int':
            cmd += parameter[parameter_name]['prefix']  + parameter_name + " " + str(int(x[i]))
        #for continous type parameter (float)
        elif parameter_type=='float':
            cmd += parameter[parameter_name]['prefix'] + parameter_name + " " + str(x[i])

    #NOTE: comment out subsampling since its not used
    #if(subsamp<1): 
    #    cmd += " --subsampling " + str(subsamp)
    #NOTE: library_type is moved into YAML config
    #cmd +=" --library_type " + library_type
    #commands for output, should be same within other software
    #TODO: maybe need to generalize later
    cmd +=" -o ./" + str(pid) + ".gtf" + " > /dev/null 2>&1"
    print(f"Run scallop with the following command: \n {cmd}")
    os.system(cmd)
    #pdb.set_trace()
    #check if output gft is started with chr
    #if the chromosome head starts with 'chr', remove it
    cmd = "grep -c '^chr' "+ str(pid) + ".gtf"
    chr_header = int(subprocess.getoutput(cmd))
    ref_cmd = "grep -c '^chr' "+ ref_file
    ref_chr_header = int(subprocess.getoutput(ref_cmd))
    print(f'number of transcript start with chr: {chr_header}')
    #NOTE: remove 'chr' if ref genome doesn't start with chr
    if ref_chr_header==0:
        cmd = "sed -i 's/^chr//' " + str(pid) + ".gtf"
        #cmd = "cut -c4- " + str(pid) + ".gtf > " + str(pid) + "_new.gtf"
        #print(cmd)
        os.system(cmd)
        #cmd = "mv " + str(pid) + "_new.gtf " + str(pid) + ".gtf"
        #os.system(cmd)
    #NOTE: evaluation program running command is also moved into YAML
    cmd = docs['gffcompare']['directory'] + docs['gffcompare']['command'] + \
         ref_file + ' ./' + str(pid) + ".gtf"
    print("Run gffcompare: \n")
    print(cmd)
    os.system(cmd)
    cmd = docs['gtfcuff']['directory'] + docs['gtfcuff']['command'] + \
         './gffcmp.' + str(pid) + ".gtf" + ".tmap " + str(num_transcripts)
    print("Run gtfcuff: \n")
    print(cmd)
    #pdb.set_trace()
    result = subprocess.getoutput(cmd)
    print(result)
    #TODO: get AUC from gtfcuff stdout, this part still needs more generalization
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
    #TODO: change this later, now only works for mixed(cag+cont+int)
    problem_type = 'mixed'
    def __init__(self, bam_file, normalize=False,boundary_fold = 0,ref_file='',library_type = 'empty',problem=''):
        super(Scallop,self).__init__(normalize)
        assert boundary_fold>=0
        self.bam_file = bam_file
        self.ref_file = ref_file
        self.library_type = library_type
        self.problem = problem

        #NOTE: read in software usage and parameter from yaml file
        #TODO:  read-in file name is still hard coded
        # YAML file store path can be further discussed
        path = os.path.abspath(os.path.join(os.getcwd(),".."))
        with open(path+'/stringtie.yml', 'r') as file:
            docs = yaml.safe_load(file)
            self.docs = docs
        #software contains a dict of basic use of specific parameters
        #parameter_bounds contain a list of tunable parameters
        software = docs['initial_option']
        parameter_bounds = docs['parameter_bounds']
        self.software = software
        self.parameter_bounds = parameter_bounds


        #set up parameter information
        self.boundary_fold = boundary_fold
        self.para_to_index = {}
        default = []
        categorical_dims = []
        continuous_dims = []
        int_constrained_dims = []
        hard_lb = []
        hard_ub = []
        
        for i in range(len(parameter_bounds)):
            parameter = parameter_bounds[i]
            parameter_name = list(parameter.keys())[0]
            parameter_type = parameter[parameter_name]['type']
            default.append(parameter[parameter_name]['default'])
            self.para_to_index[parameter_name] = i

            if (parameter_type=='cag'):
                categorical_dims.append(i)
            else:
                hard_lb.append(parameter[parameter_name]['hard_min'])
                hard_ub.append(parameter[parameter_name]['hard_max'])
                continuous_dims.append(i)
                if(parameter_type=='int'):
                    int_constrained_dims.append(i)
    
        self.categorical_dims = np.asarray(categorical_dims)
        self.continuous_dims = np.asarray(continuous_dims)
        self.default = np.array(default)
        self.dim = len(self.categorical_dims) + len(self.continuous_dims)
        self.int_constrained_dims = np.asarray(int_constrained_dims)
        self.hard_lb = np.array(hard_lb)
        self.hard_ub = np.array(hard_ub)
        
        ##specify the domain boundary
        lb = []
        ub = []
        for i in range(len(parameter_bounds)):
            parameter = parameter_bounds[i]
            parameter_name = list(parameter.keys())[0]
            parameter_type = parameter[parameter_name]['type']
            if (parameter_type!='cag'):
                if (boundary_fold==0):
                    lb.append(parameter[parameter_name]['min'])
                    ub.append(parameter[parameter_name]['max'])
                else:
                    lb.append(max(parameter[parameter_name]['hard_min'], \
                        (1-boundary_fold)*parameter[parameter_name]['default']))
                    ub.append(min(parameter[parameter_name]['hard_max'], \
                        (1+boundary_fold)*parameter[parameter_name]['default']))
        self.lb = np.array(lb)
        self.ub = np.array(ub)
        
        #TODO: now only support binary cag parameters
        #original code: self.vertices = np.array([2,2])
        self.n_vertices = np.array(len(categorical_dims)*[2])
        self.config = self.n_vertices
        
        #self.lamda = lamda
        #no normalize implementation
        self.mean = None
        self.std = None
        
        #NOTE: get number of reference transcript for future gftcuff and ca_warmup
        cmd = 'cat ' + self.ref_file + ' | awk \'{print $3}\' | grep -c transcript'
        self.num_transcripts = int(subprocess.check_output(cmd,shell=True))
        print(f'run {cmd}, get {self.num_transcripts} ref transcripts')
        #self.num_transcripts = 197649

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
                tmp_process = multiprocessing.Process(target=Scallop_base, \
                    args=(i,X[i],Y,self.ref_file,self.num_transcripts,self.bam_file,self.library_type,scallop_path,self.docs,self.problem))
                process_list.append(tmp_process)
            for process in process_list:
                process.start()
            for process in process_list:
                process.join()
            Y = list(Y)
        # return the negative AUC score as accuracy
        return -np.array(Y)

    #NOTE: make this function as Scallop method
    def read_warmup_info(self, path_name):
        #pdb.set_trace()
        paraname_to_index = {}
        for i in range(len(self.parameter_bounds)):
            parameter = self.parameter_bounds[i]
            parameter_name = list(parameter.keys())[0]
            paraname_to_index[parameter_name] = i
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
