import numpy as np
from test_funcs.base import TestFunction
import pdb,sys,os,csv,glob
from multiprocessing import Manager
import multiprocessing
import subprocess
import yaml
import copy
import threading

def Scallop_base(index,x,Result,ref_file,input_file,software_path,docs):
    pid = os.getpid()
    #NOTE: original scallop_bounds list is moved into YAML file
    #initial parameters for choosen software
    parameter_bounds = docs['parameter_bounds']
    pcmd = ''

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
                    pcmd = ' '.join([pcmd, parameter[parameter_name]['prefix'], "false"])
                else:
                    pcmd = ' '.join([pcmd, parameter[parameter_name]['prefix'], "true"])
            elif parameter[parameter_name]['usage']=='turn_on':
                if int(x[i])==1:
                    pcmd = ' '.join([pcmd, parameter[parameter_name]['prefix']])
        #for continous type parameter (int)
        elif parameter_type=='int':
            pcmd = ' '.join([pcmd, parameter[parameter_name]['prefix'], str(int(x[i]))])
        #for continous type parameter (float)
        elif parameter_type=='float':
            pcmd = ' '.join([pcmd, parameter[parameter_name]['prefix'], str(x[i])])

    #assemble command of choosen assemble software with formatting
    format_dict = {'input_file':input_file, 'parameters':pcmd, 'path':software_path}
    software = docs['testing_software']
    for option in software:
        if option == 'format':
            software_format = software[option]
        elif option == 'pid':
            format_dict[option] = pid
        else:
            format_dict[option] = software[option]
    cmd = software_format.format(**format_dict)
    print(f"Run scallop with the following command: \n {cmd}")
    os.system(cmd)
    #pdb.set_trace()

    #NOTE this chr part may be removed
    #check if output gft is started with chr, if remove it
    cmd = docs['precheck']['check_command'].format(pid)
    chr_header = int(subprocess.getoutput(cmd))
    print(f'number of transcript start with chr: {chr_header}')
    #NOTE: remove 'chr' if ref genome doesn't start with chr
    if chr_header>0:
        cmd = docs['precheck']['excute_command'].format(pid)
        print(cmd)
        os.system(cmd)

    # looping all required evaluation steps
    for val_step in docs['evaluation']:
        format_dict = {}
        for key in val_step:
            if key == 'format':
                eval_format = val_step[key]
            elif key == 'pid':
                format_dict[key] = pid
            else:  
                format_dict[key] = val_step[key]
        cmd = eval_format.format(**format_dict)
        print("Run evaluation step: \n")
        print(cmd)
        os.system(cmd)

    #pdb.set_trace()
    cmd = docs['postcheck']['auc_command'].format(pid)
    auc_val = subprocess.getoutput(cmd)
    print(auc_val)
    Result[index] = 0.0 if auc_val == '' else float(auc_val)

    # remove files with pid
    cmd = docs['postcheck']['clear_command'].format(pid)
    os.system(cmd)



class Scallop(TestFunction):
    #TODO: change this later, now only works for mixed(cag+cont+int)
    problem_type = 'mixed'
    def __init__(self, input_file, normalize=False, boundary_fold = 0,ref_file=''):
        super(Scallop,self).__init__(normalize)
        assert boundary_fold>=0
        self.input_file = input_file
        self.ref_file = ref_file

        #NOTE: read in software usage and parameter from yaml file
        #TODO:  read-in file name is still hard coded
        # YAML file store path can be further discussed
        path = os.path.abspath(os.path.join(os.getcwd(),".."))
        with open(path+'/stringtie.yml', 'r') as file:
            docs = yaml.safe_load(file)
            self.docs = docs
        #software contains a dict of basic use of specific parameters
        #parameter_bounds contain a list of tunable parameters
        parameter_bounds = docs['parameter_bounds']
        self.parameter_bounds = parameter_bounds


        #set up parameter information
        self.boundary_fold = boundary_fold
        self.para_to_index = {}
        self.step_size = {}
        self.para_type = {}
        self.parameter_values = {}
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
            self.parameter_values[parameter_name] = parameter[parameter_name]['default']
            self.step_size[parameter_name] = parameter[parameter_name]['step']
            self.para_type[parameter_name] = parameter_type
            self.para_to_index[parameter_name] = i

            if (parameter_type=='cag'):
                categorical_dims.append(i)
            else:
                hard_lb.append(float(parameter[parameter_name]['hard_min']))
                hard_ub.append(float(parameter[parameter_name]['hard_max']))
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
                    lb.append(max(float(parameter[parameter_name]['hard_min']), \
                        (1-boundary_fold)*parameter[parameter_name]['default']))
                    ub.append(min(float(parameter[parameter_name]['hard_max']), \
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
        

    def compute(self,X,normalize=False,software_path=''):
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
                    args=(i,X[i],Y,self.ref_file,self.input_file,software_path,self.docs))
                process_list.append(tmp_process)
            for process in process_list:
                process.start()
            for process in process_list:
                process.join()
            Y = list(Y)
        # return the negative AUC score as accuracy
        return -np.array(Y)

    #NOTE: not used since Perl script is discard
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
    '''
    written by Yihang
    '''
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

def coordinate_ascent_warmup_yaml(f, max_iters=60,path='',num_threads = 1):
    '''
    acheived same function of Perl script
    f: blackbox function, e.g Scallop()
    max_iters: max iteration allowed from ca_warmup
    path: system path to blackbox software
    num_threads: controls the number of parallel steps to take in each direction
    '''
    # get same variables as perl script
    parameter_values = copy.deepcopy(f.parameter_values)
    step_size = copy.deepcopy(f.step_size)
    type = copy.deepcopy(f.para_type)

 
    def check_with_one_change(param_to_change, param_value, check):
        # get param position in unsorted parameter list
        pos = f.para_to_index[param_to_change]
        # change param and check if within bounds
        if check == 'check':
            # cag param should between 0-1
            if param_value!='' and type[param_to_change]=='cag': 
                if param_value>1 or param_value<0:
                    return 0
            # int and float param should between hard lower bound and hard upper bound
            elif param_value!='' and type[param_to_change]!='cag':
                if float(param_value) > f.hard_ub[pos] or float(param_value) < f.hard_lb[pos]:
                    return 0
            # if got here and in check mode, means check is passed
            return 1
        # run software
        else:
            # ensures the floating point parameters don't get too complex
            if type[param_to_change] == 'float':
                param_value = float(format(param_value, '.2f'))
            # use index position to find the 'parameter to change'
            x_new = [param_value if p==param_to_change else parameter_values[p] for p in list(type.keys())]
            return x_new

    # initialize and gets default AUC
    print('get default auc in ca_warmup')
    ca_X = np.zeros((0,f.dim))
    ca_X = np.vstack((ca_X,f.default))
    cur_auc = float(f.compute(f.default, normalize=f.normalize, software_path=path))
    ca_Y = []
    ca_Y.append(cur_auc)

    # loops as long as you were able to decrease one of the step sizes
    print('starting ca_warmup')
    decreased_steps = 1
    while decreased_steps == 1:
        decreased_steps = 0
        made_one_change = 1

        # loops as long as you have made a change in the parameter vector
        # without a change in step size
        while made_one_change == 1:
            made_one_change = 0
            for param in sorted(list(type.keys())): 
                print("best current auc:",cur_auc)

                # loops as long as you can continue moving on this parameter (coordinate) and still increase AUC
                single_param_change = 1
                while single_param_change == 1:
                    single_param_change = 0

                    iter_num = len(ca_Y)
                    print('iterarion number: ',iter_num)
                    if max_iters > 0 and iter_num > max_iters: return (ca_X,ca_Y)

                    # non-boolen parameters are run with different step size (1,...,t)
                    sys.stderr.write(f"Updating {param}, type: {type[param]} \n")
                    if type[param] != "bool":
                        x_new_plus = np.zeros((0,f.dim))
                        x_new_minus = np.zeros((0,f.dim))
                        y_new_plus = np.zeros((0))
                        y_new_minus = np.zeros((0))
                        change_index = []
                        for t in range(1, num_threads + 1):
                            # increase parameter case
                            # check point before running scallop
                            if check_with_one_change(param, parameter_values[param] + (t * step_size[param]), "check") == 1:
                                x_new = check_with_one_change(param, parameter_values[param] + (t * step_size[param]), "")
                                x_new_plus = np.vstack((x_new_plus,x_new))
                                change_index.append(t)
                            # else check==0: something should not be run to avoid error
                            # do nothing
                        if np.shape(x_new_plus)[0] > 0: # at least one combination passed check
                            print('would try plus combination:',x_new_plus)
                            ca_X = np.vstack((ca_X, x_new_plus))
                            y_new_plus = f.compute(x_new_plus, normalize=f.normalize, software_path=path)
                        
                        for t in range(1, num_threads + 1):
                            # decrease parameter case
                            if check_with_one_change(param, parameter_values[param] - (t * step_size[param]), "check") == 1:
                                x_new = check_with_one_change(param, parameter_values[param] - (t * step_size[param]), "")
                                x_new_minus = np.vstack((x_new_minus,x_new))
                                change_index.append((-1)*t)
                            # else check ==0: something should not be run to avoid error
                            # do nothing
                        if np.shape(x_new_minus)[0] > 0: # at least one combination passed check
                            print('would try minus combination:',x_new_minus)
                            ca_X = np.vstack((ca_X, x_new_minus))
                            y_new_minus = f.compute(x_new_minus, normalize=f.normalize, software_path=path)

                        print("Num threads running: " + str(len(x_new_plus)+len(x_new_minus)))
                        
                        #update best auc
                        max_change = 0
                        y_new = np.hstack((y_new_plus, y_new_minus))
                        ca_Y = np.append(ca_Y, y_new)
                        for i in range(len(y_new)):
                            if y_new[i] < cur_auc:
                                cur_auc = y_new[i]
                                max_change = change_index[i]
                                single_param_change = 1
                                made_one_change = 1
                        #update parameter based on best auc
                        parameter_values[param] += max_change * step_size[param]


                    #boolean parameters since only one will ever need to be run.
                    else:
                        # increase cagetory type by 1
                        if check_with_one_change(param,parameter_values[param]+step_size[param],'check')==1:
                            x_new_plus = check_with_one_change(param, parameter_values[param]+step_size[param], '')
                            y_new_plus = float(f.compute(x_new_plus, normalize=f.normalize, software_path=path))
                            ca_X = np.vstack((ca_X, x_new_plus))
                            ca_Y = np.append(ca_Y, y_new_plus)
                            if y_new_plus < cur_auc:
                                cur_auc = y_new_plus
                                parameter_values[param] += step_size[param]
                                single_param_change = 1
                                made_one_change = 1
                        # decrease cagetory type by 1
                        else:
                            if check_with_one_change(param,parameter_values[param]-step_size[param],'check')==1:
                                x_new_minus = check_with_one_change(param, parameter_values[param]-step_size[param], '')
                                y_new_minus = float(f.compute(x_new_minus, normalize=f.normalize, software_path=path))
                                ca_X = np.vstack((ca_X, x_new_minus))
                                ca_Y = np.append(ca_Y, y_new_minus)
                                if y_new_minus < cur_auc:
                                    cur_auc = y_new_minus
                                    parameter_values[param] -= step_size[param]
                                    single_param_change = 1
                                    made_one_change = 1

        # corresponding to last statement in perl script
        break

        # decrease step sizes as long as you can without the steps being too small 
        # (0.01 and 1 for float and int/bool)
        for param in sorted(list(type.keys())): 
            if type[param] == 'int':
                temp = int(step_size[param] * 0.75)
                temp = temp if temp < step_size[param] - 1 else step_size[param] - 1
                if temp > 0:
                    step_size[param] = temp
                    decreased_steps = 1
            if type[param] == 'float':
                temp = format(step_size[param] * 0.75, '.2f')
                temp = float(temp) if float(temp) < step_size[param] - 0.01 else step_size[param] - 0.01
                if temp > 0:
                    step_size[param] = temp
                    decreased_steps = 1

    return (ca_X,ca_Y)
