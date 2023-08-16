import numpy as np
from test_funcs.base import TestFunction
import pdb,sys,os,csv,glob
from multiprocessing import Manager
import multiprocessing
import subprocess
import yaml
import copy

def Scallop_base(index,x,Result,input_file,software_path,docs):
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
            label_list = parameter[parameter_name]['label']
            cag_label = label_list[int(x[i])]
            if parameter[parameter_name]['usage']=='TF':
                pcmd = ' '.join([pcmd, parameter[parameter_name]['prefix'], cag_label])
            elif parameter[parameter_name]['usage']=='turn_on':
                pcmd = ' '.join([pcmd, cag_label])
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
        n_vertices = []
        
        for i in range(len(parameter_bounds)):
            #single parameter dict that stores all information for one param
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
                n_vertices.append(len(parameter[parameter_name]['label']))
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
        self.int_constrained_dims = np.asarray(int_constrained_dims) if len(int_constrained_dims)!=0 else None
        self.hard_lb = np.array(hard_lb)
        self.hard_ub = np.array(hard_ub)
        self.n_vertices = np.array(n_vertices)
        
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
        
        #TODO: now support multi cag parameters
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
        if self.int_constrained_dims is not None:
            X[:,self.int_constrained_dims] = np.round(X[:,self.int_constrained_dims])
        with Manager() as manager:
            Y = manager.list(range(N))
            process_list = []
            for i in range(N):
                tmp_process = multiprocessing.Process(target=Scallop_base, \
                    args=(i,X[i],Y,self.input_file,software_path,self.docs))
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


#-------------Coordinate Ascent class---------------
class CoordinateAscent():
    def __init__(self, f, max_iters=60, num_threads=1, path=''):
        '''
        acheived same function of Perl script
        f: blackbox function, e.g Scallop()
        max_iters: max iteration allowed from ca_warmup
        num_threads: controls the number of parallel steps to take in each direction
        path: system path to blackbox software
        '''
        self.f = f 
        self.max_iters = max_iters
        self.num_threads = num_threads
        self.path = path

        # get same variables as perl script
        self.parameter_values = copy.deepcopy(f.parameter_values)
        self.step_size = copy.deepcopy(f.step_size)
        self.type = copy.deepcopy(f.para_type)

    def check_with_one_change(self, param_to_change, param_value):
        # get param position in unsorted parameter list
        index = self.f.para_to_index[param_to_change]
        index_within_ub_lb = np.where(self.f.continuous_dims==index)
        # change param and check if within bounds
        # cag param should between 0-1
        if param_value!='' and self.type[param_to_change]=='cag': 
            if param_value>1 or param_value<0:
                return 0
        # int and float param should between hard lower bound and hard upper bound
        elif param_value!='' and self.type[param_to_change]!='cag':
            if float(param_value) > self.f.hard_ub[index_within_ub_lb] or \
                float(param_value) < self.f.hard_lb[index_within_ub_lb]:
                return 0
        # if got here and in check mode, means check is passed
        return 1

    def run_with_one_change(self, param_to_change, param_value):
        # ensures the floating point parameters don't get too complex
        if self.type[param_to_change] == 'float':
            param_value = float(format(param_value, '.2f'))
        # use index position to find the 'parameter to change'
        x_new = [param_value if p==param_to_change else self.parameter_values[p] for p in list(self.type.keys())]
        return x_new

    def compute_nonboolen_changes(self, param, change_index, direction, ca_X):
        x_new_change = np.zeros((0, self.f.dim))
        y_new_change = np.zeros((0))
        for t in range(1, self.num_threads + 1):
            # check point before running scallop
            parameter_value_new = self.parameter_values[param] + direction * (t * self.step_size[param])
            if self.check_with_one_change(param, parameter_value_new) == 1:
                x_new = self.run_with_one_change(param, parameter_value_new)
                x_new_change = np.vstack((x_new_change,x_new))
                change_index.append(direction * t)
                # else check==0: something should not be run to avoid error
                # do nothing
        if np.shape(x_new_change)[0] > 0: # at least one combination passed check
            print('would try combination:',x_new_change)
            ca_X = np.vstack((ca_X, x_new_change))
            y_new_change = self.f.compute(x_new_change, normalize=self.f.normalize, \
                                            software_path=self.path)

        return x_new_change, y_new_change, ca_X

    def update_one_param(self, param, cur_auc, ca_X, ca_Y):
        # non-boolean parameters: run values at parallel step size 
        if self.type[param] != "bool":
            change_index = []

            # increase with direction 1 and decrease with direction -1
            x_new_plus, y_new_plus, ca_X = self.compute_nonboolen_changes(param, \
                                            change_index, 1, ca_X)
            x_new_minus, y_new_minus, ca_X = self.compute_nonboolen_changes(param, \
                                            change_index, -1, ca_X)
            y_new = np.hstack((y_new_plus, y_new_minus))
            ca_Y = np.append(ca_Y, y_new)
            print("Num threads running: " + str(len(x_new_plus)+len(x_new_minus)))
            
            #update best auc
            max_change = 0
            for i in range(len(y_new)):
                if y_new[i] < cur_auc:
                    cur_auc = y_new[i]
                    max_change = change_index[i]
                    self.single_param_change = 1
                    self.made_one_change = 1
            #update parameter based on best auc
            self.parameter_values[param] += max_change * self.step_size[param]

        #boolean parameters since only one will ever need to be run.
        else:
            # increase cagetory type by 1
            parameter_value_new = self.parameter_values[param] + self.step_size[param]
            if self.check_with_one_change(param, parameter_value_new)==1:
                x_new_plus = self.run_with_one_change(param, parameter_value_new)
                y_new_plus = float(self.f.compute(x_new_plus, normalize=self.f.normalize, \
                                                software_path=self.path))
                ca_X = np.vstack((ca_X, x_new_plus))
                ca_Y = np.append(ca_Y, y_new_plus)
                # update best auc
                if y_new_plus < cur_auc:
                    cur_auc = y_new_plus
                    self.parameter_values[param] += self.step_size[param]
                    self.single_param_change = 1
                    self.made_one_change = 1

            # decrease cagetory type by 1
            else:
                parameter_value_new = self.parameter_values[param] - self.step_size[param]
                if self.check_with_one_change(param, parameter_value_new)==1:
                    x_new_minus = self.run_with_one_change(param, parameter_value_new)
                    y_new_minus = float(self.f.compute(x_new_minus, normalize=self.f.normalize, \
                                                        software_path=self.path))
                    ca_X = np.vstack((ca_X, x_new_minus))
                    ca_Y = np.append(ca_Y, y_new_minus)
                    # update best auc
                    if y_new_minus < cur_auc:
                        cur_auc = y_new_minus
                        self.parameter_values[param] -= self.step_size[param]
                        self.single_param_change = 1
                        self.made_one_change = 1

        return cur_auc, ca_X, ca_Y

    def run_param_change(self, param, cur_auc, ca_X, ca_Y):
        ''''''
        # loops as long as you can continue moving on this parameter (coordinate) and still increase AUC
        self.single_param_change = 1
        while self.single_param_change == 1:
            self.single_param_change = 0
            iter_num = len(ca_Y)
            print('iterarion number: ',iter_num)
            # if exceeds max iteration number return
            if self.max_iters > 0 and iter_num > self.max_iters: 
                self.made_one_change = 0
                return cur_auc, ca_X, ca_Y

            # compute auc based on one param change
            sys.stderr.write(f"Updating {param}, type: {self.type[param]} \n")
            cur_auc, ca_X, ca_Y = self.update_one_param(param, cur_auc, ca_X, ca_Y)
            
        return cur_auc, ca_X, ca_Y


    # main function
    def coordinate_ascent_warmup_yaml(self):
        # initialize and gets default AUC
        print('get default auc in ca_warmup')
        ca_X = np.zeros((0,self.f.dim))
        ca_X = np.vstack((ca_X,self.f.default))
        cur_auc = float(self.f.compute(self.f.default, normalize=self.f.normalize, \
                                        software_path=self.path))
        ca_Y = []
        ca_Y.append(cur_auc)

        # full loop starts here
        print('starting ca_warmup-----')
        # loops as long as you were able to decrease one of the step sizes
        self.decreased_steps = 1
        while self.decreased_steps == 1:
            self.decreased_steps = 0

            # loops as long as you have made a change in the parameter vector
            # without a change in step size
            self.made_one_change = 1
            while self.made_one_change == 1:
                self.made_one_change = 0

                # loops as long as you can continue moving on this parameter (coordinate) and still increase AUC
                for param in sorted(list(self.type.keys())): 
                    print("best current auc:", cur_auc)
                    # try to change every parameter once and accmulate changes
                    cur_auc, ca_X, ca_Y = self.run_param_change(param, cur_auc, ca_X, ca_Y)

            # statement here only to match Perl function
            break

            # decrease step sizes as long as you can without the steps being too small 
            # dead code, not used when break is applied 
            for param in sorted(list(type.keys())): 
                if self.type[param] == 'int':
                    temp = int(self.step_size[param] * 0.75)
                    temp = temp if temp < self.step_size[param] - 1 else self.step_size[param] - 1
                    if temp > 0:
                        self.step_size[param] = temp
                        self.decreased_steps = 1
                if type[param] == 'float':
                    temp = format(self.step_size[param] * 0.75, '.2f')
                    temp = float(temp) if float(temp) < self.step_size[param] - 0.01 else self.step_size[param] - 0.01
                    if temp > 0:
                        self.step_size[param] = temp
                        self.decreased_steps = 1

        return (ca_X,ca_Y)