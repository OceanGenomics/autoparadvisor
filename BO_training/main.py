from test_funcs import *
from mixed_test_func import *
from bo.optimizer import Optimizer
from bo.optimizer_mixed import MixedOptimizer
from bo.optimizer_cont import ContOptimizer
import logging
import argparse
import os,pdb
import pickle
import pandas as pd
import time, datetime
from test_funcs.random_seed_config import *
import yaml

# Set up the objective function
parser = argparse.ArgumentParser('Run Experiments')
#NOTE: 'stringTie' used here is to distinguish output file name and command and actually not used as any indicator in codes
parser.add_argument('-p', '--problem', type=str, default='pest', help='current choose can be scallop or stringtie')
parser.add_argument('--max_iters', type=int, default=150, help='Maximum number of BO iterations.')
parser.add_argument('--lamda', type=float, default=1e-6, help='the noise to inject for some problems')
parser.add_argument('--batch_size', type=int, default=1, help='batch size for BO.')
parser.add_argument('--n_trials', type=int, default=20, help='number of trials for the experiment')
parser.add_argument('--n_init', type=int, default=20, help='number of initialising random points')
parser.add_argument('--save_path', type=str, default='output/', help='save directory of the log files')
parser.add_argument('--ard', action='store_true', help='whether to enable automatic relevance determination')
parser.add_argument('--cawarmup', type=int, default=0, help='whether to use coordinate ascent to warm up the process')
parser.add_argument('-a', '--acq', type=str, default='ei', help='choice of the acquisition function.')
parser.add_argument('--random_seed_objective', type=int, default=20, help='The default value of 20 is provided also in COMBO')
parser.add_argument('-d', '--debug', action='store_true', help='Whether to turn on debugging mode (a lot of output will'
                                                               'be generated).')
parser.add_argument('--no_save', action='store_true', help='If activated, do not save the current run into a log folder.')
parser.add_argument('--seed', type=int, default=None, help='**initial** seed setting')
parser.add_argument('-k', '--kernel_type', type=str, default=None, help='specifies the kernel type')
parser.add_argument('--infer_noise_var', action='store_true')
parser.add_argument('--bamid',type=str,default="",help='ID of the bam file')
parser.add_argument('--input_file',type=str, default=None, help='the input file/files of software')
parser.add_argument('--ref_file',type=str,default='',help='reference file for assembler software, gtf format')
parser.add_argument('--software_path',type=str,default='',help='The path of the testing blackbox software')
parser.add_argument('--param_type',type=str,default='mixed',help='parameter type can be category or continuous or mixed')


args = parser.parse_args()
options = vars(args)
print(options)
assert args.max_iters >= args.n_init

if args.debug:
    logging.basicConfig(level=logging.INFO)

# Sanity checks
assert args.acq in ['ucb', 'ei', 'thompson'], 'Unknown acquisition function choice ' + str(args.acq)

#pdb.set_trace()
if not os.path.exists(args.save_path):
    os.makedirs(args.save_path)


for t in range(args.n_trials):

    kwargs = {}
    #pdb.set_trace()
    if args.random_seed_objective is not None:
        assert 1 <= int(args.random_seed_objective) <= 25
        args.random_seed_objective -= 1

    if args.problem == 'pest':
        random_seed_ = sorted(generate_random_seed_pestcontrol())[args.random_seed_objective]
        f = PestControl(random_seed=random_seed_)
        kwargs = {
            'length_max_discrete': 25,
            'length_init_discrete': 20,
        }
    elif args.problem == 'Scallop':
        f = Scallop(input_file=args.input_file, boundary_fold = 0, ref_file=args.ref_file)
        kwargs = {'failtol':18, 'guided_restart':False,'length_init_discrete':100, 'length_min':0.005}
    elif args.problem == 'stringtie':
        f = Scallop(input_file=args.input_file, boundary_fold = 0, ref_file=args.ref_file)
        kwargs = {'failtol':10, 'guided_restart':False,'length_init_discrete':500, 'length_min':0.004}
    elif args.problem == 'func2C':
        f = Func2C(lamda=args.lamda)
    elif args.problem == 'func2C_testDiscrete':
        f = Func2C_testDiscrete(lamda=args.lamda)
        kwargs = {'guided_restart':False}
    elif args.problem == 'func3C':
        f = Func3C(lamda=args.lamda)
        kwargs = {'failtol':10, 'guided_restart':False}
    elif args.problem == 'Func3C_testDiscrete':
        f = Func3C_testDiscrete(lamda=args.lamda)
        kwargs = {'guided_restart':False,'failtol':10}
    elif args.problem == 'ackley53':
        f = Ackley53(lamda=args.lamda)
        kwargs = {
            'length_max_discrete': 50,
            'length_init_discrete': 30,
        }
    elif args.problem == 'MaxSAT60':
        f = MaxSAT60()
        kwargs = {
            'length_max_discrete': 60,
        }
    elif args.problem == 'xgboost-mnist':
        f = XGBoostOptTask(lamda=args.lamda, task='mnist', seed=args.seed)
    else:
        raise ValueError('Unrecognised problem type %s' % args.problem)

    n_categories = f.n_vertices
    problem_type = args.param_type

    print('----- Starting trial %d / %d -----' % ((t + 1), args.n_trials))
    res = pd.DataFrame(np.nan, index=np.arange(int(args.max_iters*args.batch_size)),
                       columns=['Index', 'LastValue', 'BestValue', 'Time'])
    # handle random noise
    if args.infer_noise_var: 
        noise_variance = None
    else: 
        noise_variance = f.lamda if hasattr(f, 'lamda') else None

    #choose correct kernal based on input parameter types
    if args.kernel_type is None:  
        if problem_type == 'mixed':
            kernel_type = 'mixed'
        elif problem_type == 'category':
            kernel_type = 'transformed_overlap'
        elif problem_type == 'continuous':
            kernel_type = 'continuous'
        else:
            raise ValueError('cannot define kernal type without specify parameter type')
    #if kernal type is specified at the beginning
    else: 
        kernel_type = args.kernel_type

    # Perl re-implement
    if args.cawarmup > 0:
        #x_next,y_next = coordinate_ascent_warmup_yaml(f, max_iters=args.cawarmup, 
        #                        path=args.software_path, num_threads = 1)
        ca = CoordinateAscent(f, max_iters=args.cawarmup, num_threads=1, path=args.software_path)
        x_next, y_next = ca.coordinate_ascent_warmup_yaml()
        print(x_next, y_next)
        # search space is not ajust on cag parameters 
        x_next_min = x_next[np.array(y_next).argmin()]
        x_next_cont = np.array([])
        for idx, val in enumerate(f.continuous_dims):
            x_next_cont = np.append(x_next_cont, x_next_min[val])
        f.ub = np.maximum(x_next_cont*2,f.ub)
        #if upperbound was changed by cawarmup make sure is below hard ub
        f.ub = np.minimum(f.ub, f.hard_ub)


    if problem_type == 'mixed':
        optim = MixedOptimizer(f.config, f.lb, f.ub, f.continuous_dims, f.categorical_dims, int_constrained_dims=f.int_constrained_dims,
                               default_x=f.default, n_init=args.n_init, use_ard=args.ard, acq=args.acq,
                               kernel_type=kernel_type,
                               noise_variance=noise_variance,
                               **kwargs)
    # add only continuous param case
    elif problem_type == 'continuous':
        optim = ContOptimizer(f.config, f.lb, f.ub, f.continuous_dims, int_constrained_dims=f.int_constrained_dims,
                               default_x=f.default, n_init=args.n_init, use_ard=args.ard, acq=args.acq,
                               kernel_type=kernel_type,
                               noise_variance=noise_variance,
                               **kwargs)
    elif problem_type == 'category':
        optim = Optimizer(f.config, default_x=f.default, n_init=args.n_init, use_ard=args.ard, acq=args.acq,
                          kernel_type=kernel_type,
                          noise_variance=noise_variance, **kwargs)
    else:
        raise ValueError('Unsupported parameter type %s' % args.param_type)

    T_array = []
    start = time.time()
    #pdb.set_trace()
    if args.cawarmup >0:
        _ = optim.suggest(args.batch_size)
        #x_next,y_next = coordinate_ascent_warmup(f,fold_step=5,iterations=args.n_init)
        optim.observe(x_next, y_next)
        n_init_num = len(y_next)
    else:
        x_next = optim.suggest(args.batch_size)
        y_next = f.compute(x_next, normalize=f.normalize,software_path=args.software_path)
        optim.observe(x_next, y_next)
        n_init_num = args.n_init
    end = time.time()
    T_array.append(end - start)
    for i in range(args.max_iters-n_init_num):
        #pdb.set_trace()
        start = time.time()
        x_next = optim.suggest(args.batch_size)
        y_next = f.compute(x_next, normalize=f.normalize,software_path=args.software_path)
        optim.observe(x_next, y_next)
        end = time.time()
        T_array.append(end - start)
        if f.normalize:
            Y = np.array(optim.casmopolitan.fX) * f.std + f.mean
        else:
            Y = np.array(optim.casmopolitan.fX)
        if optim.casmopolitan.length <= optim.casmopolitan.length_min or optim.casmopolitan.length_discrete <= optim.casmopolitan.length_min_discrete:
            #optim.restart()
            break
        '''
        if Y[:i].shape[0]:
            # sequential
            if args.batch_size == 1:
                res.iloc[i, :] = [i, float(Y[-1]), float(np.min(Y[:i])), end-start]
            # batch
            else:
                for idx, j in enumerate(range(i*args.batch_size, (i+1)*args.batch_size)):
                    res.iloc[j, :] = [j, float(Y[-idx]), float(np.min(Y[:i*args.batch_size])), end-start]
            # x_next = x_next.astype(int)
            argmin = np.argmin(Y[:i*args.batch_size])

            print('Iter %d, Last X %s; fX:  %.4f. X_best: %s, fX_best: %.4f'
                  % (i, x_next.flatten(),
                     float(Y[-1]),
                     ''.join([str(int(i)) for i in optim.casmopolitan.X[:i * args.batch_size][argmin].flatten()]),
                     Y[:i*args.batch_size][argmin]))
        '''
        if Y[:i+args.n_init].shape[0]:
            argmin = np.argmin(Y)
            print('Iter %d, Last X %s; fX:  %.4f. X_best: %s, fX_best: %.4f'
                  % (i, x_next.flatten(),
                     float(Y[-1]),
                     ''.join([str(int(i)) for i in optim.casmopolitan.X[argmin].flatten()]),
                     Y[argmin]))
    np.save(args.save_path + "/" + args.problem + "_wall_clock_"+str(t+1)+".npy",np.array(T_array))
    np.save(args.save_path + "/" + args.problem + "_X_"+str(t+1)+".npy",optim.casmopolitan.X)
    np.save(args.save_path + "/" + args.problem + "_Y_"+str(t+1)+".npy",Y)
    print('process done!')

    if args.seed is not None:
        args.seed += 1
