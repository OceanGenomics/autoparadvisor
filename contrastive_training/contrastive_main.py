from contrastive_model import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--batch_size',type=int,default=32,help='batch size')
parser.add_argument('--lmda',default=0.2,type=float,help='lmda')
parser.add_argument('--beta',default=1.0,type=float,help='beta')
parser.add_argument('--usesubsample',default=1,type=int,help='whether use subsample')
parser.add_argument('--temperature',default=0.05,type=float,help='temperature in contrastive loss')
parser.add_argument('--epochs',default=400,type=int,help='total epochs')
parser.add_argument('--mode',default="out",type=str,help='the mode of loss function')
parser.add_argument('--weightdecay',default=0,type=int,help='whether to use l2 regularization on weights')

args = parser.parse_args()
options = vars(args)
print(options)

'''
mat_params = np.load("mash_mat_params_3_subsample.npy")

encode10_list = ["SRR307903","SRR307911","SRR315323","SRR315334","SRR387661","SRR534291","SRR534307","SRR534319","SRR545695","SRR545723"]
hash_mat = np.zeros((0,1000))
count_mat = np.zeros((0,1000))
for sample_name in encode10_list:
    mash_array = np.load("/mnt/disk18/user/yihangs/learn_and_optimize/mash_files/"+sample_name+".npy")
    hash_mat = np.vstack((hash_mat,mash_array[0]))
    count_mat = np.vstack((count_mat,mash_array[1]))
encode10_mash_normalized_mat = np.concatenate(((hash_mat-mat_params[0,0])/mat_params[0,1],(count_mat-mat_params[1,0])/mat_params[1,1]),axis=1)

hash_mat = np.zeros((0,1000))
count_mat = np.zeros((0,1000))
with open("representative_set_1300_all_hasBO_filtered_new","rt")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        mash_array = np.load("/mnt/disk18/user/yihangs/learn_and_optimize/mash_files/"+line[0]+".npy")
        hash_mat = np.vscontack((hash_mat,mash_array[0]))
        count_mat = np.vstack((count_mat,mash_array[1]))
representative_mash_normalized_mat = np.concatenate(((hash_mat-mat_params[0,0])/mat_params[0,1],(count_mat-mat_params[1,0])/mat_params[1,1]),axis=1)
'''

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#device = torch.device("cpu")
dtype = torch.float32

training_size = 1263

#input feature
#input_mash_mat_np = np.load("/mnt/disk18/user/yihangs/learn_and_optimize/test_sim_GP/input_mash_mat_zscore_3.npy")
#input_mash_mat = torch.from_numpy(input_mash_mat_np).to(dtype=dtype,device=device)
#data_size = input_mash_mat.shape[0]

input_mash_mat_subsample_np = np.load("/contrastive_training/input_features.npy")
input_mash_mat_subsample = torch.from_numpy(input_mash_mat_subsample_np).to(dtype=dtype,device=device)

input_mash_mat_np = input_mash_mat_subsample_np[np.arange(training_size)*8+7]
input_mash_mat = torch.from_numpy(input_mash_mat_np).to(dtype=dtype,device=device)
data_size = input_mash_mat.shape[0]

#similarity matrix
sim_mat_np = np.load("similarity_mat.npy")
#data_size = sim_mat_np.shape[0]
sim_mat_np = (np.ones((data_size,data_size))-np.eye(data_size))*sim_mat_np
sim_mat_np = sim_mat_np/sim_mat_np.max()
sim_mat = torch.from_numpy(sim_mat_np).to(dtype=dtype,device=device)
assert data_size == sim_mat.shape[0]

#similarity matrix with subsampling
sim_mat_subsample_np = np.load("similarity_mat_subsample.npy")
sim_mat_subsample_np = (np.ones((sim_mat_subsample_np.shape[0],sim_mat_subsample_np.shape[0]))-np.eye(sim_mat_subsample_np.shape[0]))*sim_mat_subsample_np
sim_mat_subsample_np = sim_mat_subsample_np/sim_mat_subsample_np.max()
sim_mat_subsample = torch.from_numpy(sim_mat_subsample_np).to(dtype=dtype,device=device)
assert data_size*8 == sim_mat_subsample.shape[0]
assert data_size*8 == input_mash_mat_subsample.shape[0]

use_subsample = args.usesubsample

'''
#Wasserstein similarity matrix
sim_mat_np = np.load("Wasserstein_sim_training_775_topk_20.npy")
data_size = sim_mat_np.shape[0]
sim_mat_np = (np.ones((data_size,data_size))-np.eye(data_size))*sim_mat_np
sim_mat_np = sim_mat_np/sim_mat_np.max()
sim_mat = torch.from_numpy(sim_mat_np).to(dtype=dtype,device=device)
assert data_size == sim_mat.shape[0]
'''

'''
np.random.seed(1118)
validation_set = np.random.choice(training_size,10,replace=False)
print(validation_set)
validation_input = input_mash_mat[validation_set]
print(validation_input.shape)
'''

sample_list = []
with open("../representative_sample_accesion_number","rt")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        sample_list.append(line[0])

assert len(sample_list)==training_size


#set model
lmda = args.lmda
lmda_str = "{:.2f}".format(lmda)
beta = args.beta
beta_str = "{:.2f}".format(beta)
print(lmda,beta)
model = FullyConnectedModel(input_size=2000, output_size=128)
criterion = ContrastiveLoss(temperature=args.temperature,mode=args.mode,lmda=lmda,beta=beta)

if torch.cuda.is_available():
    model = model.cuda()
    criterion = criterion.cuda()

#set optimizer
if args.weightdecay==1:
    #pdb.set_trace()
    optimizer = optim.Adam(model.parameters(),lr=1e-4,weight_decay=1e-5)
else:
    optimizer = optim.Adam(model.parameters(),lr=1e-4)

batch_size = args.batch_size
epochs = args.epochs
loss_list = []

#training loop
#pdb.set_trace()
if use_subsample==1:
    cmd = "mkdir -p ./contrastive_vali_figures_"+str(training_size)+"_lmda_"+lmda_str+"_beta_"+beta_str+"_mode_"+args.mode+"_temp_"+str(args.temperature)+"_wdecay_"+str(args.weightdecay)+"_subsample"
    output_path = "./contrastive_vali_figures_"+str(training_size)+"_lmda_"+lmda_str+"_beta_"+beta_str+"_mode_"+args.mode+"_temp_"+str(args.temperature)+"_wdecay_"+str(args.weightdecay)+"_subsample/"
else:
    cmd = "mkdir -p ./contrastive_vali_figures_"+str(training_size)+"_lmda_"+lmda_str+"_beta_"+beta_str+"_mode_"+args.mode+"_temp_"+str(args.temperature)+"_wdecay_"+str(args.weightdecay)
    output_path = "./contrastive_vali_figures_"+str(training_size)+"_lmda_"+lmda_str+"_beta_"+beta_str+"_mode_"+args.mode+"_temp_"+str(args.temperature)+"_wdecay_"+str(args.weightdecay)+"/"
os.system(cmd)
for epoch in tqdm.tqdm(range(epochs)):
    #pdb.set_trace()
    #generate input pair
    random_pair_np = np.array([np.random.choice(8,2,replace=False) for i in range(data_size)])
    random_pair = torch.from_numpy(random_pair_np).to(dtype=dtype,device=device)
    model.train()
    losses = AverageMeter()
    batch_time = AverageMeter()
    perm_idcs = torch.randperm(data_size)
    batch_num = int(data_size/batch_size)
    if batch_num*batch_size <  data_size:
        batch_num+=1
    #batch_num = 1
    print(batch_num)
    end = time.time()
    
    '''
    if((epoch)%20 ==0):
        #pdb.set_trace()
        with torch.no_grad():
            features_1 = model(torch.from_numpy(encode10_mash_normalized_mat).to(dtype=dtype,device=device))
            features_2 = model(torch.from_numpy(representative_mash_normalized_mat).to(dtype=dtype,device=device))
            features_3 = model(input_mash_mat)
            sim_cal = torch.matmul(features_1, features_2.T)
            sim_cal_2 = torch.matmul(features_1, features_3.T)
            print(sim_cal.shape)
            for i in range(10):
                if(i>3):
                    Y_2 = np.load("/mnt/disk18/user/yihangs/learn_and_optimize/test_sim_GP/"+encode10_list[i]+"/Y_best_representative_set_1300_all_hasBO_filtered_new.npy")
                    fX_test = copula_standardize(deepcopy(Y_2).ravel())
                    plt.figure(figsize=(20, 15))
                    plt.scatter(sim_cal[i],fX_test,marker='o',color="blue")
                else:
                    Y_2 = np.load("/mnt/disk18/user/yihangs/learn_and_optimize/test_sim_GP/"+encode10_list[i]+"/Y_best_representative_set_1300_all_hasBO_filtered_3.npy")
                    fX_test = copula_standardize(deepcopy(Y_2).ravel())
                    plt.figure(figsize=(20, 15))
                    plt.scatter(sim_cal_2[i],fX_test,marker='o',color="blue")
                #plt.show()
                plt.savefig("./contrastive_vali_figures_lmda_"+sys.argv[1]+"_beta_"+sys.argv[2]+"/"+encode10_list[i]+"_"+str(epoch)+".png")
                plt.close()
        with torch.no_grad():
            features_1 = model(validation_input)
            features_2 = model(input_mash_mat)
            sim_cal = torch.matmul(features_1, features_2.T)
            print(sim_cal.shape)
            for i in range(validation_input.shape[0]):
                plt.figure(figsize=(20, 15))
                plt.scatter(sim_cal[i],sim_mat_np[validation_set[i].item()],marker='o',color="blue")
                plt.savefig("./contrastive_vali_figures_lmda_"+sys.argv[1]+"_beta_"+sys.argv[2]+"/"+sample_list[validation_set[i].item()]+"_"+str(epoch)+".png")
                plt.close()
    '''
    for idx in range(batch_num):
        #pdb.set_trace()
        batch_idcs = perm_idcs[idx*batch_size:(idx+1)*batch_size]
        if use_subsample==1:
            random_pair_batch = random_pair[batch_idcs]
            batch_pair_idcs = torch.flatten(random_pair_batch+torch.unsqueeze(batch_idcs,1)*8).to(batch_idcs)
            batch_input = input_mash_mat_subsample[batch_pair_idcs]
            batch_sim_mat = sim_mat_subsample[batch_pair_idcs,:][:,batch_pair_idcs]
        else:
            batch_input = input_mash_mat[batch_idcs]
            batch_sim_mat = sim_mat[batch_idcs,:][:,batch_idcs]
        features = model(batch_input)
        loss = criterion(features,batch_sim_mat)
        losses.update(loss.item(),batch_input.shape[0])
        #pdb.set_trace()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        batch_time.update(time.time()-end)
        end = time.time()
    if((epoch)%20 ==0):
        save_model(model, optimizer, epoch, output_path+"ckpt_epoch_"+str(epoch)+".pth")
    print('epoch {}, average loss {:.4f}'.format(epoch, losses.avg))
    loss_list.append(losses.avg)

save_model(model, optimizer, epochs, output_path+"last.pth")
np.save(output_path+"training_loss.npy",np.array(loss_list))