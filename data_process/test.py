from contrastive_model import *

### load test sample vectors
test_vec = np.load("SRR1583599_normalized.npy")


### load trained models
#checkpoint = torch.load("/mnt/disk18/user/yihangs/learn_and_optimize/test_sim_GP/sra_models/l2_epoch_150.pth")
checkpoint = torch.load("l2_epoch_150.pth")
model = FullyConnectedModel(input_size=2000, output_size=128)
model.load_state_dict(checkpoint['model'])
model.eval()


### compute similarity
train_vecs_subsample_np = np.load("input_features.npy")
train_vecs_np = train_vecs_subsample_np[np.arange(int(train_vecs_subsample_np.shape[0]/8))*8+7]
train_vecs = torch.from_numpy(train_vecs_np).to(dtype=dtype,device=device)

embed_1 = model(torch.from_numpy(test_vec).to(dtype=dtype,device=device))
embed_2 = model(train_vecs)
sim_score = torch.matmul(embed_1, embed_2.T).detach().numpy()

### based on sim_score, choose top-k (e.g. top-30) from training samples
X_bests = np.load("X_bests_scallop.npy")
topk = 30

for i in range(test_vec.shape[0]):
    top_index = np.argsort(-sim_score[i])[:topk]
    recommend_params = X_bests[top_index]