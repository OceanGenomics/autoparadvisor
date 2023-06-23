import numpy as np

### z-score normalization, (x-mu)/sigma


# this normalization file should be changed if the training data is changed
normalization_params = np.load("mash_mat_params_4_subsample.npy")

hash_mat = np.zeros((0,1000))
count_mat = np.zeros((0,1000))

mash_npy = np.load("SRR1583599.npy")

hash_mat = np.vstack((hash_mat,mash_npy[0]))
count_mat = np.vstack((count_mat,mash_npy[1]))
normalized_vec = np.concatenate(((hash_mat-normalization_params[0,0])/normalization_params[0,1],(count_mat-normalization_params[1,0])/normalization_params[1,1]),axis=1)
np.save("SRR1583599_normalized.npy",normalized_vec)