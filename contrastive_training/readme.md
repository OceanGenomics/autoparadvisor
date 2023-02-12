The dataset used for the contrastive training can be downloaded from https://kilthub.cmu.edu/articles/dataset/contrastive_training_data_tar_gz/22010888. It has three files:


* `input_features.npy`: The 2000 dimensional MinHash sketch features for each representative sample and its subsampled variants. So the total dimension is 10104*2000. (10104=1263 representative samples * 8 variants)
* `similarity_mat.npy`: A 1263*1263 similarity matrix containing similarity values between every pair of representative samples. 
* `similarity_mat_subsample.npy`: A 10104*10104 similarity matrix containing similarity values between every pair of representative samples including their subsampled variants.

To train the model, put these three files and the two python scripts `contrastive_main.py` and `contrastive_model.py`, and run the following command: 

`python contrastive_main.py --batch_size 632 --epochs 150 --mode l2`