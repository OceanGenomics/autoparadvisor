### python script to generate .npy vector from .msh file


import glob,os
import numpy as np
import capnp

capnp.remove_import_hook()
#you may need to change the path of MinHash.capnp file
mash_capnp = capnp.load("MinHash.capnp")

msh_file = "SRR1583599.msh"
file_name = "SRR1583599"

f_1 = open(msh_file,"rb")
read_mash = mash_capnp.MinHash.read(f_1)
hashes = list(read_mash.referenceListOld.references[0].hashes64)
counts = list(read_mash.referenceListOld.references[0].counts32)
assert len(hashes)==len(counts)
file_npy = np.array([hashes,counts])
np.save(file_name+".npy",file_npy)
f_1.close()