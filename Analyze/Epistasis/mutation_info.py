import numpy as np

num_muts_total = 15
mutations = np.array([str(x) for x in range(1,num_muts_total+1)])
mut_names = np.array(['339','371','373','375','417','440','446','477','478','484','493','496','498','501','505'])

num_muts_H1 = 15
H1_indices = np.arange(num_muts_total,dtype=int)
H1_mutations = mutations
H1_mut_names = mut_names
