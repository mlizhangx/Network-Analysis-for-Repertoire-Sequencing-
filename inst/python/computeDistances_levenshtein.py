import numpy as np
import pandas as pd
from Levenshtein import distance
from scipy.spatial.distance import pdist, squareform
from scipy import sparse, io


strings = pd.read_csv('temp_clone_seq_list.txt', sep=' ', header=None)
input_string = strings[0].tolist()
transformed_strings = np.array(input_string).reshape(-1,1)
distance_matrix = pdist(transformed_strings,lambda x,y: distance(x[0],y[0]))
d = squareform(distance_matrix)

d = np.where(d > 1, 0 , d)
select_col_id = (d.sum(axis = 1) > 0)
d = d[select_col_id]
d = d[:,np.array(select_col_id)]


# write distance matrix into file
io.mmwrite('temp_adjacency_matrix.mtx', sparse.csr_matrix(d))
np.savetxt("select_col_id.csv", select_col_id, delimiter=",")


