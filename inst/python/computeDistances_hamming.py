import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy import sparse, io


def homogenized(strings):
    lengths = [len(s) for s in strings]
    n = max(lengths)
    for s in strings:
        k = len(s)
        yield [k] + [ord(c) for c in s] + [0] * (n - k)

strings = pd.read_csv('temp_clone_seq_list.txt', sep=' ', header=None)
input_string = strings[0].tolist()
# input_string = x

points = np.array(list(homogenized(input_string)))
d = squareform(pdist(points, metric = 'hamming'))
second_min = np.unique(d)[1]

d = np.where(d <= second_min, 99 , d)
d = np.where(d <= 1, 0 , d)
d = np.where(d == 99 , 1 , d)
select_col_id = (d.sum(axis = 1) > 1)
d = d[select_col_id]
d = d[:,np.array(select_col_id)]

# write distance matrix into file
io.mmwrite('temp_adjacency_matrix.mtx', sparse.csr_matrix(d))
np.savetxt("select_col_id.csv", select_col_id, delimiter=",")
