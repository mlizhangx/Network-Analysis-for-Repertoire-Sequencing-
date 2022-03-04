import time
start_time = time.time()

import numpy as np
import pandas as pd
from Levenshtein import distance
from scipy.spatial.distance import pdist, squareform
from scipy import sparse, io


def homogenized(strings):
    lengths = [len(s) for s in strings]
    n = max(lengths)
    for s in strings:
        k = len(s)
        yield [k] + [ord(c) for c in s] + [0] * (n - k)

strings = pd.read_csv('string_info.txt', sep=' ', header=None)
input_string = strings[0].tolist()
# input_string = x

points = np.array(list(homogenized(input_string)))
d = squareform(pdist(points, metric='hamming'))
second_min = np.unique(d)[1]

d = np.where(d <= second_min, 99 , d)
d = np.where(d <= 1, 0 , d)
d = np.where(d == 99 , 1 , d)
select_col_id = (d.sum(axis = 1) > 1)
d = d[select_col_id]
d = d[:,np.array(select_col_id)]

# write distance matrix into file
io.mmwrite('hamming_distance_martrix_degree_lt0.mtx', sparse.csr_matrix(d))
np.savetxt("select_col_id.csv", select_col_id, delimiter=",")

print(len(input_string))
print(d.shape)
print('--- %s seconds ---' % (time.time() - start_time))

