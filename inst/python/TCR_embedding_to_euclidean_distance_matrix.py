#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
# from Levenshtein import distance
from scipy.spatial.distance import pdist, squareform
from scipy import sparse, io


# In[2]:


TCR_Atchley_factor = pd.read_csv('TCR_w_Atchley_factor_tmp_file.csv',index_col = "contig_id")


# In[4]:


matrix_arry = np.array(TCR_Atchley_factor)


# In[5]:


d = squareform(pdist(matrix_arry, metric='euclidean'))


# In[7]:


d.shape


# In[8]:


d_mtx = sparse.csr_matrix(d)


# In[12]:


# write distance matrix into file
io.mmwrite('euclidean_distance_martrix_Atchley_factor_tmp_file.mtx', d_mtx)


# In[ ]:




