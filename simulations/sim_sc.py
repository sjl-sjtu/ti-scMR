#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from numpy import random
from prosstt import tree
from prosstt import simulation as sim
from prosstt import sim_utils as sut
import scipy as sp
rseed = 100 
random.seed(rseed)

import anndata as ad
# from scanpy import pp
# from scanpy.tl import diffmap

import pandas as pd

import os


# In[2]:


N = 500
T = 10
P = 10000
top = [["A", "B"]]
branches = np.unique(np.array(top).flatten())
time = {b: T for b in branches}
G = 100
lineage = tree.Tree(topology=top, G=G, time=time, num_branches=len(branches), branch_points=0, modules=8)

# prior hyperparameter for h
a = min(1/lineage.modules,0.05)

# assign genes to programs to calculate relative average gene expression
# uMs, Ws, Hs = sim.simulate_lineage(lineage, a=a)
uMs, Ws, Hs = sim.simulate_lineage(lineage, a=a)
#uMs: rel_means (T*G); Ws: programs (T*modules); Hs:modules*G


# In[3]:


dat = pd.read_csv("geno.csv")
dat


# In[4]:


# simulate eqtl effect
ntr = list(np.random.randint(15,25,G))
def getlambda(i):
    tr = np.random.choice(np.arange(1,P+1),ntr[i-1],replace=False)  #np.arange(ntr*(i-1)+1, ntr*i+1)
    beta = np.random.gamma(shape=3.6, scale=1/12, size=ntr[i-1])
    lambda_ = np.dot(dat.iloc[:, tr-1], beta)
    return tr,lambda_

result_list = list(map(getlambda, range(1, G+1)))
tr = [result[0] for result in result_list]
lambda_ = np.array([result[1] for result in result_list]) #G*N


# In[5]:


# true eqtls
arr = np.zeros((G,P), dtype=int)
for i, sublist in enumerate(tr):
    arr[i, sublist-1] = 1
dfiv = pd.DataFrame(np.transpose(arr))
dfiv.columns = ["gene" + str(g) for g in range(1, G+1)]
dfiv['snp'] = ["snp" + str(p) for p in range(1, P+1)]
dfiv


# In[6]:


# scale relative average expression to absolute values
# base expression
gene_scale = sut.simulate_base_gene_exp(lineage, uMs)  #G*1


# time-invariant SNP-gene

# In[7]:


os.chdir("stable_eqtl")


# In[8]:


lineage_list = []
dfmu_list = []


# In[9]:


import copy
for i in range(N):
    lineage1 = copy.deepcopy(lineage)
#     gene_scale1 = lambda_[:,i]*gene_scale+gene_scale
#     lineage1.add_genes(uMs, gene_scale1)
    lineage1.add_genes(uMs, gene_scale)
    for j in list(lineage1.means.keys()):
        ts = np.arange(lineage1.branch_times()[j][0],lineage1.branch_times()[j][1]+1)
        # lambda_1 = np.outer(np.sqrt(2)*np.abs(np.sin(ts)),lambda_[:,i]) #time-varying
        lambda_1 = np.outer(ts*0+1,lambda_[:,i]) #time-invariant
        lineage1.means[j] = lineage1.means[j]+lineage1.means[j]*lambda_1
    lineage_list.append(lineage1)
    dfmu1 = pd.DataFrame(np.concatenate(list(lineage1.means.values())))
    dfmu1.columns = ["gene" + str(g) for g in range(1, G+1)]
    dfmu1["pseudotime"] = np.array([np.arange(i,j+1) for i,j in lineage1.branch_times().values()]).flatten()
    dfmu1['state'] = np.array([[i]*j for i,j in lineage1.time.items()]).flatten()
    dfmu1['id'] = i+1
    dfmu_list.append(dfmu1)


# In[10]:


dfmu = pd.concat(dfmu_list)


# In[11]:


dfmu


# In[12]:


dfmu.to_csv("true_sc.csv",index=0)
dfiv.to_csv("dfiv.csv",index=0)


# In[13]:


# sampling varaince hyperparameters
alpha = np.exp(random.normal(loc=np.log(0.2), scale=np.log(1.5), size=lineage.G))
beta = np.exp(random.normal(loc=np.log(1), scale=np.log(1.5), size=lineage.G)) + 1


# In[14]:


df_list = []
for i in range(N):
    lineage1 = lineage_list[i]
    
    n = np.random.randint(8,12)
    X, labs, brns, scalings = sim.sample_whole_tree(lineage1, n, alpha=alpha, beta=beta)
    
    df1 = pd.DataFrame(X)
    df1.columns = list(map(lambda x: "gene"+str(x),range(1,G+1)))
    df1['pseudotime'] = labs
    df1['state'] = brns
    df1 = df1.sort_values(["pseudotime"])
    df1['id'] = i+1
    df_list.append(df1)


# In[15]:


df = pd.concat(df_list)


# In[16]:


df['cellid'] = df['id'].astype(str) + '_' + (df.groupby('id').cumcount()).astype(str)


# In[17]:


df


# In[18]:


df.to_csv("complete_sample/sim_sc.csv",index=0)


# In[19]:


df_list = []
for i in range(N):
    lineage1 = lineage_list[i]
    
    maxT = lineage1.get_max_time()
    n = np.random.randint(10*maxT, 15*maxT)

    series_points = list(np.random.choice(np.arange(0,maxT),5,replace=False))
    point_std = list(np.random.uniform(0.5,2,5))
    X, labs, brns, scalings = sim.sample_pseudotime_series(lineage1, cells=n,
                                                           series_points=series_points,
                                                           point_std=point_std,
                                                           alpha=alpha, beta=beta)

    # X = (X.transpose() / scalings).transpose()
    df1 = pd.DataFrame(X)
    df1.columns = list(map(lambda x: "gene"+str(x),range(1,G+1)))
    df1['pseudotime'] = labs
    df1['state'] = brns
    df1 = df1.sort_values(["pseudotime"])
    df1['id'] = i+1
    df_list.append(df1)


# In[20]:


df = pd.concat(df_list)
df['cellid'] = df['id'].astype(str) + '_' + (df.groupby('id').cumcount()).astype(str)
df


# In[21]:


df.to_csv("mixed_sample/sim_sc.csv",index=0)


# time-varying SNP-gene

# In[22]:


os.chdir("../dyn_eqtl")


# In[23]:


lineage_list = []
dfmu_list = []
import copy
for i in range(N):
    lineage1 = copy.deepcopy(lineage)
#     gene_scale1 = lambda_[:,i]*gene_scale+gene_scale
#     lineage1.add_genes(uMs, gene_scale1)
    lineage1.add_genes(uMs, gene_scale)
    for j in list(lineage1.means.keys()):
        ts = np.arange(lineage1.branch_times()[j][0],lineage1.branch_times()[j][1]+1)
        lambda_1 = np.outer(np.sqrt(2)*np.abs(np.sin(ts)),lambda_[:,i]) #time-varying
        # lambda_1 = np.outer(ts*0+1,lambda_[:,i]) #time-invariant
        lineage1.means[j] = lineage1.means[j]+lineage1.means[j]*lambda_1
    lineage_list.append(lineage1)
    dfmu1 = pd.DataFrame(np.concatenate(list(lineage1.means.values())))
    dfmu1.columns = ["gene" + str(g) for g in range(1, G+1)]
    dfmu1["pseudotime"] = np.array([np.arange(i,j+1) for i,j in lineage1.branch_times().values()]).flatten()
    dfmu1['state'] = np.array([[i]*j for i,j in lineage1.time.items()]).flatten()
    dfmu1['id'] = i+1
    dfmu_list.append(dfmu1)
    
dfmu = pd.concat(dfmu_list)


# In[24]:


dfmu.to_csv("true_sc.csv",index=0)
dfiv.to_csv("dfiv.csv",index=0)


# In[25]:


# sampling varaince hyperparameters
alpha = np.exp(random.normal(loc=np.log(0.2), scale=np.log(1.5), size=lineage.G))
beta = np.exp(random.normal(loc=np.log(1), scale=np.log(1.5), size=lineage.G)) + 1


# In[26]:


df_list = []
for i in range(N):
    lineage1 = lineage_list[i]
    
    n = np.random.randint(8,12)
    X, labs, brns, scalings = sim.sample_whole_tree(lineage1, n, alpha=alpha, beta=beta)
    
    df1 = pd.DataFrame(X)
    df1.columns = list(map(lambda x: "gene"+str(x),range(1,G+1)))
    df1['pseudotime'] = labs
    df1['state'] = brns
    df1 = df1.sort_values(["pseudotime"])
    df1['id'] = i+1
    df_list.append(df1)

df = pd.concat(df_list)
df['cellid'] = df['id'].astype(str) + '_' + (df.groupby('id').cumcount()).astype(str)
df.to_csv("complete_sample/sim_sc.csv",index=0)


# In[27]:


df_list = []
for i in range(N):
    lineage1 = lineage_list[i]
    
    maxT = lineage1.get_max_time()
    n = np.random.randint(10*maxT, 15*maxT)

    series_points = list(np.random.choice(np.arange(0,maxT),5,replace=False))
    point_std = list(np.random.uniform(0.5,2,5))
    X, labs, brns, scalings = sim.sample_pseudotime_series(lineage1, cells=n,
                                                           series_points=series_points,
                                                           point_std=point_std,
                                                           alpha=alpha, beta=beta)

    # X = (X.transpose() / scalings).transpose()
    df1 = pd.DataFrame(X)
    df1.columns = list(map(lambda x: "gene"+str(x),range(1,G+1)))
    df1['pseudotime'] = labs
    df1['state'] = brns
    df1 = df1.sort_values(["pseudotime"])
    df1['id'] = i+1
    df_list.append(df1)


# In[28]:


df = pd.concat(df_list)
df['cellid'] = df['id'].astype(str) + '_' + (df.groupby('id').cumcount()).astype(str)
df.to_csv("mixed_sample/sim_sc.csv",index=0)


# In[ ]:




