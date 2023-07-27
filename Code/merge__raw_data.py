# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 15:37:14 2023

@author: Benson
"""

import scanpy as sc
import pandas as pd
import numpy as np

import anndata # merge anndata



def read_data(filename):
    adata = sc.read_mtx(f"C:/Users/User/Desktop/論文/PBMC_Data/GSE165080/{filename}/matrix.mtx")
    adata = adata.T
    
    gene = pd.read_table(f"C:/Users/User/Desktop/論文/PBMC_Data/GSE165080/{filename}/features.tsv", header = None)
    gene = np.array(gene[1])
    
    barcode = pd.read_table(f"C:/Users/User/Desktop/論文/PBMC_Data/GSE165080/{filename}/barcodes.tsv", header = None)
    barcode = np.array(barcode[0])
    
    adata.obs_names = barcode
    adata.var_names = gene
    adata.var_names_make_unique()
    
    return adata
    
# Read Dataset information
data = pd.read_csv("SraRunTable.txt", delimiter=",")
data = data[["Run", "Age","gender"]]
#data.loc[data['Run'] == 'SRR13482553', 'Age'].values[0]
Dataset = data["Run"]
    
    
# 把每個Raw data合併到 H5ad file 中
first = True
for i in Dataset:
    print(f"Start {i}！\n")
    
    if first == False: # not the first round
        # read raw data in h5ad
        Raw = sc.read_h5ad("Raw.h5ad")
    
    # read file
    adata = read_data(i)
    
    # Add age, gender and dataset
    adata.obs["Age"] = [data.loc[data['Run'] == i, 'Age'].values[0]]*adata.n_obs
    adata.obs["Gender"] = [data.loc[data['Run'] == i, 'gender'].values[0]]*adata.n_obs
    adata.obs["Dataset"] = [i]*adata.n_obs
    
    if first == True: # not the first round
        first = False
        adata.write_h5ad("Raw.h5ad")
    else: # after first round
        # merge anndata
        merged_adata = anndata.concat([Raw, adata], axis=0)
        # write file
        merged_adata.write_h5ad("Raw.h5ad")  
    
    print(f"Finished {i} merge！\n")