"""Broadly useful python packages"""
import pandas as pd
import os
import numpy as np
import pickle
from copy import deepcopy
from shutil import move
import warnings
from metric_compute import cluster_acc, cluster_purity, evaluate
import os
import csv
from time import time as get_time
import tracemalloc
"""Machine learning and single cell packages"""
import sklearn.metrics as metrics
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi
import scanpy as sc
from anndata import AnnData

"""CarDEC Package"""
from CarDEC import CarDEC_API 

import sys
from pathlib import Path
import time

dataset_name = sys.argv[1] 
n_cluster = int(sys.argv[2])
Benchamark_Result_csv = sys.argv[3] 
benchmark_randomseed = int(sys.argv[4])
protein_num_value = int(sys.argv[5])
Benchamark_data_path = sys.argv[6]

x = pd.read_csv(Benchamark_data_path+'/ADTCount.csv', header=0, index_col=0).T.values
cell_name = pd.read_csv(Benchamark_data_path+'/celltypes.csv', usecols=['celltype']).values
_, y = np.unique(cell_name, return_inverse=True)

adata = sc.AnnData(x,dtype=np.float32)

current_time = time.time()
time_str = time.strftime('%Y-%m-%d_%H:%M:%S', time.localtime(current_time))
save_dir = 'weights_dir/'+dataset_name+'_'+time_str

tracemalloc.start()
start_time = get_time()

CarDEC = CarDEC_API(adata, weights_dir = save_dir, n_high_var = protein_num_value, LVG = True)
CarDEC.build_model(n_clusters = n_cluster, random_seed = benchmark_randomseed)
CarDEC.make_inference()
CarDEC.model_counts()
q = deepcopy(CarDEC.dataset.obsm['cluster memberships']) #The cluster membership numpy array

end_time = get_time()
Running_time = end_time - start_time
memory_current, memory_peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
memory_peak = memory_peak/(1024*1024)

cardec_labels = np.argmax(q, axis=1)
true_label = np.array(y)
pred_label = np.array(cardec_labels)

acc, purity, nmi, ari= evaluate(true_label,pred_label)
type = "ADT"
method_name = "CarDEC-ADT"
df_pred = pd.DataFrame(pred_label, columns=['pred_label'])
directory_path = os.path.dirname(Benchamark_Result_csv)
pred_file_path = os.path.join(directory_path, "{}_pred.csv".format(method_name))
df_pred.to_csv(pred_file_path, index=False, header=True)
with open(Benchamark_Result_csv, 'a+', newline="") as f:
    csv_write = csv.writer(f)
    data_row = [dataset_name, type, method_name, ari, nmi, acc, purity, memory_peak, Running_time]
    csv_write.writerow(data_row)
