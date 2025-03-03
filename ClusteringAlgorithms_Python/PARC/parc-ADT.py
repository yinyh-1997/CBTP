import scanpy as sc
import numpy as np
import parc
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import datasets
from metric_compute import cluster_acc, cluster_purity, evaluate
import os
import csv
from time import time as get_time
import tracemalloc
import sys
from pathlib import Path

dataset_name = sys.argv[1] 
n_cluster = int(sys.argv[2])
Benchamark_Result_csv = sys.argv[3] 
benchmark_randomseed = int(sys.argv[4])
Benchamark_data_path = sys.argv[5]

x = pd.read_csv(Benchamark_data_path+'/ADT.csv', header=0, index_col=0).T.values
cell_name = pd.read_csv(Benchamark_data_path+'/celltypes.csv', usecols=['celltype']).values
_, y = np.unique(cell_name, return_inverse=True)

adata = sc.AnnData(x,dtype=np.float32)
# sc.pp.highly_variable_genes(adata,n_top_genes=2000,subset=True,inplace=True)

tracemalloc.start()
start_time = get_time()

p = parc.PARC(adata.X,random_seed=benchmark_randomseed)
p.run_PARC()

end_time = get_time()
Running_time = end_time - start_time
memory_current, memory_peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
memory_peak = memory_peak/(1024*1024)

parc_labels = p.labels
true_label = np.array(y)
pred_label = np.array(parc_labels)

acc, purity, nmi, ari= evaluate(true_label,pred_label)
type = "ADT"
method_name = "PARC-ADT"
df_pred = pd.DataFrame(pred_label, columns=['pred_label'])
directory_path = os.path.dirname(Benchamark_Result_csv)
pred_file_path = os.path.join(directory_path, "{}_pred.csv".format(method_name))
df_pred.to_csv(pred_file_path, index=False, header=True)

with open(Benchamark_Result_csv, 'a+', newline="") as f:
    csv_write = csv.writer(f)
    data_row = [dataset_name, type, method_name, ari, nmi, acc, purity, memory_peak, Running_time]
    csv_write.writerow(data_row)
