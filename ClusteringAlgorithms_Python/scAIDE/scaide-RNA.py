"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy.api as sc
from sklearn.utils import check_array

from aide import AIDE, AIDEConfig
from rph_kmeans import RPHKMeans
from sklearn.cluster import KMeans

from metric_compute import cluster_acc, cluster_purity, evaluate
import os
import csv
from time import time as get_time
import tracemalloc
import sys
from pathlib import Path


def to_adata(features, labels=None):
	adata = sc.AnnData(features)
	if labels is not None:
		adata.obs['Group'] = labels
	return adata


def from_adata(adata):
	features, labels = adata.X, None
	if 'Group' in adata.obs:
		labels = adata.obs['Group'].values
	return features, labels


def preprocess(X, y_true=None, x_dtype=np.float32, y_dtype=np.int32):
	"""
	Args:
		X (np.ndarray or sp.csr_matrix): raw counts data
			dtype = (np.float32 | np.float64 | np.int32 | np.int64)
			shape = (cell_num, gene_num)
		y_true (np.ndarray or None): labels of X
			dtype = (np.int32 | np.int64)
			shape = (cell_num,)
		x_dtype (type): dtype of preprocessed X; should be one of (np.float32, np.float64)
		y_dtype (type): dtype of preprocessed Y; should be one of (np.int32, np.int64)
	Returns:
		np.ndarray or sp.csr_matrix: preprocessed X
			type = type(X)
			dtype = x_dtype
			shape = (filtered_cell_num, filtered_gene_num)
		np.ndarray or None: preprocessed y_true
			dtype = y_dtype
			shape = (filtered_cell_num,)
	"""
	assert x_dtype == np.float32 or x_dtype == np.float64
	assert y_dtype == np.int32 or y_dtype == np.int64
	adata = to_adata(X, y_true)
	sc.pp.filter_cells(adata, min_counts=1)
	sc.pp.filter_genes(adata, min_counts=1)
	sc.pp.normalize_per_cell(adata)
	sc.pp.log1p(adata)
	X, y_true = from_adata(adata)

	X = check_array(X, accept_sparse="csr", order='C', dtype=[x_dtype])
	if sp.issparse(X):
		X.sorted_indices()
	if y_true is not None:
		y_true = check_array(y_true, ensure_2d=False, order='C', dtype=[y_dtype])
	return X, y_true


if __name__ == '__main__':
	import h5py
	import os
	import itertools
	from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

	DEMO_FOLDER = os.path.dirname(os.path.realpath(__file__))

	dataset_name = sys.argv[1]
	n_cluster = int(sys.argv[2])
	Benchamark_Result_csv = sys.argv[3]
	Benchamark_data_path = sys.argv[4]

	x = pd.read_csv(Benchamark_data_path + '/RNA.csv', header=0, index_col=0).T.values
	cell_name = pd.read_csv(Benchamark_data_path + '/celltypes.csv', usecols=['celltype']).values
	_, y_true = np.unique(cell_name, return_inverse=True)
	adata = sc.AnnData(x, dtype=np.float32)
	
	temp_top_genes=2000
	if adata.X.shape[1]<2000:
		temp_top_genes = adata.X.shape[1]
	
	sc.pp.highly_variable_genes(adata, n_top_genes=temp_top_genes, subset=True, inplace=True)
	X = adata.X

	tracemalloc.start()
	start_time = get_time()

	aide_model = AIDE(name='aide_for_bladder', save_folder='aide_for_bladder')
	config = AIDEConfig()   # Run with default config
	embedding = aide_model.fit_transform(X, config)
	print(f'Type of embedding = {type(embedding)}; Shape of embedding = {embedding.shape}; Data type of embedding = {embedding.dtype}')
	clt = RPHKMeans(n_clusters=n_cluster, n_init=1, verbose=0)
	y_pred = clt.fit_predict(embedding)

	end_time = get_time()
	Running_time = end_time - start_time
	memory_current, memory_peak = tracemalloc.get_traced_memory()
	tracemalloc.stop()
	memory_peak = memory_peak/(1024*1024)
	
	scaide_labels = y_pred
	true_label = np.array(y_true)
	pred_label = np.array(scaide_labels)

	acc, purity, nmi, ari= evaluate(true_label,pred_label)
	type = "RNA"
	method_name = "scAIDE-RNA"
	df_pred = pd.DataFrame(pred_label, columns=['pred_label'])
	directory_path = os.path.dirname(Benchamark_Result_csv)
	pred_file_path = os.path.join(directory_path, "{}_pred.csv".format(method_name))
	df_pred.to_csv(pred_file_path, index=False, header=True)

	with open(Benchamark_Result_csv, 'a+', newline="") as f:
		csv_write = csv.writer(f)
		data_row = [dataset_name, type, method_name, ari, nmi, acc, purity, memory_peak, Running_time]
		csv_write.writerow(data_row)
