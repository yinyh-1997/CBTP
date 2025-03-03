from pathlib import Path
import os
import subprocess
import time
import csv
import pandas as pd
import numpy as np

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

current_timestamp = time.time()
current_time_str = time.strftime('%Y-%m-%d_%H:%M:%S', time.localtime(current_timestamp))
dataset_path = Path('/home/yhyin/code/Benchmarking/IntegrationDataset')
result_Path = Path('/home/yhyin/code/Benchmarking/BenchmarkResult_integration')
dataset_list = sorted(os.listdir(dataset_path))

benchmark_randomseed = 2024
file_path = '/home/yhyin/code/Benchmarking/Datasets_information.xlsx'
df = pd.read_excel(file_path, engine='openpyxl')

for dataset_ in dataset_list:

    filtered_row = df[df['Dataset'] == dataset_]
    Integration_num_value = int(filtered_row['NumIntegrationFeatures'].values[0])
    Integration_num_value_CarDEC = Integration_value-1

    dataset = os.path.join(dataset_path, dataset_)
    cell_name = pd.read_csv(dataset + '/celltypes.csv', usecols=['celltype']).values
    unique_values, _ = np.unique(cell_name, return_inverse=True)
    n_cluster = len(unique_values)
    result_csv_dir = os.path.join(result_Path, current_time_str, dataset_)
    if not os.path.exists(result_csv_dir):
        os.makedirs(result_csv_dir)
    Benchamark_Result_csv = os.path.join(result_csv_dir, "Result.csv")
    with open(Benchamark_Result_csv, 'w+', newline="") as f:
        csv_write = csv.writer(f)
        csv_head = ["Dataset","Type","Method","ARI","NMI","ACC","Purity","Memory_peak","Running_time"]
        csv_write.writerow(csv_head)

    ##Louvain
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/Louvain')
    subprocess.run(f'conda run -n louvain-env python louvain-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {benchmark_randomseed} {dataset}', shell=True)

    ##Leiden
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/Leiden')
    subprocess.run(f'conda run -n leiden-env python leiden-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {benchmark_randomseed} {dataset}', shell=True)

    ##FFC
    PCA_dim = 50
    if(Integration_value <50):
        if(Integration_value==2):
            PCA_dim = Integration_value
        else:
            PCA_dim = Integration_value - 1
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/FFC')
    subprocess.run(f'conda run -n ffc-env python ffc-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {benchmark_randomseed} {PCA_dim} {dataset}', shell=True)

    ##PARC
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/PARC')
    subprocess.run(f'conda run -n parc-env python parc-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {benchmark_randomseed} {dataset}', shell=True)

    ## DESC
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/DESC')
    subprocess.run(f'conda run -n desc-env python desc-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {benchmark_randomseed} {dataset}', shell=True)

    ##CarDEC
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/CarDEC')
    subprocess.run(f'conda run -n cardec-env python cardec-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {benchmark_randomseed} {Integration_value_CarDEC} {dataset}', shell=True)

    ##SCHNEL
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/SCHNEL')
    subprocess.run(f'conda run -n schnel-env python schnel-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {dataset}', shell=True)

    ##PhenoGraph
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/PhenoGraph')
    subprocess.run(f'conda run -n phenograph-env python phenograph-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {dataset}', shell=True)
    
    ##scAIDE
    os.chdir('/home/yhyin/code/Benchmarking/ClusteringAlgorithm/scAIDE')
    subprocess.run(f'conda run -n scaide-env python scaide-integration_feature.py {dataset_} {n_cluster} {Benchamark_Result_csv} {dataset}', shell=True)
