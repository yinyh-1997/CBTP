import numpy as np
from sklearn.metrics.cluster import *
from sklearn import metrics
from scipy.optimize import linear_sum_assignment as linear_assignment
from time import time as get_time
import os


def cluster_acc(y_true, y_pred):
    """
    Calculate clustering accuracy. Require scikit-learn installed
    # Arguments
        y: true labels, numpy.array with shape `(n_samples,)`
        y_pred: predicted labels, numpy.array with shape `(n_samples,)`
    # Return
        accuracy, in [0,1]
    """
    y_true = y_true.astype(np.int64)
    # assert y_pred.size == y_true.size
    assert len(y_pred) == len(y_true)
    # D = max(y_pred.max(), y_true.max()) + 1
    D = max(max(y_pred), max(y_true)) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1

    #from sklearn.utils.linear_assignment_ import linear_assignment
    ind = linear_assignment(w.max() - w)
    ind = np.array((ind[0], ind[1])).T
    # pdb.set_trace()
    return sum([w[i, j] for i, j in ind]) * 1.0 / y_pred.size

def cluster_purity(y_true, y_pred):
    # Preprocessing: Ensure y_true and y_pred are numpy arrays and have the same length
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    assert y_true.shape == y_pred.shape

    # Get the number of predicted clusters
    unique_clusters = np.unique(y_pred)
    n_clusters = len(unique_clusters)

    # If there are no clusters, return 0
    if n_clusters == 0: 
        return 0

    # Initialize an array to store the labels for each cluster
    cluster_labels = np.zeros(n_clusters)
    for i, cluster in enumerate(unique_clusters):
        # For each cluster, find the most common true class
        cluster_mask = y_pred == cluster
        # If the cluster is empty, skip
        if np.sum(cluster_mask) == 0:
            continue
        cluster_labels[i] = np.bincount(y_true[cluster_mask]).argmax()

    # Calculate the number of correctly classified samples
    n_correct = 0
    for i in range(len(y_true)):
        # If the predicted cluster is empty, skip
        if y_pred[i] not in unique_clusters:
            continue
        n_correct += (y_true[i] == cluster_labels[np.where(unique_clusters == y_pred[i])[0][0]])

    # Calculate purity
    purity = n_correct / len(y_true)

    return purity

def evaluate(y_true, y_pred):
    acc = cluster_acc(y_true, y_pred)
    purity = cluster_purity(y_true, y_pred)
    nmi = normalized_mutual_info_score(y_true, y_pred)
    ari = adjusted_rand_score(y_true, y_pred)
    return acc, purity, nmi, ari