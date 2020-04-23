#!/usr/bin/env python
"""
Author: Yuki Fujishima <yfujishima1001@gmail.com>

Method developed by
    Maris, E., & Oostenveld, R. (2007)
    Nonparametric statistical testing of EEG- and MEG-data.
    Journal of Neuroscience Methods, 164(1), 177–190.
    https://doi.org/10.1016/J.JNEUMETH.2007.03.024

""" 

import json
import matplotlib.pyplot as plt

import numpy as np
from numpy import ndarray
import random
from scipy.stats import ttest_ind

from pathlib import Path
from typing import List, Tuple, Dict


def cluster_perm_test(X: List[ndarray], threshold: float,
                      n_permutations: int = 1024):
    """Perform the cluster-based permutation test.
    Args:
        data (List[ndarray]): Each array in X should contain the observations
                              for one group.
    Returns:
        t_val (ndarray): Array of T-values of each datapoint.
        clusters (List[Dict]): Each content corresponds to a cluster.
            Each dictionary has 'start', 'stop', 'sum_tval'.
        cluster_pv (ndarray): Cluster p-values.
        maxstats (mdarray): The maximum sums of T-values from each permutation
    """
    t_val, clusters = make_clusters(X, threshold)
    if len(clusters) == 0:
        return t_val, clusters, [], []
    else:
        maxstats: ndarray = maxstats_distribution(X, threshold, n_permutations)
        orig = np.abs([cluster['sum_tval'] for cluster in clusters]).max()
        maxstats = np.insert(maxstats, 0, orig)
        cluster_pv = cluster_pvals(clusters, maxstats)
        return t_val, clusters, cluster_pv, maxstats


def cluster_pvals(clusters: List[Dict],
                  maxstats: ndarray):
    """Check if each cluster is significanct.
    Args:
        clusters (List[Dict]): Each content corresponds to a cluster.
            Each dictionary has 'start', 'stop', 'sum_tval'.
        maxstats (mdarray): The maximum sums of T-values from each permutation
    Returns:
        p_vals: ndarray
    """
    clus_stats = [cluster['sum_tval'] for cluster in clusters]
    p_vals = np.array([np.mean(abs(maxstats) >= abs(t)) for t in clus_stats])
    return p_vals


def make_clusters(X: List[ndarray],
                  threshold: float = None) -> Tuple[ndarray, List[Dict]]:
    """Make clusters by t-test.
    Args:
        X (List[ndarray]): Each array in X should contain the observations
                              for one group.
        threshold (float): Threshold with which clusters will be made.
        
    Returns:
        clusters[List[Dict]]: Each content corresponds to a cluster.
            Each dictionary has 'start', 'stop', 'sum_tval'.
    """
    # test at each data point.
    t_val, _ = ttest_ind(X[0], X[1], axis=0)
    return t_val, _make_clusters(t_val, threshold)


def _make_clusters(t_val: ndarray, threshold: float = None) -> List[Dict]:
    """Make clusters by putting neighbours together.
    Args:
        t_val (ndarray): Array of T-values of each datapoint.
        threshold (float): Threshold with which clusters will be made.
    Returns:
        clusters[List[Dict]]: Each content corresponds to a cluster.
            Each dictionary has 'start', 'stop', 'sum_tval'.
    """
    clusters: List[Dict] = []
    prev_sig: bool = False
    for pos, t in enumerate(t_val):
        if prev_sig is False:
            clus_idx: List[int] = []
            clus: List[float] = []
            if np.abs(t) >= threshold:
                clus_idx.append(pos)
                clus.append(t)
                prev_sig = True
            elif np.abs(t) < threshold:
                pass
            else:
                raise ValueError
        elif prev_sig is True:
            if np.abs(t) >= threshold:
                clus_idx.append(pos)
                clus.append(t)
            elif np.abs(t) < threshold:
                clusters.append({'start': clus_idx[0],
                                 'stop': clus_idx[-1] + 1,
                                 'sum_tval': np.sum(clus)})
                prev_sig = False
            else:
                raise ValueError
    return clusters


def find_maxstat(clusters: List[Dict]) -> float:
    """Find the biggest sum of t-vales in a given set of clusters.
    Returns:
        maxstat: float = The biggest sum of t-values.
    """
    if len(clusters) == 0:
        return 0
    return np.max([clus['sum_tval'] for clus in clusters])


def maxstats_distribution(X: List[ndarray], threshold: float,
                          n_permutations: int) -> ndarray:
    """Permutate over N times and make MaxStats distribution.
    Args:
        X (List[ndarray]): Each array in X should contain the observations
                           for one group.
        threshold (float): Threshold with which clusters will be made.
    Returns:
        maxstats: mdarray: The maximum sums of T-values from each permutation
    """
    maxstats: List[float] = []
    for i in range(n_permutations - 1):
        shuffled_X = _shuffle_arrays(X)
        _, clusters = make_clusters(shuffled_X, threshold)
        maxstats.append(find_maxstat(clusters))
    return np.array(maxstats) 


def _permutate(X: List[ndarray], threshold: float = None):
    shuffled_X = _shuffle_arrays(X)
    _, clusters = make_clusters(shuffled_X, threshold)
    return find_maxstat(clusters)


def _shuffle_arrays(X: List[ndarray]) -> List[ndarray]:
    size1 = int(X[0].shape[0])
    con: ndarray = np.concatenate(X)
    size = int(con.shape[0])
    indices = list(range(0, size))
    random.shuffle(indices)
    return [con[indices[0:size1]], con[indices[size1:]]]
