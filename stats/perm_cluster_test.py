#!/usr/bin/env python3
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
from progressbar import ProgressBar, Percentage, Bar, ETA
from time import sleep

import numpy as np
from numpy import ndarray
import random
from scipy.stats import ttest_ind

from pathlib import Path
from typing import List, Tuple, Dict


def cluster_perm_test(X: List[ndarray], threshold: float,
                      n_permutations: int = 1024) -> Tuple[ndarray, List[Dict],
                                                           ndarray, ndarray]:
    """Perform the cluster-based permutation test.
    Args:
        X (List[ndarray]): Each array in X should contain the observations
                           for one group.
    Returns:
        t_val (ndarray): Array of T-values of each datapoint.
        clusters (List[Dict]): Each content corresponds to a cluster.
            Each dictionary has 'start', 'stop', 'sum_tval'.
        cluster_pv (ndarray): Cluster p-values.
        maxstats (mdarray): The maximum sums of T-values from each permutation
    """
    # Run t-test on each datapoint, make clusters,
    # calculate the sums of t-values within each cluster.
    t_val, clusters = make_clusters(X, threshold)
    # If there is no cluster, the algorithm ends here.
    if len(clusters) == 0:
        return t_val, clusters, [], []
    else:
        # Permuate n_permutation - 1 & make MaxStat distribution.
        # MaxStat: The maximum of the sums of T-values within clusters
        #          made in a given permutation.
        maxstats: ndarray = maxstats_distribution(X, threshold, n_permutations)
        # Add the original MaxStat.
        orig = np.abs([cluster['sum_tval'] for cluster in clusters]).max()
        maxstats = np.insert(maxstats, 0, orig)
        # Calculate cluster-pvalue from the percentile of the original
        # cluster-stats in the MaxStats distribution.
        cluster_pv = cluster_pvals(clusters, maxstats)
        return t_val, clusters, cluster_pv, maxstats


def cluster_pvals(clusters: List[Dict],
                  maxstats: ndarray):
    """Check if each cluster is significanct. Cluster-pvalue is calculated by
       acquiring the percentile of the original cluster_stats (sum of T-values)
       in the MaxStats distribution.
    Args:
        clusters (List[Dict]): Each content corresponds to a cluster.
            Each dictionary has 'start', 'stop', 'sum_tval'.
        maxstats (mdarray): The maximum sums of T-values from each permutation
    Returns:
        p_vals: ndarray
    """
    # Array of cluster stats.
    clus_stats = [cluster['sum_tval'] for cluster in clusters]
    # For each cluster (cluster stat), take the ratio of the maxstats in the 
    # MaxStats distribution that exceeds the cluster stat. 
    # This is for two-tailed test as we are looking at the absolute values.
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
    # Run t-test at each data point.
    t_val, _ = ttest_ind(X[0], X[1], axis=0)
    # If two or more neighbouring t-values are above the threshold,
    # they make a cluster.
    return t_val, _make_clusters(t_val, threshold)


def _make_clusters(t_val: ndarray, threshold: float = None) -> List[Dict]:
    """Make clusters by putting neighbours together.
       If two or more neighbouring t-values are above the threshold,
       they make a cluster
    Args:
        t_val (ndarray): Array of T-values of each datapoint.
        threshold (float): Threshold with which clusters will be made.
    Returns:
        clusters[List[Dict]]: Each content corresponds to a cluster.
            Each dictionary has 'start', 'stop', 'sum_tval'.
            start: The index at which the cluster starts.
            stop: The index at which the cluster ends.
            sum_tval: Cluster stat. Sum of T-values within the cluster.
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
    # Absolute
    return np.max([abs(clus['sum_tval']) for clus in clusters])


def maxstats_distribution(X: List[ndarray], threshold: float,
                          n_permutations: int) -> ndarray:
    """Permutate over n_permutation - 1 times and make MaxStats distribution.
    Args:
        X (List[ndarray]): Each array in X should contain the observations
                           for one group.
        threshold (float): Threshold with which clusters will be made.
    Returns:
        maxstats: mdarray: The maximum sums of T-values from each permutation
    """
    maxstats: List[float] = []
    # Permutate n_permutation - 1 times.
#   bar = ProgressBar(
#                     maxval=30,
#                     widgets=[Bar('=', '[', ']'),
#                              ' ',
#                              Percentage()])
#                          #   ' ',
#                          #   ETA()])
#   bar.start()
    for i in range(n_permutations - 1):
        # Shuffle the arrays.
        shuffled_X = _shuffle_arrays(X)
        # Find & make clusters. Get the MaxStat.
        _, clusters = make_clusters(shuffled_X, threshold)
        maxstats.append(find_maxstat(clusters))
#     # if i % 1000 == 0:
#       bar.update(i+1)
#       sleep(0.1)
#   bar.finish()
    # Return MaxStats as ndarray
    return np.array(maxstats)


def _shuffle_arrays(X: List[ndarray]) -> List[ndarray]:
    """Shuffle the data arrays for permutation.
    Args:
        X (List[ndarray]): Each array in X should contain the observations
                           for one group.
    Returns:
        shuffled_X (List[ndarray]): Each array in X should contain randomly
                                    chosen observations.
    """
    # Shuffle the data array
    con: ndarray = np.concatenate(X)
    size = int(con.shape[0])
    indices = list(range(0, size))
    random.shuffle(indices)
    size1 = int(X[0].shape[0])
    # Take the first size1 data as the first group. The rest is the second.
    # Keep the number of subjects inside each group the same as the original.
    return [con[indices[0:size1]], con[indices[size1:]]]
