"""Cluster-based permutation statistics for NIRS channel-level simulations."""

import numpy as np
import mne

def ftest_rel_no_p_custom(*args):
    """
    Compute a paired-samples t-statistic and convert it to F.

    """
    condition1, condition2 = args[0], args[1]
    diff = condition1 - condition2
    var = np.var(diff, axis=0, ddof=1)

    sigma = 0.0
    limit = sigma * np.max(var) if var.size else 0.0
    var = var + limit

    mean_diff = np.mean(diff, axis=0)
    t_stat = mean_diff / np.sqrt(var / diff.shape[0])
    return np.square(t_stat)


def run_cluster_permutation_test(
    X,
    adjacency_sparse,
    *,
    n_permutations=1024,
    alpha=0.05,
    stat_fun=ftest_rel_no_p_custom,
    n_jobs=1,
    seed=None,
):
    """
    Run an MNE cluster-based permutation test and identify significant channels.

    Parameters
    ----------
    X : list of ndarray
        Two input groups, each shaped (n_subjects, n_channels).
    adjacency_sparse : sparse matrix
        Channel adjacency matrix used for clustering.
    n_permutations : int
        Number of label permutations.
    alpha : float
        Cluster significance threshold.
    stat_fun : callable
        Statistic function used by MNE.
    n_jobs : int
        Number of MNE worker jobs. Keep this at 1 inside process workers.
    seed : int or None
        Random seed used for MNE's internal permutations.
    """
    t_obs, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_test(
        X,
        n_permutations=n_permutations,
        adjacency=adjacency_sparse,
        stat_fun=stat_fun,
        n_jobs=n_jobs,
        seed=seed,
    )

    sig_clusters = np.where(cluster_p_values < alpha)[0]
    clusters_temp = np.array([cluster[0] for cluster in clusters], dtype=object)

    is_valid = bool(clusters_temp[sig_clusters].size) and any(
        subarray.size for subarray in clusters_temp[sig_clusters]
    )

    significant_channels = (
        np.concatenate(clusters_temp[sig_clusters]) if is_valid else None
    )

    return t_obs, clusters, cluster_p_values, H0, significant_channels
