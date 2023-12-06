import numpy as np
import pandas as pd
import sklearn
from scipy import stats
from sklearn.metrics import pairwise_distances


def get_cluster_centroids(data: np.ndarray, labels: np.ndarray):
    """
    Calculate the centroids of each cluster.

    Args:
        data (np.ndarray): The data points.
        labels (np.ndarray): The cluster labels for each data point.

    Returns:
        np.ndarray: The centroids of each cluster.
    """
    centroids = []
    for label in labels.unique():
        ind = list(np.nonzero(labels == label)[0])
        centroids.append(np.mean(data[ind, :], axis=0))
    centroids = np.row_stack(centroids)
    return centroids


def compute_centroid_distances(data: np.ndarray, labels: np.ndarray, metric: str):
    """
    Compute the pairwise distances between the cluster centroids.

    Parameters:
        data (np.ndarray): The data space used to compute the centroids.
        labels (np.ndarray): Labels of the data points.
        metric (str, optional): The distance metric to be used.

    Returns:
        ndarray: The pairwise distances between the cluster centroids.
    """
    if metric == "angular":
        metric = "cosine"
    assert metric in sklearn.metrics.pairwise.PAIRWISE_DISTANCE_FUNCTIONS
    centroids = get_cluster_centroids(data=data, labels=labels)
    centroid_distances = pairwise_distances(centroids, metric=metric)
    return centroid_distances


def compute_centroid_correlations(hd_centroid_distances, ld_centroid_distances):
    """
    Compute the centroid correlations between high-dimensional and low-dimensional centroid distances.

    Parameters:
    - hd_centroid_distances (np.ndarray): Array of high-dimensional centroid distances.
    - ld_centroid_distances (np.ndarray): Array of low-dimensional centroid distances.

    Returns:
    - corr (np.ndarray): Array of centroid correlations.
    """
    corr = np.empty((hd_centroid_distances.shape[0],))
    for i in range(hd_centroid_distances.shape[0]):
        corr[i] = stats.spearmanr(
            a=hd_centroid_distances[i, :],
            b=ld_centroid_distances[i, :],
            axis=0,
        ).correlation
    return corr


def get_hd_centroid_distance_df(
    labels: pd.Series, label_key: str, hd_centroid_distances: np.ndarray
):
    """
    Store the distances between the centroids of each cluster.

    Args:
        labels (pd.Series): Cluster labels.
        cluster_unique (np.ndarray): Unique cluster labels.
        hd_centroid_distances (np.ndarray): High-dimensional centroid distances.
    """
    labels_unique = labels.unique()
    distance_dict = {label_key: labels_unique}
    distance_dict.update(
        dict(
            (f"{labels_unique[i]}_HD_centroid_dist", hd_centroid_distances[i])
            for i in np.argsort(labels_unique)
        )
    )
    return pd.DataFrame(distance_dict)
