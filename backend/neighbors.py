from annoy import AnnoyIndex
import numpy as np
import os
import time
import sklearn
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm


def build_annoy_index(
    data: np.ndarray,
    metric: str,
    ntrees: int = 100,
    filepath: str = None,
):
    """
    Build an Annoy index for nearest neighbor search.

    Args:
        data (np.ndarray): The data points to build the index on.
        metric (str): The distance metric to use for nearest neighbor search.
            Possible values are "angular", "euclidean", "manhattan", "hamming", or "dot".
        ntrees (int, optional): The number of trees to build in the index. Defaults to 100.
        filepath (str, optional): The filepath to save the index. If not provided, a default filepath
            will be generated based on the distance metric and current timestamp.

    Returns:
        str: The filepath where the index is saved.
    """
    if filepath is None:
        filepath = os.path.join(
            ".", f"annoy_{metric}_{time.strftime('%Y%m%d_%H%M%S')}.ann"
        )

    if not os.path.isfile(filepath):
        data = np.asarray(data)
        n, d = data.shape
        t = AnnoyIndex(d, metric)
        for i in range(n):
            t.add_item(i, data[i, :])
        t.build(ntrees)
        t.save(filepath)
    return filepath


def get_nearest_neighbors(
    data: np.ndarray,
    indices: np.ndarray,
    k: int,
    metric: str,
    filepath: str = None,
    exact: bool = False,
):
    """
    Get the nearest neighbors for given data points.

    Args:
        data (np.ndarray): The data points to find nearest neighbors for.
        indices (np.ndarray): The indices of the data points to find nearest neighbors for.
        k (int): The number of nearest neighbors to retrieve.
        metric (str): The distance metric to use for finding nearest neighbors.
            Supported distance metrics (Annoy) are ['angular', 'euclidean', 'manhattan', 'hamming', 'dot'].
            Supported distance metrics for exact neighbors (KD-tree) are
            ['euclidean', 'l2', 'minkowski', 'p', 'manhattan', 'cityblock', 'l1', 'chebyshev', 'infinity']
        filepath (str, optional): The file path to the pre-built index. If None, a new index will be built.
        exact (bool, optional): Whether to use exact nearest neighbor search or approximate search.
            Defaults to False.

    Returns:
        np.ndarray: An array of shape (len(indices), k) containing the indices of the k nearest neighbors for each data point.
    """
    if exact:
        return get_exact_neighbors(
            data,
            indices,
            k,
            metric,
        )
    else:
        if filepath is None or not os.path.isfile(filepath):
            filepath = build_annoy_index(
                data,
                metric,
            )

        data = np.asarray(data)
        d = data.shape[1]
        u = AnnoyIndex(d, metric)
        u.load(filepath)

        if len(indices) == 1:
            return u.get_nns_by_item(indices[0], k + 1)[1:]
        else:
            assert k < data.shape[0]
            nbrs = np.empty((len(indices), k), dtype=int)
            for i in tqdm(indices):
                nbrs[i, :] = u.get_nns_by_item(i, k + 1)[1:]
            return nbrs


def get_exact_neighbors(
    data: np.ndarray,
    indices: np.ndarray,
    k: int,
    metric: str,
):
    """
    Get the exact neighbors for the given indices in the data array.

    Args:
        data (np.ndarray): The input data array.
        indices (np.ndarray): The indices for which to find neighbors.
        k (int): The number of neighbors to retrieve.
        metric (str): The distance metric to use for neighbor search.
            Supported distance metrics (KD-tree) are
            ['euclidean', 'l2', 'minkowski', 'p', 'manhattan', 'cityblock', 'l1', 'chebyshev', 'infinity']

    Returns:
        np.ndarray: The indices of the exact neighbors for the given input indices.
    """
    if metric not in sklearn.neighbors.KDTree.valid_metrics():
        raise ValueError(f"Distance metric {metric} not supported by KDTree")

    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm="auto", metric=metric).fit(data)
    neighbor_indices = nbrs.kneighbors(data[indices, :], return_distance=False)[:, 1:]
    return neighbor_indices


def neighborhood_preservation(
    hd_data: np.ndarray,
    ld_data_arr: list[np.ndarray],
    k: int,
    hd_metric: str,
    ld_metric: str,
    hd_neighbors: np.ndarray = None,
    hd_annoy_filepath: str = None,
):
    """
    Computes the neighborhood preservation quality for given high-dimensional data and a list of low-dimensional embeddings.

    Args:
        hd_data (np.ndarray): The high-dimensional dataset.
        ld_data_arr (list[np.ndarray]): The list of low-dimensional datasets.
        k (int): The number of neighbors to consider.
        hd_metric (str): The distance metric to use for the high-dimensional dataset.
        ld_metric (str): The distance metric to use for the low-dimensional datasets.
        hd_neighbors (np.ndarray, optional): The precomputed nearest neighbors for the high-dimensional dataset.
            If the hd_neighbors are provided, the hd_annoy_filepath is not needed. Defaults to None.
        hd_annoy_filepath (str, optional): The file path to the precomputed Annoy index for the high-dimensional dataset. Defaults to None.

    Returns:
        list[np.ndarray]: The neighborhood preservation quality in [0,1] for each point and low-dimensional dataset.
    """
    n = hd_data.shape[0]

    if hd_neighbors is None:
        if hd_annoy_filepath is None or not os.path.isfile(hd_annoy_filepath):
            hd_annoy_filepath = build_annoy_index(
                hd_data,
                hd_metric,
            )

        hd_annoy = AnnoyIndex(hd_data.shape[1], hd_metric)
        hd_annoy.load(hd_annoy_filepath)
    else:
        assert hd_neighbors.shape[1] >= k

    emb_quality = [np.zeros(n) for _ in ld_data_arr]

    # compute exact ld neighbors
    emb_neighbors = []
    for e, embedding in enumerate(ld_data_arr):
        emb_neighbors.append(
            get_exact_neighbors(
                data=embedding,
                indices=np.arange(n),
                k=k,
                metric=ld_metric,
            )
        )

    # compute set intersection size
    for i in tqdm(range(n)):
        if hd_neighbors is not None:
            hd_nn = hd_neighbors[i, 0:k]
        else:
            hd_nn = hd_annoy.get_nns_by_item(i, n=k + 1)[1:]
        hd_set = set(hd_nn)

        for e in range(len(ld_data_arr)):
            ld_nn = emb_neighbors[e][i]
            emb_quality[e][i] = len(hd_set.intersection(ld_nn)) / k

    return emb_quality
