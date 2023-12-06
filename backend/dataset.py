import numpy as np
import pandas as pd
import os
import neighbors
import time
import glasbey
import anndata as ad
from collections import defaultdict
import sklearn
import centroid_correlation


from abc import ABC, abstractmethod


class Dataset(ABC):
    """
    Embeddings are stored in adata.obsm[methodName_methodIndex]
    and should be normalized to [-1,1].

    The quality scores are stored in adata.uns[methodName_methodIndex][quality_score_name]

    The dictionary adata.uns['methods'] contains the available embeddings
    as keys and the number of different embeddings from that method as values.

    """

    def __init__(
        self,
        filepath: str,
        name: str,
        hd_label_key: str,
        hd_metric: str,
        description: str,
    ):
        self.filepath = filepath
        self.name = name
        self.hd_label_key = hd_label_key
        self.hd_metric = hd_metric
        self.description = description
        self.hd_annoy_filepath = dict()

        if not os.path.isfile(self.filepath):
            raise FileNotFoundError(f"File {self.filepath} not found.")
        self.adata = ad.read_h5ad(self.filepath)
        if "hd_neighbors" not in self.adata.uns_keys():
            self.adata.uns["hd_neighbors"] = dict()

        print(
            f"Loaded {self.name} data with embeddings {list(self.adata.uns['methods'].keys())}"
        )

        if self.adata.uns["hd_neighbors"].get(self.hd_metric, None) is None:
            print("Computing global HD neighbors")
            self.precompute_HD_neighbors(maxK=200, metric=self.hd_metric)

        # which quality measures are available?
        self.quality_features = set()
        emb_names = self.get_embedding_names()
        for name in emb_names:
            if name in self.adata.uns_keys():
                self.quality_features.union(self.adata.uns[name].keys())
            else:
                self.adata.uns[name] = dict()

        # global correlations
        if "quality centroid_corr" not in self.quality_features:
            print(
                f"Computing global quality: centroid distance correlation with {self.hd_label_key} labels"
            )
            self.compute_centroid_distances(self.hd_metric, self.hd_label_key)
            self.quality_features.add("quality centroid_corr")

        # make sure the local quality is computed for all embeddings
        print("Checking local quality scores.")
        self.compute_quality(k=50, hd_metric=self.hd_metric)
        self.compute_quality(k=200, hd_metric=self.hd_metric)
        self.quality_features.add("quality qnx@50")
        self.quality_features.add("quality qnx@200")
        self.quality_features = list(self.quality_features)

    @abstractmethod
    def get_metadata_features(self):
        """
        List of column names that can be used to color points.
        """
        pass

    @abstractmethod
    def get_HD_data(self):
        """Return the ndarray to use for HD neighbor computation.
        e.g. obsm["X_pca"]
        """
        pass

    def save_adata(self, overwrite: bool = False):
        """
        Save the AnnData object to a file.

        Parameters:
            overwrite (bool): If True, overwrite the existing file. If False, create a new file with a timestamp in the filename.
        """
        if overwrite:
            self.adata.write(self.filepath, compression="gzip")
        else:
            new_filename = (
                os.path.splitext(os.path.basename(self.filepath))[0]
                + f"_{time.strftime('%Y%m%d_%H%M%S')}.h5ad"
            )
            self.adata.write(
                os.path.join(os.path.dirname(self.filepath), new_filename),
                compression="gzip",
            )
            print(f"Saved dataset to {new_filename}")

    def build_annoy_index(self, metric: str):
        """
        Builds an Annoy index for the dataset using the specified metric.

        Parameters:
            metric (str): The metric to be used for building the index.
        """
        self.hd_annoy_filepath[metric] = neighbors.build_annoy_index(
            self.get_HD_data(),
            metric=metric,
            filepath=os.path.join(
                os.path.dirname(self.filepath),
                self.name.replace(" ", "_") + "_" + metric + "_annoy.ann",
            ),
        )

    def get_embedding_name(self, embName, embScale):
        """
        Defines the string that describes the embedding of a certain scale.
        This name is used to store the embedding in self.adata.obsm

        Args:
            embName (str): The name of the embedding.
            embScale (int): The scale of the embedding.

        Returns:
            str: The string of the embedding name and scale.
        """
        return embName + "_" + str(embScale)

    def cleanup(self):
        """
        Clean up the dataset by removing any existing HD_annoy files.

        This method iterates over the HD_annoy filepaths and removes the corresponding files if they exist.
        """
        if len(self.hd_annoy_filepath) > 0:
            for fpath in self.hd_annoy_filepath.values():
                if fpath is not None and os.path.isfile(fpath):
                    os.remove(fpath)

    def get_embedding_options(self):
        """
        Returns a dictionary stating how many embeddings are available for each method.

        Returns:
            dict: A dictionary where the keys are the method names and the values are the number of embeddings available for each method.
                  Example: {"tSNE": 2, "UMAP": 5}
        """
        # must ensure that the numbers are integers and not numpy types
        for k, v in self.adata.uns["methods"].items():
            self.adata.uns["methods"][k] = int(v)
        return self.adata.uns["methods"]

    def get_embedding_names(self):
        """
        Get the names of all embeddings in the dataset.

        Returns:
            list: A list of embedding names to be used to index into adata.obsm.
        """
        names = []

        for embedding_method in self.adata.uns["methods"]:
            for embedding_scale in range(self.adata.uns["methods"][embedding_method]):
                names.append(self.get_embedding_name(embedding_method, embedding_scale))
        return names

    def get_category_colors(self, fname):
        """
        Get the colors associated with each category in the specified feature.

        Args:
            fname (str): The name of the feature.

        Returns:
            dict: A dictionary mapping each category to its corresponding color.
        """
        categories = list(
            self.adata.obs[fname].value_counts(ascending=False, sort=True).keys()
        )
        if fname + "_colors" in self.adata.uns_keys():
            colors = self.adata.uns[fname + "_colors"][: len(categories)]
        else:
            colors = glasbey.create_palette(
                palette_size=len(categories), grid_size=32, colorblind_safe=True
            )
        return dict(zip(categories, colors))

    def print_quality(self):
        """
        Print the quality of the dataset for each method.

        This method iterates over the methods in the dataset and prints the quality
        information for each method. It calculates the mean value and the percentage
        of zeros for each embedding in the dataset.

        Returns:
            None
        """
        for method in self.adata.uns["methods"]:
            print(f"Quality for {method}:")
            for i in range(self.adata.uns["methods"][method]):
                name = self.get_embedding_name(method, i)
                for k in self.adata.uns[name].keys():
                    print(
                        f"  {i}: {self.adata.uns[name][k].mean():.3f}, "
                        f"{100*np.sum(self.adata.uns[name][k] == 0.0)/self.adata.n_obs:.3f}% zero"
                    )

    def precompute_HD_neighbors(self, maxK: int, metric: str, exact=False):
        """
        Precomputes the high-dimensional (HD) neighbors for each data point.

        Args:
            maxK (int): The maximum number of neighbors to compute.
            metric (str): The distance metric to use for computing neighbors.
            exact (bool, optional): Whether to compute exact neighbors or use an approximate algorithm. Defaults to False.
        """
        if (
            metric in self.adata.uns["hd_neighbors"]
            and self.adata.uns["hd_neighbors"][metric].shape[1] < maxK
        ) or (metric not in self.adata.uns["hd_neighbors"]):
            print(f"Computing all HD neighbors until k = {maxK}")
            start = time.time()
            self.adata.uns["hd_neighbors"][metric] = neighbors.get_nearest_neighbors(
                data=self.get_HD_data(),
                indices=range(self.adata.n_obs),
                k=maxK,
                metric=metric,
                filepath=self.hd_annoy_filepath.get(metric, None),
                exact=exact,
            )
            print(f"Computing neighbors in {time.time() - start:.2f} seconds")
            self.save_adata()

    def get_HD_neighbors(self, k: int, metric: str, indices: np.ndarray):
        """
        Get the high-dimensional neighbors for the given indices.

        Parameters:
            k (int): The number of neighbors to retrieve.
            metric (str): The distance metric to use for neighbor search.
            indices (np.ndarray): The indices of the data points for which to find neighbors.

        Returns:
            np.ndarray: An array of shape (len(indices), k) with the high-dimensional neighbors
            for the given indices.
        """
        if (
            metric in self.adata.uns["hd_neighbors"]
            and self.adata.uns["hd_neighbors"][metric].shape[1] >= k
        ):
            return self.adata.uns["hd_neighbors"][metric][indices, :k]
        elif indices.shape[0] > self.adata.n_obs / 2:
            # precompute all neighbors
            self.precompute_HD_neighbors(maxK=k, metric=metric)
            return self.adata.uns["hd_neighbors"][metric][indices, :k]
        else:
            # compute neighbors on the fly
            return neighbors.get_nearest_neighbors(
                data=self.get_HD_data(),
                indices=indices,
                k=k,
                metric=metric,
                filepath=self.hd_annoy_filepath.get(metric, None),
            )

    def compute_quality(self, k: int, hd_metric: str, ld_metric: str = "euclidean"):
        """
        Computes the quality scores for embeddings based on neighborhood preservation.
        It stores the quality scores in the adata.uns dictionary of every embedding.

        Args:
            k (int): The number of nearest neighbors to consider.
            hd_metric (str): The metric to use for high-dimensional space.
            ld_metric (str, optional): The metric to use for low-dimensional space. Defaults to "euclidean".
        """
        if (
            hd_metric not in self.adata.uns["hd_neighbors"]
            or self.adata.uns["hd_neighbors"][hd_metric].shape[1] < k
        ):
            print("Computing HD neighbors before calculating quality")
            self.precompute_HD_neighbors(maxK=k, hd_metric=hd_metric)

        embedding_names = []
        embeddings = []
        quality_key = f"quality qnx@{k}"
        for name in self.get_embedding_names():
            if (
                name in self.adata.uns_keys()
                and self.adata.uns[name].get(quality_key, None) is None
            ):
                print(f"adding embedding {name}")
                embedding_names.append(name)
                embeddings.append(self.adata.obsm[name])

        if len(embeddings) > 0:
            print(f"Computing quality for {embedding_names}")
            quality_scores = neighbors.neighborhood_preservation(
                hd_data=self.get_HD_data(),
                ld_data_arr=embeddings,
                k=k,
                hd_metric=hd_metric,
                ld_metric=ld_metric,
                hd_annoy_filepath=self.hd_annoy_filepath.get(hd_metric, None),
                hd_neighbors=self.adata.uns["hd_neighbors"][hd_metric],
            )

            for name, quality_result in zip(embedding_names, quality_scores):
                self.adata.uns[name][quality_key] = quality_result
            self.save_adata()

    def compute_centroid_distances(self, hd_metric: str, label_key: str):
        """
        Compute the correlation between centroid distances in HD and LD space.

        Parameters:
            hd_metric (str): The metric used to compute distances between HD data points.
            label_key (str): Key for cluster labels in adata.obs.
        """
        hd_centroid_distances = centroid_correlation.compute_centroid_distances(
            data=self.get_HD_data(),
            labels=self.adata.obs[label_key],
            metric=hd_metric,
        )

        # store distances in adata.obs
        self.adata.obs = pd.merge(
            self.adata.obs,
            centroid_correlation.get_hd_centroid_distance_df(
                labels=self.adata.obs[label_key],
                label_key=label_key,
                hd_centroid_distances=hd_centroid_distances,
            ),
            on=label_key,
            how="left",
        )

        # compute correlations to ld centroids
        embedding_names = self.get_embedding_names()
        for name in embedding_names:
            ld_centroid_distances = centroid_correlation.compute_centroid_distances(
                data=self.adata.obsm[name],
                labels=self.adata.obs[label_key],
                metric="euclidean",
            )
            corr = centroid_correlation.compute_centroid_correlations(
                hd_centroid_distances, ld_centroid_distances
            )
            centroid_correlations = pd.merge(
                pd.DataFrame({"cluster": self.adata.obs[label_key]}),
                pd.DataFrame(
                    {
                        "cluster": self.adata.obs[label_key].unique(),
                        "centroid_correlation": corr,
                    }
                ),
                on="cluster",
                how="left",
            )["centroid_correlation"].values
            self.adata.uns[name]["quality centroid_corr"] = centroid_correlations
