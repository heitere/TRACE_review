from dataset import Dataset


class HumanImmune(Dataset):
    def __init__(self):
        super().__init__(
            "../data/immune/immune_with_embeddings_new_20231031_172941.h5ad",
            "Human Immune",
            hd_label_key="celltype",
            hd_metric="euclidean",
            description="Human Immune, Luecken et al. (2022)",
        )

    def get_metadata_features(self):
        cols = [
            "celltype",
            "batch",
            "chemistry",
            "species",
            "study",
            "tissue",
        ]
        return cols

    def get_HD_data(self):
        return self.adata.obsm["X_pca"]


class MouseCD45neg(Dataset):
    def __init__(self):
        super().__init__(
            filepath="../data/mouse_fibro/CD45_PCA_init_quality.h5ad",
            name="Mouse CD45neg",
            hd_label_key="louvain",
            hd_metric="angular",
            description="Mouse cells, CD45neg (from Guilliams M. et al. 2022)",
        )

    def get_metadata_features(self):
        cols = ["annot", "louvain", "type", "sampleName"]
        cols = self.adata.obs_keys()
        cols.remove("annotID")
        cols.remove("sample")
        cols.remove("sampleID")
        cols.remove("cell")
        return cols

    def get_HD_data(self):
        return self.adata.obsm["X_totalVI"]


class MouseFibro(Dataset):
    def __init__(self):
        super().__init__(
            filepath="../data/mouse_fibro/mouseFibroblasts.h5ad",
            name="Mouse Fibroblasts",
            hd_label_key="louvain",
            hd_metric="angular",
            description="Mouse cells, Fibroblasts subset (from Guilliams M. et al. 2022)",
        )

    def get_metadata_features(self):
        cols = ["annot", "louvain", "type", "sampleName"]
        cols = self.adata.obs_keys()
        cols.remove("annotID")
        cols.remove("sample")
        cols.remove("sampleID")
        cols.remove("cell")
        return cols

    def get_HD_data(self):
        return self.adata.obsm["X_totalVI"]


class Mammoth(Dataset):
    def __init__(self):
        super().__init__(
            filepath="../data/mammoth/mammoth.h5ad",
            name="Mammoth",
            hd_label_key="louvain",
            hd_metric="euclidean",
            description="Mammoth dataset from Wang et al. (2021)",
        )

    def get_HD_data(self):
        return self.adata.X

    def get_metadata_features(self):
        return self.adata.obs_keys()


class GaussLine(Dataset):
    def __init__(self):
        super().__init__(
            filepath="../data/gauss_line/gauss_line.h5ad",
            name="GaussLine",
            hd_label_key="labels_cat",
            hd_metric="euclidean",
            description="Gaussian clusters shifted along a line from Böhm et al. (2022)",
        )

    def get_HD_data(self):
        return self.adata.X

    def get_metadata_features(self):
        return self.adata.obs_keys()


datasets = {
    "Mammoth": Mammoth,
    "Mouse CD45neg": MouseCD45neg,
    #"Mouse Fibroblasts": MouseFibro,
    #"Human Immune": HumanImmune,
    "GaussLine": GaussLine,
}
