from typing import List
from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
import numpy as np
from data import datasets
from contextlib import asynccontextmanager
import pandas as pd
from pandas.api.types import CategoricalDtype
import continuous_palettes

dataset = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    print("Starting server")
    yield
    if dataset is not None:
        dataset.cleanup()
    print("Shutdown server")


app = FastAPI(title="Dashboard API", lifespan=lifespan)

origins = [
    "http://localhost:3000",
    "http://localhost:8000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class NeighborRequest(BaseModel):
    k: int
    points: List[int]
    hd_metric: str
    

@app.get("/backend/datasetOptions")
async def getDatasetNames():
    return {"result": list(datasets.keys())}


@app.get("/backend/loadDataset")
async def loadDataset(datasetName: str):
    """
    Load a dataset by name.

    Args:
        datasetName (str): The name of the dataset to load. 

    Returns:
        dict: A dictionary containing the result, hd_metric, and dataset_info 
            to be displayed in the dashboard.
    """
    global dataset

    if dataset is not None:
        dataset.cleanup()

    if datasetName in datasets.keys():
        try:
            dataset = datasets[datasetName]()
            print(f"Loaded dataset {datasetName}")
        except FileNotFoundError as e:
            print(f"Error loading dataset {datasetName}: {e}")
            dataset = datasets["GaussLine"]()
            datasetName = "GaussLine"
    else:
        print(f"Error loading dataset {datasetName}. Not in dataset options {list(datasets.keys())}")
        dataset = datasets["GaussLine"]()
    return {
        "dataset_name": datasetName,
        "hd_metric": dataset.hd_metric,
        "dataset_info": dataset.description,
    }


@app.get("/backend/metadataNames")
async def getMetadataNames():
    return {"result": dataset.get_metadata_features()}


@app.get("/backend/featureNames")
async def getFeatureNames():
    return {"result": dataset.adata.var.index.tolist()}


@app.get("/backend/qualityNames")
async def getQualityNames():
    return {"result": dataset.quality_features}


@app.get("/backend/pointColorOptions")
async def getPointColorOptions():
    options = dataset.get_metadata_features()
    options = options + dataset.quality_features + dataset.adata.var.index.tolist()
    return {"result": options}


@app.get("/backend/featureValues/{fname}")
async def getFeatureValues(fname: str):
    if fname in dataset.adata.var.index:
        expression = dataset.adata.obs_vector(fname).tolist()
        return {"result": expression}
    else:
        raise ValueError(f"feature {fname} not in variable names of dataset")


@app.get("/backend/embeddingOptions")
async def getEmbeddingOptions():
    return dataset.get_embedding_options()


@app.get("/backend/embedding")
async def getEmbedding(embName: str, embScale: int):
    """Retrieve embedding from anndata.obsm
    
    Args:
        embName (str): The name of the embedding.
        embScale (int): The scale of the embedding.
    
    Returns:
        dict: A dictionary containing the x and y coordinates of the embedding.
    """
    embNameScale = dataset.get_embedding_name(embName, embScale)
    if embNameScale in dataset.adata.obsm_keys():
        result = {
            "x": dataset.adata.obsm[embNameScale][:, 0].tolist(),
            "y": dataset.adata.obsm[embNameScale][:, 1].tolist(),
        }
    else:
        print(
            f"Embedding {embName} with scale {embScale} not found in obsm_keys {dataset.adata.obsm_keys()}"
        )
        result = {
            "x": np.zeros((dataset.adata.n_obs,), dtype=int).tolist(),
            "y": np.zeros((dataset.adata.n_obs,), dtype=int).tolist(),
        }
    return result


@app.get("/backend/pointColor/{fname}")
async def getPointColors(fname: str, embeddingName: str, embeddingIndex: int):
    """Retrieve values for metadata or feature, normalize,
    and compute colormaps.

    Args:
        fname (str): name of the feature to be found in the anndata.obs_keys,
            anndata.var.index, or anndata.obsm[embeddingName].colnames()
        embeddingName (str): name of embedding when retrieving
            embedding specific values (such as quality scores).
            Defaults to None.
        embeddingIndex (int): index of embedding (scale)

    Returns:
        dict: Dictionary with the following keys:
            "values": list,
            "colorMap": dict with labels as keys and colors as values
            "type": continuous | categorical,
    """
    name = dataset.get_embedding_name(embeddingName, embeddingIndex)

    # metadata
    if fname in dataset.adata.obs_keys():
        fvalues = dataset.adata.obs[fname]

        if dataset.adata.obs[fname].dtype.name != "category":
            ftype = "continuous"
            range = [float(fvalues.min()), float(fvalues.max())]
            encoded_fvalues = (fvalues - range[0]) / (range[1] - range[0])
            colors = continuous_palettes.palettes["viridis"]
            colorticks = np.linspace(
                start=range[0], stop=range[1], num=len(colors), endpoint=True
            )
            colorticks = [f"{x:.2}" for x in colorticks]
            colors = dict(zip(colorticks, colors))
        else:
            ftype = "categorical"
            cat_dtype = CategoricalDtype(
                categories=list(fvalues.value_counts(ascending=False, sort=True).keys())
            )
            encoded_fvalues = pd.Series(fvalues.values.tolist()).astype(cat_dtype)
            encoded_fvalues = encoded_fvalues.cat.codes
            colors = dataset.get_category_colors(fname)
            range = [0, 1]

    # features
    elif fname in dataset.adata.var.index:
        if "norm_genes" in dataset.adata.layers:
            layer = "norm_genes"
        else:
            layer = "X"
        fvalues = dataset.adata.obs_vector(fname, layer=layer)
        range = [fvalues.min().astype(float), fvalues.max().astype(float)]
        encoded_fvalues = (fvalues - range[0]) / (range[1] - range[0])
        colors = continuous_palettes.palettes["viridis"]
        colorticks = np.linspace(
            start=range[0], stop=range[1], num=len(colors), endpoint=True
        )
        colorticks = [f"{x:.2}" for x in colorticks]
        colors = dict(zip(colorticks, colors))
        ftype = "continuous"

    # quality score for specific dataset
    elif name in dataset.adata.uns_keys() and fname in dataset.adata.uns[name]:
        fvalues = dataset.adata.uns[name][fname]

        # if correlation scale from [-1, 1] to [0,1]
        if "corr" in fname:
            feature_range = [-1, 1]
        else:
            # other quality features are within [0,1] (neighborhood preservation)
            feature_range = [0, 1]

        # only scale if outside of [0,1] range
        if feature_range[0] < 0 or feature_range[1] > 1:
            encoded_fvalues = (fvalues - feature_range[0]) / (
                feature_range[1] - feature_range[0]
            )

        else:
            encoded_fvalues = fvalues

        colors = continuous_palettes.palettes["continuous_PpGn"]
        colorticks = np.linspace(
            start=feature_range[0],
            stop=feature_range[1],
            num=len(colors),
            endpoint=True,
        )
        colorticks = [f"{x:.2}" for x in colorticks]
        colors = dict(zip(colorticks, colors))
        ftype = "continuous"
        print(
            f"quality score {fname} has encoded range min {encoded_fvalues.min()} max {encoded_fvalues.max()} mean {encoded_fvalues.mean()}"
        )

    # fallback
    else:
        print(f"{fname} is not in obs_names {dataset.adata.obs_keys()}")
        fvalues = encoded_fvalues = np.zeros((dataset.adata.n_obs,), dtype=int)
        colors = {"none": "#444444"}
        ftype = "categorical"

    res = {
        "values": fvalues.tolist(),
        "encoded_values": encoded_fvalues.tolist(),
        "colorMap": colors,
        "type": ftype,
    }
    return res


@app.post("/backend/precomputeAllNeighbors")
async def precomputeAllNeighbors(maxK: int, hd_metric: str):
    dataset.precompute_HD_neighbors(maxK, hd_metric)
    return {"result": True}


@app.post("/backend/intrusions")
async def getIntrusions(item: NeighborRequest):
    """
    Get intrusions based on the given NeighborRequest. Intrusions are defined
    as points that are in the selection item.points, but not in the k nearest
    neighbors of these points.

    Args:
        item (NeighborRequest): The NeighborRequest object containing the parameters for intrusion detection.

    Returns:
        dict: A dictionary containing the computed results.
            - "result": A list of computed intrusions.
            - "binary": A list representing the binary representation of the intrusions.
    """
    neighbors = dataset.get_HD_neighbors(
        k=item.k, metric=item.hd_metric, indices=np.asarray(item.points)
    )
    neighbors = neighbors.flatten()
    binary = np.zeros((dataset.adata.n_obs,), dtype=int)
    binary[item.points] = 1
    binary[neighbors] = 0
    return {"result": binary.nonzero()[0].tolist(), "binary": binary.tolist()}


@app.post("/backend/computeHDNeighbors")
async def computeHDNeighbors(item: NeighborRequest):
    """
    Compute high-dimensional neighbors based on the given parameters.

    Args:
        item (NeighborRequest): The request object containing the parameters.

    Returns:
        dict: A dictionary containing the computed results.
            - "result": A list of computed neighbors.
            - "binary": A list representing the binary representation of the neighbors.
    """
    neighbors = dataset.get_HD_neighbors(
        k=item.k, metric=item.hd_metric, indices=np.asarray(item.points)
    )
    neighbors = neighbors.flatten()
    binary = np.zeros((dataset.adata.n_obs,), dtype=int)
    binary[neighbors] = 1
    return {"result": binary.nonzero()[0].tolist(), "binary": binary.tolist()}


@app.post("/backend/computeReverseHDNeighbors")
async def computeReverseHDNeighbors(item: NeighborRequest):
    """
    Compute reverse HD neighbors of a point selection. Reverse HD neighbors are
    defined as points that have at least one of the points in the selection as
    one of their k nearest neighbors.

    Args:
        item (NeighborRequest): The NeighborRequest object containing the parameters.

    Returns:
        dict: A dictionary containing the indices and a binary encoding of the
            reverse neighbors of the selection.
    """
    neighbors = dataset.get_HD_neighbors(
        k=item.k, metric=item.hd_metric, indices=np.arange(dataset.adata.n_obs)
    )
    query_points = np.asarray(item.points)

    # having the second array as the smaller one is faster
    if len(query_points) > item.k:
        bool_result = np.apply_along_axis(
            lambda a: np.isin(query_points, a, assume_unique=True).any(),
            1,
            neighbors,
        )
    else:
        bool_result = np.apply_along_axis(
            lambda a: np.isin(a, query_points, assume_unique=True).any(),
            1,
            neighbors,
        )

    binary = bool_result.astype(int)
    return {"result": binary.nonzero()[0].tolist(), "binary": binary.tolist()}
