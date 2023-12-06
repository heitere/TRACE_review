# Interactively Exploring Embedding Quality with TRACE

Analyzing **global and local quality 🕵🏽‍♀️** of two-dimensional embeddings, based on [Regl-scatterplot](https://github.com/flekschas/regl-scatterplot).

## Installation

### Option 1: Using Docker 🐋

<details>
    <summary>Did you ran it without Docker before? 💡</summary>
    Make sure you change the <u>destination</u> in <i>frontend/next.config.js</i> to <code>destination: "http://backend:8000/backend/:path*"</code>.
</details>

```bash
docker-compose build
docker-compose up
```
This will mount the /frontend, /backend, and /data directories into the repective containers.
The /frontend/.next folder will be recreated every time the frontend is started. One might have to delete this folder before running the app without docker (as it is owned by 'root'). 

Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

### Option 2: Without Docker

#### Required packages
**Backend**: Install the required python packages for the backend, tested with Python 3.11 from `backend/pip_requirements.txt` or `backend/conda_requirements.yml`. 

**Frontend**: Install the packages in `frontend/package.json` using e.g. `npm install`.


#### Usage

<details>
    <summary>First time trying without Docker? 💡</summary>
        Make sure your user has write access to /frontend/.next or delete this folder. 
        Change the <u>destination</u> in <i>frontend/next.config.js</i> to <code>destination: "http://localhost:8000/backend/:path*"</code>.
</details>


First, start the backend within the right python evironment:
```bash
conda activate backend_env/
uvicorn main:app --reload
```

Then start the frontend development server:
```bash
npm run dev
# or
yarn dev
# or
pnpm dev
```
Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

## Data

The datasets are stored in `data` and should be preprocessed to the [Anndata](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) format and include the following fields:

* `adata.X` high-dimensional data
* `adata.obs`: metadata e.g. cluster labels
* `adata.obsm` low-dimensional embeddings, one entry for each embedding, e.g. `adata.obsm["t-SNE_0"]` for the first t-SNE embedding. 
* `adata.uns` unstructured data:

    * `adata.uns["methods"]`: a dictionary with the embedding methods and the number of embeddings for each method (i.e. for the slider to display an embedding sequence)
    * `adata.uns["neighbors"]`: an `nxk` array of the k-nearest high-dimensional neighbors of each point
    * `adata.uns["t-SNE_0"]`: a data frame with quality measures for each embedding

### Gaussian Line 🟢 🟠 🟣
A small example dataset that is included in the repository. UMAP and t-SNE are already precomputed, but the HD neighbors and quality measures are computed when the dashboard is loaded. This data will also be shown when the files for other dataset cannot be shown.

### Mammoth 🦣
This dataset is from Wang et al. and can be downloaded from their [PaCMAP](https://github.com/YingfanWang/PaCMAP/blob/master/data/mammoth_3d_50k.json) repository.

It then needs to be processed using the `mammoth.ipynb` notebook. 

### Single-Cell Data 🐁
The processed dataset of gene expressions from [Guilliams et al.](https://pubmed.ncbi.nlm.nih.gov/35021063/) is not available online, please reach out if you are interested. A raw version is available under [GSE192742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192742).


## References
This dashboard is based on the awesome [Regl-Scatterplot](https://github.com/flekschas/regl-scatterplot) library.
