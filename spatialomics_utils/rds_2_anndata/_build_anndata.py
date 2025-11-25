import os
import scipy
import anndata
import pandas as pd

from sys import argv

INPUT_DIR = argv[1]
ADATA_PATH = argv[2]


def _read_metadata(input_dir):
    print ("Reading the metadata")
    return pd.read_csv(
        os.path.join(input_dir, "metadata.csv"), 
        index_col=0
    )


def _read_embedding(input_dir):
    print ("Reading the embeddings")
    embedding_dir = os.path.join(input_dir, "embeddings")
    embeddings = {}
    for file in os.listdir(embedding_dir):
        if file.endswith(".csv"):
            key = file[:-4]
            df = pd.read_csv(
                os.path.join(embedding_dir, file), 
                index_col=0
            )
            embeddings[key] = df.values
    return embeddings


def _read_var(input_dir):
    print ("Reading the var")
    return pd.read_csv(
        os.path.join(input_dir, "var.csv"), 
        index_col=0
    )


def _read_matrix(input_dir):
    print ("Reading the matrix")
    matrix_path = os.path.join(input_dir, "X.mtx")
    return scipy.io.mmread(matrix_path).T


def _add_layers(adata, input_dir):
    print ("Reading the layers")
    for file in os.listdir(input_dir):
        if file.endswith(".mtx") and file != "X.mtx":
            key = file[:-4]
            matrix_path = os.path.join(input_dir, file)
            adata.layers[key] = scipy.io.mmread(matrix_path).T
    return adata


def build_anndata(input_dir, adata_path):
    matrix = _read_matrix(input_dir)
    var = _read_var(input_dir)
    metadata = _read_metadata(input_dir)
    embedding = _read_embedding(input_dir)

    adata = anndata.AnnData(
        X=matrix, 
        obs=metadata, 
        var=var
    )
    for key, value in embedding.items():
        adata.obsm[key] = value

    adata = _add_layers(adata, input_dir)

    adata.write_h5ad(adata_path)
    
    return adata_path