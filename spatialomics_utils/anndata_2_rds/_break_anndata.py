import os
import scipy
import scanpy
import pandas as pd

from sys import argv

ADATA_PATH = argv[1]
OUTPUT_DIR = argv[2]

def _write_metadata(adata, output_dir):
    print ("Writing the adata.obs")
    adata.obs.to_csv(
        os.path.join(output_dir, "metadata.csv"), 
        index=True
    )
    return

def _write_embedding(adata, output_dir):
    os.makedirs(os.path.join(output_dir, "embeddings"), exist_ok=True)
    for key in adata.obsm.keys():
        print (f"Writing the adata.obsm[{key}]")
        df = pd.DataFrame(
            adata.obsm[key], 
            index=adata.obs_names, 
            columns=[f"{key}_{i}" for i in range(adata.obsm[key].shape[1])]
        )
        df.to_csv(
            os.path.join(output_dir, "embeddings", f"{key}.csv"), 
            index=True
        )
    return

def _write_var(adata, output_dir):
    print ("Writing the adata.var")
    adata.var.to_csv(
        os.path.join(output_dir, "var.csv"), 
        index=True
    )
    return

def _write_matrix(adata, output_dir):
    print ("Writing the adata.X")
    os.makedirs(
        os.path.join(output_dir, "layers"),
        exist_ok=True   
    )
    scipy.io.mmwrite(
        os.path.join(output_dir, "X.mtx"), 
        adata.X.T
    )
    for layer in adata.layers.keys():
        print (f"Writing the adata.layers[{layer}]")
        scipy.io.mmwrite(
            os.path.join(output_dir, "layers", f"{layer}.mtx"), 
            adata.layers[layer].T
        )
    return

def break_anndata(adata, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    _write_matrix(adata, output_dir)
    _write_metadata(adata, output_dir)
    _write_embedding(adata, output_dir)
    _write_var(adata, output_dir)
    return

adata = scanpy.read_h5ad(ADATA_PATH)
break_anndata(adata, OUTPUT_DIR)