import os
import subprocess
import sys


def _ensure_package(module_name, package_name=None):
    package_name = package_name or module_name
    try:
        __import__(module_name)
    except ImportError:
        print(f"Installing missing Python package: {package_name}")
        subprocess.check_call([
            sys.executable,
            "-m",
            "pip",
            "install",
            package_name,
        ])


_ensure_package("scipy")
_ensure_package("scanpy")
_ensure_package("pandas")

import scipy
import scanpy
import pandas as pd

from sys import argv

if len(argv) != 3:
    raise SystemExit("Usage: break_anndata.py YOUR_ADATA_PATH YOUR_OUTPUT_DIR")

ADATA_PATH = argv[1]
OUTPUT_DIR = argv[2]

def _write_metadata(adata, output_dir):
    print ("Writing the adata.obs")
    adata.obs.to_csv(
        os.path.join(output_dir, "obs.csv"), 
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
            columns=[f"{key}_{i + 1}" for i in range(adata.obsm[key].shape[1])]
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

    try:
        if adata.raw.X.shape == adata.X.shape:
            print ("Writing the adata.raw.X")
            scipy.io.mmwrite(
                os.path.join(output_dir, "layers", "raw_X.mtx"), 
                adata.raw.X.T
            )
    except AttributeError:
        print("No raw data found.")

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
    print ("Done!")
    return

print ('Reading AnnData object')
adata = scanpy.read_h5ad(ADATA_PATH)
break_anndata(adata, OUTPUT_DIR)