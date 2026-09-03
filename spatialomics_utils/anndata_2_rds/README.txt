USAGE:

One-command conversion:
bash anndata_2_rds.sh YOUR_ADATA_PATH OUT_RDS_PATH

#######
Internal steps:

Step 1: Separate AnnData into components
python3 break_anndata.py YOUR_ADATA_PATH YOUR_OUTPUT_FOL

Step 2: Create a new Seurat object from your output folder
Rscript build_rds.R YOUR_OUTPUT_FOL OUT_RDS_PATH