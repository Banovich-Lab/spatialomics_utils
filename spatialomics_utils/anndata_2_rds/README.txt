USAGE:

Step 1: Separate anndata into components
python3 break_adata.py YOUR_ADATA_PATH YOUR_OUTPUT_FOL

Step 2: Create a new seurat object from your output_fol
Rscript anndata_2_rds.R YOUR_OUTPUT_FOL OUT_RDS_PATH