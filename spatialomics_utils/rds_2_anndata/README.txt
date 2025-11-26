USAGE:

Step 1: Separate rds into components
Rscript break_rds.R YOUR_RDS_PATH YOUR_OUTPUT_FOL

Step 2: Create a new anndata object from your output_fol
python3 rds_2_anndata.py YOUR_OUTPUT_FOL OUT_H5AD_PATH