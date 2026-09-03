USAGE:

###### One-command conversion: ###### 
bash rds_2_anndata.sh YOUR_RDS_PATH OUT_H5AD_PATH

####### Internal steps: ###### 
# Step 1: Separate rds into components
Rscript break_rds.R YOUR_RDS_PATH YOUR_OUTPUT_FOL

# Step 2: Create a new anndata object from your output_fol
python3 build_anndata.py YOUR_OUTPUT_FOL OUT_H5AD_PATH
########

Required packages are installed automatically the first time they are missing:
Python: scipy, anndata, pandas
R: Seurat, Matrix