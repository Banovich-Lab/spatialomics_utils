# $1 : input anndata .h5ad file
# $2 : output directory for broken down anndata components
# $3 : output RDS file path
# Usage: bash anndata_2_rds.sh input.h5ad output_dir output.rds

python3 _break_anndata.py $1 $2
Rscript _build_rds.R $2 $3