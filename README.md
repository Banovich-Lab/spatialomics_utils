# spatialomics_utils
Functions for spatial omics data analysis

## Includes:
### 1. Tissue-background segmentation
- Simply define a threshold to distinguish between background and tissue region. 
- Overview of all algorithms:

![Skimage Algorithm](img/skimage_algo.png)


### 2. Mask dilation
- Broaden the "mask" to the black pixels within the defined radius

![segment and dilate](img/demo.gif)

### 3. Convert Anndata <-> Rds
- Mitigate dependency on other pre-build packages
- Break anndata/rds object into small components (.mtx, .csv), then build up rds/anndata