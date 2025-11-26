library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
INPUT_RDS <- args[1]
OUTPUT_DIR <- args[2]

write_features <- function(features_slot, output_dir) {
    print ("Writing the var metadata")
    write.csv(
        features_slot,
        file.path(output_dir, "var.csv"),
        row.names=FALSE,
        quote=FALSE
    )
    return
}

write_metadata <- function(meta_slot, output_dir) {
    print ("Writing the obs metadata")
    write.csv(
        meta_slot,
        file.path(output_dir, "obs.csv"),
        row.names=TRUE,
        quote=FALSE
    )
    return
}

write_matrix <- function(matrix_slot, output_dir) {
    print ("Writing the main matrix")
    Matrix::writeMM(
        matrix_slot,
        file.path(output_dir, "X.mtx")
    )
    return
}

write_embeddings <- function(seurat_obj, output_dir) {
    print ("Writing the embeddings")
    for (reduc_name in names(seurat_obj@reductions)) {
        embedding_data <- Embeddings(seurat_obj, reduction=reduc_name)
        write.csv(
            embedding_data,
            file.path(output_dir, paste0(reduc_name, ".csv")),
            row.names=TRUE,
            quote=FALSE
        )
    }
    return
}

write_assay5_matrix <- function(seurat_obj, assay_name, output_dir) {
    print (paste0("Writing the assay layers for assay: ", assay_name))
    for (s_assay in names(seurat_obj@assays[[assay_name]]@layers)) {
        s_name <- paste0(assay_name, '_', s_assay)
        Matrix::writeMM(
            seurat_obj@assays[[assay_name]]@layers[[s_assay]],
            file.path(output_dir, paste0(s_name, ".mtx"))
        )
    }
    return
}

write_assay_matrix <- function(seurat_obj, assay_name, output_dir) {
    print (paste0("Writing the assay layers for assay: ", assay_name))
    try (
        Matrix::writeMM(
            seurat_obj@assays[[assay_name]]@counts,
            file.path(
                output_dir, 
                paste0(assay_name, '_counts.mtx')
            )
        )
    )
    try (
        Matrix::writeMM(
            seurat_obj@assays[[assay_name]]@data,
            file.path(
                output_dir, 
                paste0(assay_name, '_data.mtx')
            )
        )
    )
    try (
        Matrix::writeMM(
            seurat_obj@assays[[assay_name]]@scale.data,
            file.path(
                output_dir, 
                paste0(assay_name, '_scale.data.mtx')
            )
        )
    )
    return
}

break_rds <- function(input_rds, output_dir) {
    dir.create(output_dir, showWarnings = FALSE)
    
    seurat_obj <- readRDS(input_rds)
    ## Write matrix of RNA counts
    a_name <- "RNA"
    a_class = class(seurat_obj@assays[[a_name]])[[1]]
    n_cells <- ncol(seurat_obj@assays[[a_name]])

    if (a_class == 'Assay5') {
        write_assay5_matrix(seurat_obj, a_name, output_dir)
    } else {
        write_assay_matrix(seurat_obj, a_name, output_dir)
    }
    write_features(
        row.names(seurat_obj@assays[[a_name]]), 
        output_dir
    )
    write_metadata(seurat_obj@meta.data, output_dir)

    ## Write all other matrices
    all_assays <- names(seurat_obj@assays)
    all_assays <- all_assays[all_assays != "RNA"]
    for (a_name in all_assays) {
        layers_fol <- file.path(output_dir, "layers")
        dir.create(layers_fol, showWarnings = FALSE)
        
        if (ncol(seurat_obj@assays[[a_name]]) != n_cells) {
            next
        }
        
        a_class = class(seurat_obj@assays[[a_name]])[[1]]
        if (a_class == 'Assay5') {
            write_assay5_matrix(seurat_obj, a_name, layers_fol)
        } else {
            write_assay_matrix(seurat_obj, a_name, layers_fol)
        }
    }

    ## Write embeddings
    embeddings_fol <- file.path(output_dir, "embeddings")
    dir.create(embeddings_fol, showWarnings = FALSE)
    write_embeddings(seurat_obj, embeddings_fol)

    return
}

break_rds(INPUT_RDS, OUTPUT_DIR)