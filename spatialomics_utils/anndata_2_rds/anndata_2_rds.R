library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
output_rds <- args[2]

read_meta <- function(input_dir) {
    print ("Reading obs metadata")
    obs <- read.csv(file.path(input_dir, "obs.csv"), row.names=1)
    return (obs)
}

read_var <- function(input_dir) {
    print ("Reading var metadata")
    var <- read.csv(file.path(input_dir, "var.csv"), row.names=1)
    return (var)
}

read_embedding <- function(input_dir) {
    print ("Reading embeddings")
    all_embeddings <- list()
    for (file in list.files(
            file.path(input_dir, "embeddings"), pattern="*.csv"
    )) {
        embedding_name <- gsub(".csv", "", file)
        embedding_data <- read.csv(
            file.path(input_dir, "embeddings", file), row.names=1
        )
        all_embeddings[[embedding_name]] <- as.matrix(embedding_data)
    }
    return(all_embeddings)
}

read_matrix <- function(input_dir) {
    print ("Reading matrices")
    all_matrices <- list()
    all_matrices[["X"]] <- readMM(file.path(input_dir, "X.mtx"))
    for (file in list.files(
        file.path(input_dir, "layers"), pattern="*.mtx"
    )) {
        layer_name <- gsub(".mtx", "", file)
        layer_data <- readMM(file.path(input_dir, "layers", file))
        all_matrices[[layer_name]] <- layer_data
    }
    return (all_matrices)
}

build_rds <- function(input_dir, output_rds) {
    matrices <- read_matrix(input_dir)
    embeddings <- read_embedding(input_dir)
    obs <- read_meta(input_dir)
    var <- read_var(input_dir)

    row.names(matrices[["X"]]) <- row.names(var)
    colnames(matrices[["X"]]) <- row.names(obs)

    seurat_obj <- CreateSeuratObject(
        counts = matrices[["X"]],
        meta.data = obs,
        assay = "RNA"
    )
    
    for (layer_name in names(matrices)) {
        if (layer_name != "X") {
            row.names(matrices[[layer_name]]) <- row.names(var)
            colnames(matrices[[layer_name]]) <- row.names(obs)
            seurat_obj[[layer_name]] <- CreateAssayObject(counts = matrices[[layer_name]])
        }
    }
    
    for (embedding_name in names(embeddings)) {
        seurat_obj[[embedding_name]] <- CreateDimReducObject(
            embeddings = embeddings[[embedding_name]],
            assay = "RNA"
        )
    }

    print ("Saving RDS file")
    saveRDS(seurat_obj, file = output_rds)
}

build_rds(input_dir, output_rds)