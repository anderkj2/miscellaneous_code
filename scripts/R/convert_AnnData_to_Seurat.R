# Date created: 2/1/24
# convert_AnnData_to_Seurat.R - Creates a Seurat object from an AnnData object
# where data filtering and normalization have been performed, and possibly
# PCA and UMAP performed as well
# NOTE: output object will be generated in same directory as input
# NOTE: run script in py3 conda environment


require(optparse)
require(Seurat)
require(SeuratDisk)
require(anndataR)
require(dplyr)
require(pracma)


##### Input arguments #####
# Set command line parser arguments
option_list <- list(
  make_option(c("-d","--dir"), action="store", default=NULL, type="character",
              help="Directory containing input AnnData object file (without trailing '/')"),
  make_option(c("-i","--input"), action="store", default=NULL, type="character",
              help="File path for input AnnData object")
)

args <- parse_args(OptionParser(option_list=option_list))
dir <- args$dir
input <- args$input

# Set working directory
setwd(dir)
###########################


# Load AnnData object in memory
adata <- read_h5ad(input, to = "InMemoryAnnData")

# Convert AnnData object to Seurat object
sdata <- adata$to_Seurat()

# NOTE: anndataR may not properly transfer dimensionality reductions
# Check if any reductions exist in converted Seurat object
sreductions <- !isempty(sdata@reductions)


### If dimensionality reductions don't exist in converted Seurat object,
### manually transfer them from the AnnData object
if (sreductions == FALSE) {
  ## Transfer PCA dim reduction object
  # Check for existence of PCA embeddings in AnnData object
  if (!is.null(dim(adata$obsm$X_pca))) {
    # Extract PCA embeddings
    pc_embeddings <- adata$obsm$X_pca
    
    # Format row and column names
    rownames(pc_embeddings) <- rownames(sdata@meta.data)
    colnames(pc_embeddings) <- paste0("PC_",1:ncol(adata$obsm$X_pca))
    
    # Extract PCA loadings
    pc_loadings <- adata$varm$PCs
    
    # Format row and column names
    rownames(pc_loadings) <- rownames(sdata@assays$RNA@meta.features)
    colnames(pc_loadings) <- paste0("PC_",1:ncol(adata$obsm$X_pca))
    
    # Add PCA data to AnnData object
    sdata@reductions$pca <- CreateDimReducObject(
      embeddings = pc_embeddings,
      loadings = pc_loadings,
      key = "PC_",
      assay = "RNA"
    )
  }
  
  
  ## Transfer UMAP dim reduction object
  # Check for existence of UMAP embeddings in AnnData object
  if (!is.null(dim(adata$obsm$X_umap))) {
    # Extract UMAP embeddings
    umap_embeddings <- adata$obsm$X_umap
    
    # Format row and column names
    rownames(umap_embeddings) <- rownames(sdata@meta.data)
    colnames(umap_embeddings) <- paste0("UMAP_",1:ncol(adata$obsm$X_umap))
    
    # Add UMAP data to AnnData object
    sdata@reductions$umap <- CreateDimReducObject(
      embeddings = umap_embeddings,
      key = "UMAP_",
      assay = "RNA"
    )
  }
}
###


# Save updated Seurat object in H5Seurat format
SaveH5Seurat(sdata, gsub(".h5ad",".h5seurat",input), overwrite = TRUE)



