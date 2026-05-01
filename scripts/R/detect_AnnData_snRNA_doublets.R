# Date created: 6/15/24
# Date last modified (for public viewing): 5/1/26
# detect_AnnData_snRNA_doublets.R - Converts an AnnData object containing raw counts 
# with normalized cluster classifications and Scrublet-identified doublets in metadata 
# to Seurat format and performs doublet detection using scDblFinder
# NOTE: run script in py3 conda environment


require(optparse)
require(Seurat)
require(SeuratData)
require(SeuratDisk)
require(anndataR)
require(dplyr)
require(pracma)
#require(pipeComp)
#source(system.file("extdata", "scrna_alternatives.R", package="pipeComp"))
require(SingleCellExperiment)
#require(variancePartition)


##### Input arguments #####
# Set command line parser arguments
option_list <- list(
  make_option(c("-i","--indir"), action="store", default=NULL, type="character",
              help="Directory containing input AnnData object file (with trailing '/')"),
  make_option(c("-s","--sample"), action="store", default=NULL, type="character",
              help="Individual sample ID"),
  make_option(c("-c","--cohort"), action="store", default=NULL, type="character",
              help="Sample cohort ID"),
  make_option(c("-o","--outdir"), action="store", default=NULL, type="character",
              help="Output directory")
)

args <- parse_args(OptionParser(option_list=option_list))
indir <- args$indir
sample <- args$sample
cohort <- args$cohort
outdir <- args$outdir


# Set working directory
setwd(outdir)

# Set input files
h5ad_file <- paste0(indir,sample,"_snRNA_raw_expression-updated_Pair#,doublets_detected,cell_typed.h5ad")
meta_file <- paste0(indir,sample,"_snRNA_raw_expression-updated_Pair#,doublets_detected,cell_typed_metadata.csv")
features_file <- paste0(indir,sample,"_snRNA_raw_expression-updated_Pair#,doublets_detected,cell_typed_features.csv")
###########################


##### Load AnnData object and convert to SingleCellExperiment format #####
# Load AnnData object in memory
adata <- read_h5ad(h5ad_file, to = "InMemoryAnnData")

# Convert to a SingleCellExperiment object
sce <- adata$to_SingleCellExperiment()


### Reformat converted SCE object for consistency with scDblFinder functions
# Rename the default assay
names(assays(sce)) <- "counts"

# NOTE: The use of the default assay data type results in the error
#   no method or default for coercing “dgRMatrix” to “dgCMatrix”
# Need to change the counts matrix data type
sce@assays@data[[1]] <- as(sce@assays@data[[1]], 'CsparseMatrix')
###
##########################################################################


##### Run scDblFinder using the fast clustering approach, which computes a PCA and #####
##### KNN graph and generates artificial doublets by combining random cells from   #####
##### different clusters, proportionally to the cluster sizes                      #####
# Run sclDblFinder, first running the fast clustering approach
sce <- scDblFinder::scDblFinder(sce, clusters = TRUE)

# Export updated metadata with doublet detection results as a CSV table
write.csv(sce@colData,
          paste0(outdir,sample,"_snRNA_raw_expression-updated_Pair#,doublets_detectedx2,cell_typed_metadata.csv"),
          quote = FALSE, col.names = NA)
########################################################################################


##### Generate updated Seurat object and export for downstream analysis #####
# Convert updated SCE object to AnnData format 
newAdata <- from_SingleCellExperiment(sce, output_class = "InMemory")

# Write updated AnnData object to file
to_HDF5AnnData(newAdata, 
               paste0(outdir,sample,"_snRNA_raw_expression-updated_Pair#,doublets_detectedx2,cell_typed.h5ad"),
               compression = "none")


### Reformat AnnData object components for Seurat object generation
# Extract AnnData components for Seurat generation
acounts <- newAdata$X
ameta <- newAdata$obs

# Update the rownames to include barcode information rather than indices
rownames(acounts) <- ameta$index
rownames(ameta) <- ameta$index

# Transpose matrix so that genes are rows and cells are columns
acounts <- t(acounts)

# NOTE: Seurat object can't accept underscores in input
# Convert underscores to dashes
rownames(acounts) <- gsub("__","_",rownames(acounts))
rownames(acounts) <- gsub("_","-",rownames(acounts))

# NOTE: Feature reformatting results in duplicates, which Seurat can't accept
# Fix duplicated features
if (any(duplicated(rownames(acounts)))) {
  # Identify features with duplicates
  duplicates <- rownames(acounts)[which(duplicated(rownames(acounts)))]
  
  # Loop through duplicates and make second entry unique
  for (i in 1:length(duplicates)) {
    # Find all features with same name
    locs <- grep(duplicates[i],rownames(acounts))
    
    # Append a ".1" to second instance of feature
    rownames(acounts)[grep(duplicates[i],rownames(acounts))[2]] <- 
      paste0(rownames(acounts)[grep(duplicates[i],rownames(acounts))[2]],".1")
  }
}
###


# Manually create a Seurat object from AnnData components
sdata <- CreateSeuratObject(
  acounts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = ameta,
  project = cohort,
  min.cells = 0,
  min.features = 0
)

# Save Seurat object in H5Seurat format
SaveH5Seurat(sdata,
             paste0(outdir,sample,"_snRNA_raw_expression-updated_Pair#,doublets_detectedx2,cell_typed.h5Seurat"),
             overwrite = TRUE)

# NOTE: counts are stored in sdata@assays$RNA@counts
#############################################################################


