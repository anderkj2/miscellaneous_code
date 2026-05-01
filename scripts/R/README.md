# R Scripts

`convert_AnnData_to_Seurat.R` - Creates a Seurat object from an AnnData object where data filtering and normalization have been performed, and possibly PCA and UMAP performed as well, and transfers results of PCA and/or UMAP if present

`detect_AnnData_snRNA_doublets.R` - Converts an AnnData object containing raw counts with normalized cluster classifications and Scrublet-identified doublets in metadata to Seurat format and performs doublet detection using scDblFinder

