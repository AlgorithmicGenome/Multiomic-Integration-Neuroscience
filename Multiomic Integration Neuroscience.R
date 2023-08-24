
# Set Directory
setwd("C:/Users/pgand/OneDrive/Documents/Bioinformatics Journey/Multiomic Integration Neuroscience Application/")

# Install Packages

install.packages("Signac")

# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("EnsDb.Hsapiens.v86")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biovizBase")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnsDb.Mmusculus.v79")

install.packages("hdf5r")

# Install Libraries
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(cowplot)
library(hdf5r)
library(ggplot2)

# the 10x hdf5 file contains both data types. 
alzheimers_obj <- Read10X_h5(filename = 'Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_filtered_feature_bc_matrix.h5')
# str(alzheimers_obj)
frag.file <- "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_atac_fragments.tsv.gz"

# 1. Extract RNA and ATAC data
rna <- alzheimers_obj$`Gene Expression`
atac <- alzheimers_obj$Peaks

# 2. Get gene annotations for mm10. Add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac <- atac[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# 3. Create Seurat object
mca <- CreateSeuratObject(counts = rna, assay = "RNA")

# 4. Create Seurat object and add the assays
mca[["ATAC"]] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
mca
# Check to see if files are loaded correctly
# file.exists(frag.file)
# file.exists(paste0(frag.file, '.tbi'))


# 5. Quality Control
DefaultAssay(mca) <- "ATAC"
mca <- NucleosomeSignal(mca)
mca <- TSSEnrichment(mca)
mca$blacklist_fraction <- FractionCountsInRegion(
  object = mca,
  assay = 'ATAC',
  regions = blacklist_mm10
)

Idents(mca) <- "all" # group all cells together, rather than by replicate
VlnPlot(mca,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "blacklist_fraction"),
  ncol = 5,
  log = TRUE,
  pt.size = 0) + NoLegend()

# 6. Filter out low quality cells
mca <- subset(
  x = mca,
  subset = nCount_ATAC < 7e4 &
    TSS.enrichment > 1 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    blacklist_fraction < 1
)
mca

# 7. DNA Accessibility Data Processing
# RNA analysis
DefaultAssay(mca) <- "RNA"
mca <- FindVariableFeatures(mca, nfeatures = 3000)
mca <- NormalizeData(mca)
mca <- ScaleData(mca)
mca <- RunPCA(mca, npcs = 30)
mca <- RunUMAP(mca, dims = 1:30, reduction.name = "umap.rna")
mca <- FindNeighbors(mca, dims = 1:30)
mca <- FindClusters(mca, resolution = 0.5, algorithm = 3)

# Plot 1
p1 <- DimPlot(mca, label = TRUE) + NoLegend() + ggtitle("RNA UMAP")
# mca <- SCTransform(mca, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# 8. ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mca) <- "ATAC"
mca <- RunTFIDF(mca)
mca <- FindTopFeatures(mca, min.cutoff = '10')
mca <- RunSVD(mca)
mca <- RunUMAP(mca, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Plot 2
p2 <- DimPlot(mca, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")

# Plot 1 + Plot 2
p1 + p2

# 9. Integration with scRNA-seq data
# label transfer from Allen brain
allen_rna <- readRDS("allen_brain.rds")
allen_rna
allen_rna <- FindVariableFeatures(
  object = allen_rna,
  nfeatures = 5000)

# 10.use the RNA assay in the mca-seq data for integration with scRNA-seq
DefaultAssay(mca) <- 'RNA'
transfer.anchors <- FindTransferAnchors(
  reference = allen_rna,
  query = mca,
  dims = 1:50,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_rna$subclass,
  weight.reduction = mca[['lsi']],
  dims = 2:50
)

mca <- AddMetaData(object = mca, metadata = predicted.labels)

# 11. Visualize Clustering Based on gene expression, ATAC-seq, or WNN Analysis
plot1 <- DimPlot(allen_rna, group.by = 'subclass', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(mca, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
p1 + p2

# 12 Replace each label with its most likely prediction
# replace each label with its most likely prediction
for(i in levels(mca)) {
  cells_to_reid <- WhichCells(mca, idents = i)
  newid <- names(which.max(table(mca$predicted.id[cells_to_reid])))
  Idents(mca, cells = cells_to_reid) <- newid
}
