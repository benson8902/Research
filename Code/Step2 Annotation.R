library(anndata)
library(scater)
library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(pbmcref.SeuratData)

# confirm directory location
getwd()

# setting directory location
setwd("C:\\Users\\User\\Desktop\\論文\\Research\\Code\\Data\\Annotation")

# Convert anndata to h5seurat
Convert("C:\\Users\\User\\Desktop\\論文\\Research\\Code\\Data\\GSE165080_healthy.h5ad", dest = "h5seurat", overwrite = T)

# Load h5seurat
pbmc = LoadH5Seurat("C:\\Users\\User\\Desktop\\論文\\Research\\Code\\Data\\GSE165080_healthy.h5seurat", assays = "RNA", meta.data = FALSE, misc = FALSE)

# Cell type annotation
pbmc = RunAzimuth(pbmc, reference = "pbmcref")


pbmc = DietSeurat(pbmc)

# Save and convert to h5ad
SaveH5Seurat(pbmc, filename = "GSE165080_healthy_annotation.h5seurat", overwrite = T)
Convert("GSE165080_healthy_annotation.h5seurat", dest = "h5ad", overwrite = T)



