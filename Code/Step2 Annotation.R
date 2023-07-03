library(anndata)
library(scater)
library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(pbmcref.SeuratData)

# confirm directory location
getwd()

# setting directory location
setwd("C:\\Users\\User\\Desktop\\論文\\PBMC_Data\\HCA\\Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19")

# Convert anndata to h5seurat
Convert("benson_pbmc_healthy_donor.h5ad", dest = "h5seurat", overwrite = T)

# Load h5seurat
pbmc = LoadH5Seurat("benson_pbmc_healthy_donor.h5seurat", assays = "RNA", meta.data = FALSE, misc = FALSE)

# Cell type annotation
pbmc = RunAzimuth(pbmc, reference = "pbmcref")


pbmc = DietSeurat(pbmc)

# Save and convert to h5ad
SaveH5Seurat(pbmc, filename = "benson_pbmc_post_healthy_donor.h5seurat", overwrite = T)
Convert("benson_pbmc_post_healthy_donor.h5seurat", dest = "h5ad", overwrite = T)



