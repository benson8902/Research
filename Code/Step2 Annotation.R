library(anndata)
library(scater)
library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(pbmcref.SeuratData)

# confirm directory location
getwd()

# setting directory location
setwd("C:/Users/benson/Desktop/Research/Annotation")

# Convert anndata to h5seurat
Convert("C:/Users/benson/Desktop/Research/GSE165080_healthy.h5ad", dest = "h5seurat", overwrite = T)

# Load h5seurat
pbmc = LoadH5Seurat("C:/Users/benson/Desktop/Research/GSE165080_healthy.h5seurat", assays = "RNA", meta.data = FALSE, misc = FALSE)

# Cell type annotation
pbmc = RunAzimuth(pbmc, reference = "pbmcref")


pbmc = DietSeurat(pbmc)

# Save and convert to h5ad
SaveH5Seurat(pbmc, filename = "GSE174188_healthy_annotation.h5seurat", overwrite = T)
Convert("GSE174188_healthy_annotation.h5seurat", dest = "h5ad", overwrite = T)





#如果package有問題就必須執行這個去安裝
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("scater")
#BiocManager::install("pbmcref.SeuratData")




# 確認套件是否被安裝
#if ("pbmcref.SeuratData" %in% rownames(installed.packages())) {
#  print("pbmcref.SeuratData 套件已經安裝")
#} else {
#  print("pbmcref.SeuratData 套件尚未安裝")
#}
# 還是裝不了就
# 如果Azimuth已經安裝但是SHA1沒有做改變，想要強制重新安裝，可以用這個
# devtools::install_github("satijalab/seurat-data", force = TRUE)
# 之後用
# library(SeuratData)
# AvaliableData()
# 看pbmcref.SeuratData版本多少
#install.packages("https://seurat.nygenome.org/src/contrib/pbmcref.SeuratData_1.0.0.tar.gz", repos = NULL, type = "source")
#如果不行就
#install.packages("/路徑/pbmcref.SeuratData_1.0.0.tar.gz", repos = NULL, type = "source")



# 更改library路徑，因為pbmcref.SeuratData會因為權限問題而無法安裝
# .libPaths("C:/Users/benson/R/win-library/")
