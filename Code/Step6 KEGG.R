# 加載必要的R包
library(clusterProfiler)
# 用於可以刪除我不要的列
library(dplyr)
library(KEGG.db)
library(org.Hs.eg.db)
# 載入ggplot2套件
library(ggplot2)

# Github token until Wed, Jan 17 2024
#Sys.setenv(GITHUB_PAT = "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")
# 安裝createKEGGdb包（用於創建KEGG.db的包）
#remotes::install_github("YuLab-SMU/createKEGGdb")
# 下載KEGG.db數據
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
createKEGGdb::create_kegg_db('hsa')


# confirm directory location
getwd()
# setting directory location
setwd("C:/Users/benson/Desktop/Research/")
# 讀取.csv文件
negative <- read.csv('negative_celltype_genes.csv', na.strings = c("None",""))

# 提取cell type列表
celltype <- negative$X

# 初始化一個空的列表來存儲結果
genes <- list()

# 對每個細胞類型進行操作
genes <- lapply(celltype, function(x) {
  temp <- negative[which(celltype == x), 2:ncol(negative[which(celltype == x),])]
  temp <- temp %>% select_if(~!any(is.na(.))) # 刪除包含NA的列
  temp <- unlist(temp) # 將temp轉換為向量
  return(temp)})

# 給列表的每個元素命名
names(genes) <- celltype

# 確認資料型態
class(genes)

# 對每種細胞類型的基因列表進行轉換
ids <- lapply(genes, function(gene_list) {
  bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
})


# 進行KEGG富集分析
kegg_results <- list()
for(i in names(ids)){
  kegg_results[[i]] <- enrichKEGG(gene         = ids[[i]]$ENTREZID,
                                  organism     = 'hsa',
                                  pvalueCutoff = 0.05)}

# 對每種細胞類型的結果進行繪圖
lapply(names(kegg_results), function(celltype) {
  # 提取結果
  result <- kegg_results[[celltype]]
  # 繪製點狀圖
  p = dotplot(result, showCategory = 10) + ggtitle(paste0(celltype, " KEGG enrichment"))
  # 指定保存位置
  filepath <- paste0("C:/Users/benson/Desktop/Research/Result_picture/GSE174188/KEGG/", celltype, "_KEGG_enrichment.png")
  # 保存圖片
  ggsave(filepath, plot = p)
})
