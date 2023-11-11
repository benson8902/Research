# 加載必要的R包
library(clusterProfiler)
# 用於可以刪除我不要的列
library(dplyr)
library(org.Hs.eg.db)
library(GO.db)
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
negative <- read.csv('pearson_negative_celltype_genes_top200.csv', na.strings = c("None",""))
# 讀取.csv文件
positive <- read.csv('pearson_positive_celltype_genes_top200.csv', na.strings = c("None",""))




# Negative
# 提取cell type列表
celltype_N <- negative$X
# 提取cell type列表
celltype_P <- positive$X
# 初始化一個空的列表來存儲結果
genes_N <- list()
# 初始化一個空的列表來存儲結果
genes_P <- list()

# 對每個細胞類型進行操作
genes_N <- lapply(celltype_N, function(x) {
  temp <- negative[which(celltype_N == x), 2:ncol(negative[which(celltype_N == x),])]
  temp <- temp %>% select_if(~!any(is.na(.))) # 刪除包含NA的列
  temp <- unlist(temp) # 將temp轉換為向量
  return(temp)})
# 對每個細胞類型進行操作
genes_P <- lapply(celltype_P, function(x) {
  temp <- positive[which(celltype_P == x), 2:ncol(positive[which(celltype_P == x),])]
  temp <- temp %>% select_if(~!any(is.na(.))) # 刪除包含NA的列
  temp <- unlist(temp) # 將temp轉換為向量
  return(temp)})

# 給列表的每個元素命名
names(genes_N) <- celltype_N
# 給列表的每個元素命名
names(genes_P) <- celltype_P

# 確認資料型態
class(genes_N)
# 確認資料型態
class(genes_P)

# 對每種細胞類型的基因列表進行轉換
ids_N <- lapply(genes_N, function(gene_list) {
  bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)})
# 對每種細胞類型的基因列表進行轉換
ids_P <- lapply(genes_P, function(gene_list) {
  bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)})

# 進行GO(Gene Ontology) analysis
go_results_N <- list()
for(i in names(ids_N)){
  go_results_N[[i]] <- enrichGO(gene         = ids_N[[i]]$ENTREZID,
                              OrgDb         = org.Hs.eg.db,
                              ont           = "BP",  # 進行 Biological Process 的分析
                              pAdjustMethod = "BH",  # 使用 Benjamini & Hochberg (BH) 方法來調整 p 值
                              pvalueCutoff  = 0.05,  # 設定 p 值的閾值為 0.05
                              qvalueCutoff  = 0.05,  # 設定 q 值的閾值為 0.05
                              readable      = TRUE)}  # 將 ENTREZ ID 轉換為可讀的基因名稱
# 進行GO(Gene Ontology) analysis
go_results_P <- list()
for(i in names(ids_P)){
  go_results_P[[i]] <- enrichGO(gene         = ids_P[[i]]$ENTREZID,
                              OrgDb         = org.Hs.eg.db,
                              ont           = "BP",  # 進行 Biological Process 的分析
                              pAdjustMethod = "BH",  # 使用 Benjamini & Hochberg (BH) 方法來調整 p 值
                              pvalueCutoff  = 0.05,  # 設定 p 值的閾值為 0.05
                              qvalueCutoff  = 0.05,  # 設定 q 值的閾值為 0.05
                              readable      = TRUE)}  # 將 ENTREZ ID 轉換為可讀的基因名稱

# 進行KEGG Enrichment analysis
kegg_results_N <- list()
for(i in names(ids_N)){
  kegg_results_N[[i]] <- enrichKEGG(gene         = ids_N[[i]]$ENTREZID,
                                  organism     = 'hsa',
                                  pvalueCutoff = 0.05)}
# 進行KEGG Enrichment analysis
kegg_results_P <- list()
for(i in names(ids_P)){
  kegg_results_P[[i]] <- enrichKEGG(gene         = ids_P[[i]]$ENTREZID,
                                  organism     = 'hsa',
                                  pvalueCutoff = 0.05)}

# 定義一個函數來繪製並保存結果
plot_results_N <- function(results, result_type) {
  lapply(names(results), function(celltype_N) {
    # 提取結果
    result <- results[[celltype_N]]
    # 檢查結果是否有效
    if (is.null(result) || nrow(result) == 0) {
      warning(paste("No valid", result_type, "results for", celltype_N))
      return()
    }
    # 繪製點狀圖
    p = dotplot(result, showCategory = 10) + ggtitle(paste0(celltype_N, " Negative ", result_type, " enrichment"))
    # 指定保存位置
    filepath <- paste0("C:/Users/benson/Desktop/Research/Result_picture/GSE174188/Pearson/Top200/", result_type, "/Negative/", celltype_N, "_negative_", result_type, "_enrichment.png")
    # 保存圖片
    ggsave(filepath, plot = p)})}

# 定義一個函數來繪製並保存結果
plot_results_P <- function(results, result_type) {
  lapply(names(results), function(celltype_P) {
    # 提取結果
    result <- results[[celltype_P]]
    # 檢查結果是否有效
    if (is.null(result) || nrow(result) == 0) {
      warning(paste("No valid", result_type, "results for", celltype_N))
      return()
    }
    # 繪製點狀圖
    p = dotplot(result, showCategory = 10) + ggtitle(paste0(celltype_P, " Positive ", result_type, " enrichment"))
    # 指定保存位置
    filepath <- paste0("C:/Users/benson/Desktop/Research/Result_picture/GSE174188/Pearson/Top200/", result_type, "/Positive/", celltype_P, "_positive_", result_type, "_enrichment.png")
    # 保存圖片
    ggsave(filepath, plot = p)})}



# 對 GO 結果進行繪圖
plot_results_N(go_results_N, "GO")
# 對 KEGG 結果進行繪圖
plot_results_N(kegg_results_N, "KEGG")

# 對 GO 結果進行繪圖
plot_results_P(go_results_P, "GO")
# 對 KEGG 結果進行繪圖
plot_results_P(kegg_results_P, "KEGG")







#ID：GO term的ID。
#Description：GO term的描述，說明這個GO term代表的生物學過程、細胞組分或分子功能。
#GeneRatio：在你的基因集中，有多少基因和這個GO term相關。
#BgRatio：在背景基因集中，有多少基因和這個GO term相關。
#pvalue：p值，用來評估這個GO term是否顯著富集在你的基因集中。
#p.adjust：經過多重檢定校正後的p值。
#qvalue：q值，另一種經過多重檢定校正後的p值。
#Count：在你的基因集中，有多少基因和這個GO term相關。
#取出GO term的部分
head(go_results_N)
# 初始化一個空的向量來存儲結果
go_term_Description_N <- c()
# 將每個cell type pvalue最小的10個go term合併在一起
for(celltype in names(go_results_N)){
  #從enrichResult物件中提取結果
  result <- go_results_N[[celltype]]@result
  # 根據p值對GO term進行排序，並取出p值最小的前10個
  go_term <- result[order(result$pvalue, decreasing = FALSE), ][1:10, ]
  # 將description挑出來
  descriptions <- go_term$Description
  go_term_Description_N <- c(go_term_Description_N, descriptions)
}
# 把重複的description移除
go_term_Description_N <- unique(go_term_Description_N)

# 初始化一個空的向量來存儲結果
go_term_Description_P <- c()
# 將每個cell type pvalue最小的10個go term合併在一起
for(celltype in names(go_results_P)){
  #從enrichResult物件中提取結果
  result <- go_results_P[[celltype]]@result
  # 根據p值對GO term進行排序，並取出p值最小的前10個
  go_term <- result[order(result$pvalue, decreasing = FALSE), ][1:10, ]
  # 將description挑出來
  descriptions <- go_term$Description
  go_term_Description_P <- c(go_term_Description_P, descriptions)
}
# 把重複的description移除
go_term_Description_P <- unique(go_term_Description_P)

# 合併Negative & Positive 的 Description
go_term_Description <- unique(c(go_term_Description_N, go_term_Description_N))



# 建立一個list要放N個cell type的GeneRatio
go_term_GeneRatio_N <- list()
# 把這個description與go_result match的GeneRatio挑出來
for(celltype in names(go_results_N)){
  result <- go_results_N[[celltype]]@result
  go_term_GeneRatio_N[[celltype]] <- result$GeneRatio[match(go_term_Description,result$Description)]
  # 把分數轉成小數點
  go_term_GeneRatio_N[[celltype]] <- as.numeric(sub("/.*", "", go_term_GeneRatio_N[[celltype]])) / as.numeric(sub(".*/", "", go_term_GeneRatio_N[[celltype]]))
  # 把NA轉成0
  go_term_GeneRatio_N[[celltype]][is.na(go_term_GeneRatio_N[[celltype]])] <- 0
}
# 建立一個list要放N個cell type的GeneRatio
go_term_GeneRatio_P <- list()
# 把這個description與go_result match的GeneRatio挑出來
for(celltype in names(go_results_P)){
  result <- go_results_P[[celltype]]@result
  go_term_GeneRatio_P[[celltype]] <- result$GeneRatio[match(go_term_Description,result$Description)]
  # 把分數轉成小數點
  go_term_GeneRatio_P[[celltype]] <- as.numeric(sub("/.*", "", go_term_GeneRatio_P[[celltype]])) / as.numeric(sub(".*/", "", go_term_GeneRatio_P[[celltype]]))
  # 把NA轉成0
  go_term_GeneRatio_P[[celltype]][is.na(go_term_GeneRatio_P[[celltype]])] <- 0
}

# 合併Negative & Positive 的 GeneRatio
go_term_GeneRatio <- c(go_term_GeneRatio_N, go_term_GeneRatio_P)


# 開始做heatmap的處理
matrix_GO <- do.call(rbind, go_term_GeneRatio)
class(matrix_GO)
colnames(matrix_GO) <- go_term_Description
# 畫heatmap
pheatmap(t(matrix_GO), cluster_rows = FALSE, cluster_cols = FALSE)


matrix_GO
