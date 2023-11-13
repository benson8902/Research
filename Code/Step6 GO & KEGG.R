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


#========================================================Function============================================================
# 對每個細胞類型進行操作
get_genes <- function(celltype, data) {
  genes <- lapply(celltype, function(x) {
    temp <- data[which(celltype == x), 2:ncol(data[which(celltype == x),])]
    temp <- temp %>% select_if(~!any(is.na(.))) # 刪除包含NA的列
    temp <- unlist(temp) # 將temp轉換為向量
    return(temp)
  })
  return(genes)
}

# 對每種cell type的gene list進行轉換id
convert_gene_ids <- function(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) {
  ids <- lapply(genes, function(gene_list) {
    bitr(gene_list, fromType = fromType, toType = toType, OrgDb = OrgDb)
  })
  return(ids)
}

# 進行GO(Gene Ontology) analysis
perform_GO_analysis <- function(ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE) {
  go_results <- list()
  for(i in names(ids)){
    go_results[[i]] <- enrichGO(gene         = ids[[i]]$ENTREZID,
                                OrgDb        = OrgDb,
                                ont          = ont,  
                                pAdjustMethod= pAdjustMethod,  
                                pvalueCutoff = pvalueCutoff,  
                                qvalueCutoff = qvalueCutoff,  
                                readable     = readable)
  }
  return(go_results)
}

# 進行KEGG Enrichment analysis
perform_KEGG_analysis <- function(ids, organism = 'hsa', pvalueCutoff = 0.05){
  kegg_result <- list()
  for(i in names(ids)){
    kegg_result[[i]] <- enrichKEGG(gene         = ids[[i]]$ENTREZID,
                                   organism     = organism,
                                   pvalueCutoff = pvalueCutoff)}
  return(kegg_result)
}

# 定義一個function來繪製並保存結果
plot_results <- function(results, result_type, type) {
  lapply(names(results), function(celltype) {
    # 提取結果
    result <- results[[celltype]]
    # 檢查結果是否有效
    if (is.null(result) || nrow(result) == 0) {
      warning(paste("No valid", result_type, "results for", celltype))
      return()
    }
    # 繪製dotplot
    p = dotplot(result, showCategory = 10, font.size = 10) + ggtitle(paste0(celltype, " ", type, " ", result_type, " enrichment"))
    # 指定保存位置
    filepath <- paste0("C:/Users/benson/Desktop/Research/Result_picture/GSE174188/Pearson/Top200/", result_type, "/", type, "/", celltype, "_", tolower(type), "_", result_type, "_enrichment.png")
    # 保存圖片
    ggsave(filepath, plot = p, width = 5, height = 6)
  })
}


# GO & KEGG enrichment result
    #ID：GO term的ID。
    #Description：GO term的描述，說明這個GO term代表的生物學過程、細胞組分或分子功能。
    #GeneRatio：在你的基因集中，有多少基因和這個GO term相關。
    #BgRatio：在背景基因集中，有多少基因和這個GO term相關。
    #pvalue：p值，用來評估這個GO term是否顯著富集在你的基因集中。
    #p.adjust：經過多重檢定校正後的p值。
    #qvalue：q值，另一種經過多重檢定校正後的p值。
    #Count：在你的基因集中，有多少基因和這個GO term相關。

# 將top200 gene的term提取出來取聯集看一下GeneRatio的部分抓出來去做heatmap看一下negative和positive之間有沒有可證明性
create_heatmap <- function(result_N, result_P){
  # 初始化一個空的向量來存儲結果
  term_Description_N <- c()
  # 將每個cell type pvalue最小的10個go term合併在一起
  for(celltype in names(result_N)){
    #從enrichResult物件中提取結果
    result <- result_N[[celltype]]@result
    # 根據p值對GO term進行排序，並取出p值最小的前10個
    term <- result[order(result$pvalue, decreasing = FALSE), ][1:10, ]
    # 將description挑出來
    descriptions <- term$Description
    term_Description_N <- c(term_Description_N, descriptions)
  }
  # 把重複的description移除
  term_Description_N <- unique(term_Description_N)
  
  # 初始化一個空的向量來存儲結果
  term_Description_P <- c()
  # 將每個cell type pvalue最小的10個go term合併在一起
  for(celltype in names(result_P)){
    #從enrichResult物件中提取結果
    result <- result_P[[celltype]]@result
    # 根據p值對GO term進行排序，並取出p值最小的前10個
    term <- result[order(result$pvalue, decreasing = FALSE), ][1:10, ]
    # 將description挑出來
    descriptions <- term$Description
    term_Description_P <- c(term_Description_P, descriptions)
  }
  # 把重複的description移除
  term_Description_P <- unique(term_Description_P)
  
  # 合併Negative & Positive 的 Description
  term_Description <- unique(c(term_Description_N, term_Description_P))
  
  # 建立一個list要放N個cell type的GeneRatio
  term_GeneRatio_N <- list()
  
  # 把這個description與go_result match的GeneRatio挑出來
  for(celltype in names(result_N)){
    result <- result_N[[celltype]]@result
    term_GeneRatio_N[[celltype]] <- result$GeneRatio[match(term_Description,result$Description)]
    # 把分數轉成小數點
    term_GeneRatio_N[[celltype]] <- as.numeric(sub("/.*", "", term_GeneRatio_N[[celltype]])) / as.numeric(sub(".*/", "", term_GeneRatio_N[[celltype]]))
    # 把NA轉成0
    term_GeneRatio_N[[celltype]][is.na(term_GeneRatio_N[[celltype]])] <- 0
  }
  
  # 建立一個list要放N個cell type的GeneRatio
  term_GeneRatio_P <- list()
  # 把這個description與go_result match的GeneRatio挑出來
  for(celltype in names(result_P)){
    result <- result_P[[celltype]]@result
    term_GeneRatio_P[[celltype]] <- result$GeneRatio[match(term_Description,result$Description)]
    # 把分數轉成小數點
    term_GeneRatio_P[[celltype]] <- as.numeric(sub("/.*", "", term_GeneRatio_P[[celltype]])) / as.numeric(sub(".*/", "", term_GeneRatio_P[[celltype]]))
    # 把NA轉成0
    term_GeneRatio_P[[celltype]][is.na(term_GeneRatio_P[[celltype]])] <- 0
  }
  
  # 合併Negative & Positive 的 GeneRatio
  term_GeneRatio <- c(term_GeneRatio_N, term_GeneRatio_P)
  
  # 開始做heatmap的處理
  matrix_heatmap <- do.call(rbind, term_GeneRatio)
  class(matrix_heatmap)
  colnames(matrix_heatmap) <- term_Description
  # 畫heatmap
  pheatmap(t(matrix_heatmap), cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 8, angle_col = 45)
}




#============================================================Main============================================================
# confirm directory location
getwd()
# setting directory location
setwd("C:/Users/benson/Desktop/Research/")


# 讀取.csv文件
negative <- read.csv('pearson_negative_celltype_genes_top200.csv', na.strings = c("None",""))
# 讀取.csv文件
positive <- read.csv('pearson_positive_celltype_genes_top200.csv', na.strings = c("None",""))


# 提取cell type列表
celltype_N <- negative$X
# 提取cell type列表
celltype_P <- positive$X
# 初始化一個空的列表來存儲結果
genes_N <- list()
# 初始化一個空的列表來存儲結果
genes_P <- list()


# 取得每個cell type的gene
genes_N <- get_genes(celltype_N, negative)
genes_P <- get_genes(celltype_N, positive)


# list的index放上 cell type
names(genes_N) <- celltype_N
# list的index放上 cell type
names(genes_P) <- celltype_P


# 確認資料型態
class(genes_N)


# 對每種細胞類型的基因列表進行轉換
ids_N <- convert_gene_ids(genes_N)
# 對每種細胞類型的基因列表進行轉換
ids_P <- convert_gene_ids(genes_P)


# 進行GO(Gene Ontology) analysis
go_results_N <- perform_GO_analysis(ids_N)
# 進行GO(Gene Ontology) analysis
go_results_P <- perform_GO_analysis(ids_P)


# 進行KEGG Enrichment analysis
kegg_results_N <- perform_KEGG_analysis(ids_N)
# 進行KEGG Enrichment analysis
kegg_results_P <- perform_KEGG_analysis(ids_P)


# 對 GO 結果進行繪圖
plot_results(go_results_N, "GO", "Negative")
# 對 GO 結果進行繪圖
plot_results(go_results_P, "GO", "Positive")


# 對 KEGG 結果進行繪圖
plot_results(kegg_results_N, "KEGG", "Negative")
# 對 KEGG 結果進行繪圖
plot_results(kegg_results_P, "KEGG", "Positive")


# GO heatmap
create_heatmap(go_results_N, go_results_P) # 先放negative,再放 positive
# KEGG heatmap
create_heatmap(kegg_results_N, kegg_results_P) # 先放negative,再放 positive


