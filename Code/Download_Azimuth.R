# 更改library路徑，因為pbmcref.SeuratData會因為權限問題而無法安裝
# .libPaths("C:/Users/benson/R/win-library/")


install.packages("devtools")

install.packages("httpuv")

install.packages("xtable")

devtools::install_github("satijalab/azimuth")


# 如果Azimuth已經安裝但是SHA1沒有做改變，想要強制重新安裝，可以用這個

# devtools::install_github("satijalab/azimuth", force = TRUE)


library(Azimuth)


# 如果還是不成功
# 確認Azimuth是否有成功安裝
installed.packages()

# 如果有存在，確認一下library的部分是否有大小寫的問題