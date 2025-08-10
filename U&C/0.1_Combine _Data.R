##%######################################################%##
#                                                          #
####         Organise PREDICTS data for insects         ####
#                                                          #
##%######################################################%##

# This script takes the complete PREDICTS database, selects those entries for
# insects, and organises the data for analysis.

# load required libraries
library(predictsFunctions)
library(ggplot2)
library(purrr)

# show the R environment
sessionInfo()

# Load required library (for data manipulation)
library(dplyr)

# Read the two CSV files into data frames
data1 <- readRDS("/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/0_data/database_01.rds") 
data2 <- readRDS("/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/0_data/database_02.rds")

# 3. 检查两个数据集的列名差异
common_cols <- intersect(names(data1), names(data2))  # 共有的列名
unique_to_data1 <- setdiff(names(data1), names(data2))  # data1独有的列名
unique_to_data2 <- setdiff(names(data2), names(data1))  # data2独有的列名

# 打印列名差异
cat("共有的列名：", paste(common_cols, collapse = ", "), "\n")
cat("resource_01.csv独有的列名：", paste(unique_to_data1, collapse = ", "), "\n")
cat("resource_02.csv独有的列名：", paste(unique_to_data2, collapse = ", "), "\n")

# 获取两个数据框共有的列
common_cols <- intersect(names(data1), names(data2))

# 检测类型不一致的列
conflict_cols <- common_cols[map_lgl(common_cols, function(col) {
  class(data1[[col]])[1] != class(data2[[col]])[1]
})]

# 打印冲突列信息
cat("发现以下列的类型冲突:
")
for(col in conflict_cols) {
  cat(col, ": data1是", class(data1[[col]])[1], 
      ", data2是", class(data2[[col]])[1], "
")
}

# 将冲突列都转换为字符型
data1 <- data1 %>%
  mutate(across(all_of(conflict_cols), as.character))

data2 <- data2 %>%
  mutate(across(all_of(conflict_cols), as.character))

# 合并数据
combined_data <- bind_rows(data1, data2)

cat("已成功合并数据
")

# 保存数据框为 RDS 文件
saveRDS(combined_data, file = "database.rds")

data3 <- readRDS("/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/0_data/sites_01.rds") 
data4 <- readRDS("/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/0_data/sites_02.rds")
# 合并数据
combined_data1 <- bind_rows(data3, data4)

cat("已成功合并数据
")
saveRDS(combined_data1, file = "sites.rds")
