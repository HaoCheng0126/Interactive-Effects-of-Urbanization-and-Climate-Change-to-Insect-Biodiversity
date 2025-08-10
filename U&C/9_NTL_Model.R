#设置工作目录并准备好我们需要配对的数据框
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "15_NTL_Biodiversity/"
if(!dir.exists(outDir)) dir.create(outDir)
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))
# 添加 Year 列
predictsSites <- predictsSites %>%
  mutate(Year = as.integer(substr(Sample_end_latest, 1, 4)))

#读取TIF数据
# 安装必要的包（如果尚未安装）
# install.packages(c("terra", "sf", "raster", "tidyverse", "tmap", "RColorBrewer"))

# 加载包
library(terra)      # 新一代栅格数据处理包，比raster更快
library(sf)         # 处理矢量数据
library(raster)     # 栅格数据处理（与旧代码兼容）
library(tidyverse)  # 数据处理和可视化
library(tmap)       # 专题地图制作
library(RColorBrewer) # 色彩方案
library(lme4)
library(MuMIn)
library(predictsFunctions)
source("Functions.R")


# 设置工作目录路径
folder_path <- "/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/NTL Database/"

# 列出文件夹中的所有tif文件
dmsp_files <- list.files(path = folder_path, 
                        pattern = "\\.tif$", 
                        full.names = TRUE)

# 打印找到的文件列表
print(dmsp_files)

# 从文件名中提取年份 (根据您的文件命名格式调整)
get_year <- function(file_path) {
  # 提取文件名，假设格式如"DMSP_XXXX.tif"或"XXXX.tif"，其中XXXX是年份
  file_name <- basename(file_path)
  # 提取年份，调整正则表达式以匹配您的文件名格式
  year <- as.numeric(gsub(".*?([0-9]{4}).*\\.tif$", "\\1", file_name))
  return(year)
}

# 获取所有文件的年份
dmsp_years <- sapply(dmsp_files, get_year)
names(dmsp_files) <- as.numeric(dmsp_years)

# 检查年份是否正确提取
print(data.frame(File = basename(dmsp_files), Year = dmsp_years))

# 添加NTL列
predictsSites$NTL_DMSP <- NA

# 按年份提取NTL值
library(raster)
library(sp)

# 按年份提取NTL值
for(year in unique(predictsSites$Year)) {
  if(year %in% names(dmsp_files)) {
    # 读取该年份的栅格
    r <- raster(dmsp_files[as.character(year)])
    
    # 找到该年份的站点
    idx <- which(predictsSites$Year == year)
    
    # 提取坐标
    coords <- cbind(predictsSites$Longitude[idx], predictsSites$Latitude[idx])
    
    # 创建SpatialPoints对象
    sp_points <- SpatialPoints(coords, proj4string = CRS(proj4string(r)))
    
    # 创建1km的缓冲区
    buffer <- gBuffer(sp_points, width = 1000, byid = TRUE)
    
    # 提取缓冲区内的NTL值并计算平均值
    avg_ntl_values <- sapply(1:length(buffer), function(i) {
      masked_raster <- mask(r, buffer[i, ])
      mean_value <- cellStats(masked_raster, stat = 'mean', na.rm = TRUE)
      return(mean_value)
    })
    
    # 将平均值赋值给predictsSites
    predictsSites$NTL_DMSP[idx] <- avg_ntl_values
  }
}

# 画箱线图
# 保证UI2顺序与颜色一致（如已是该顺序可省略）
predictsSites$UI2 <- factor(
  predictsSites$UI2,
  levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High")
)

mycols <- c("#009E73", "#0072B2", "#E69F00", "#D55E00")

p <- ggplot(predictsSites, aes(x = UI2, y = NTL_DMSP, fill = UI2)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "grey") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() +
  xlab("Land Use Type") +
  ylab("Nighttime Light") +
  ggtitle("NTL Distribution by Land Use Type") +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 11, angle = 30, vjust = 0.9, hjust = 0.85, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),    # 纵坐标标题加粗
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 大标题加粗居中
    legend.position = "none"
  )

print(p)



##%######################################################%##
#                                                          #
####        Running the model selection process         ####
#                                                          #
##%######################################################%##
# 1. Abundance, NTL

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ] # 1984

NTLModelAbund <- glmer(LogAbund ~ UI2 + NTL_DMSP + 
                                UI2:NTL_DMSP + 
                                (1|SS) + (1|SSB),
                              data = model_data,
                       family=gaussian)

r_squared <- r.squaredGLMM(NTLModelAbund)
print(r_squared)

library(lmerTest)  # 提供F值和p值计算
library(car)       # 提供Anova()函数
# 1. 计算温度异常值与城市化强度(UI2)的交互作用
# 丰度模型
anova_abund_UI <- anova(NTLModelAbund)
print(anova_abund_UI)  # 查看"UI2:poly(StdTmeanAnomalyRS, 1)"行的F值和p值

# 或使用car包的Anova函数获取III型方差分析结果
anova_abund_UI_type3 <- Anova(NTLModelAbund, type = "III")
print(anova_abund_UI_type3)

# 步骤2: 拟合简化模型（去掉UI2和NTL_DMSP及其交互项）
NTLModelAbundNull <- glmer(LogAbund ~ 1 + 
                             (1|SS) + (1|SSB),
                           data = model_data,
                           family = gaussian)

# 步骤3: 使用anova()比较两个模型
anova_results <- anova(NTLModelAbundNull, NTLModelAbund)

# 输出结果
print(anova_results)
# 2. Richness, NTL

model_data <- predictsSites[!is.na(predictsSites$Species_richness), ] # 2170
# 去除UI2为"Primary vegetation"的数据
# model_data <- model_data[model_data$UI2 != "Primary vegetation", ] # 484
# model_data$UI2 <- droplevels(model_data$UI2) # 去除Primary vegetation的level

NTLModelRich <- glmer(Species_richness ~ UI2 + NTL_DMSP + 
                               UI2:NTL_DMSP + 
                               (1|SS) + (1|SSB),
                             data = model_data,
                     family = poisson)
r_squared <- r.squaredGLMM(NTLModelRich)
print(r_squared)
# 物种丰富度模型
anova_rich_UI <- anova(NTLModelRich)
print(anova_rich_UI)  # 查看"UI2:poly(StdTmeanAnomalyRS, 1)"行的F值和p值
# 步骤2: 拟合简化模型（去掉UI2和NTL_DMSP及其交互项）
NTLModelRichNull <- glmer(Species_richness ~ 1 + 
                            (1|SS) + (1|SSB),
                          data = model_data,
                          family = poisson)

# 步骤3: 使用anova()比较两个模型
anova_results <- anova(NTLModelRichNull, NTLModelRich)

# 输出结果
print(anova_results)
##%######################################################%##
#                                                          #
####                 Manuscript Figures                 ####
#                                                          #
##%######################################################%##

library(ggplot2)
library(dplyr)

# 设置分位数
exclQuantiles <- c(0.025, 0.975)
model_data <- predictsSites[!is.na(predictsSites$LogAbund), ] # 1984
# 创建预测数据框
nd <- expand.grid(
  NTL_DMSP = seq(from = min(model_data$NTL_DMSP),
                          to = max(model_data$NTL_DMSP),
                          length.out = 100),
  UI2 = factor(c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"),
               levels = levels(model_data$UI2)))

# 反向转换预测变量（如果你有这个函数的话）
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$NTL_DMSP,
  originalX = predictsSites$NTL_DMSP)

# 设置丰度为0（将被预测）
nd$LogAbund <- 0

# 设置随机效应为NA（群体水平预测）
nd$SS <- NA
nd$SSB <- NA  

# 参考行：主要植被且温度异常接近0
refRow <- which((nd$UI2 == "Primary vegetation") & 
                  (nd$NTL_DMSP == nd$NTL_DMSP[which.min(abs(nd$NTL_DMSP))]))
print(refRow)
# 计算各土地利用类型的分位数
QPV <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Primary vegetation"],
                probs = exclQuantiles)
QSV <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Urban_Low"],
                probs = exclQuantiles)
QAL <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Urban_Moderate"],
                probs = exclQuantiles)
QAH <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Urban_High"],
                probs = exclQuantiles)

# 检查并打印分位数范围
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "
")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "
")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "
")
cat("Urban_High range:", QAH[1], "to", QAH[2], "
")

# 如果你有 PredictGLMERRandIter 函数，使用它

a.preds.tmean <- PredictGLMERRandIter(model = NTLModelAbund, data = nd)


# 反向转换丰度值（从对数尺度转换回原始尺度）
a.preds.tmean <- exp(a.preds.tmean) - 0.01

# 转换为相对于参考的值（相对于主要植被在零异常温度下的值）
a.preds.tmean <- sweep(x = a.preds.tmean, MARGIN = 2, 
                       STATS = a.preds.tmean[refRow, ], FUN = '/')

# 移除超出分位数范围的值（确保预测在观测数据范围内）
a.preds.tmean[which((nd$UI2 == "Primary vegetation") & 
                      (nd$NTL_DMSP < QPV[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Primary vegetation") & 
                      (nd$NTL_DMSP > QPV[2])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Low") & 
                      (nd$NTL_DMSP < QSV[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Low") & 
                      (nd$NTL_DMSP > QSV[2])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Moderate") & 
                      (nd$NTL_DMSP < QAL[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Moderate") & 
                      (nd$NTL_DMSP > QAL[2])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_High") & 
                      (nd$NTL_DMSP < QAH[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_High") & 
                      (nd$NTL_DMSP > QAH[2])), ] <- NA

# 计算中位数和置信区间（转换为百分比变化）
nd$PredMedian <- ((apply(X = a.preds.tmean, MARGIN = 1, 
                         FUN = median, na.rm = TRUE)) * 100) - 100
nd$PredUpper <- ((apply(X = a.preds.tmean, MARGIN = 1, 
                        FUN = quantile, probs = 0.975, na.rm = TRUE)) * 100) - 100
nd$PredLower <- ((apply(X = a.preds.tmean, MARGIN = 1, 
                        FUN = quantile, probs = 0.025, na.rm = TRUE)) * 100) - 100

# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd$UI2, !is.na(nd$PredMedian)))


# 设置因子水平
nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Urban_Low", 
                                    "Urban_Moderate", "Urban_High"))

# 加载ggplot2
library(ggplot2)

# 绘图
p1 <- ggplot(data = nd, aes(x = NTL_DMSP, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,20,40,60), limits = c(0,63)) +
  scale_y_continuous(breaks = c( -75, -50, -25, 0, 25, 50,75), limits = c(-75, 75)) +
  ylab("Change in total abundance (%)") +
  xlab("Nighttime Light Intensity") +
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(), 
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2)) +
  ggtitle("a")

# 显示图形
print(p1)

##%######################################################%##
#                                                          #
####                 Manuscript Figures                 ####
#                                                          #
##%######################################################%##
## now the species richness plot
model_data <- predictsSites[!is.na(predictsSites$Species_richness), ] # 2170
exclQuantiles <- c(0.025,0.975)
nd2 <- expand.grid(
  NTL_DMSP = seq(from = min(model_data$NTL_DMSP),
                          to = max(model_data$NTL_DMSP),
                          length.out = 100),
  UI2 = factor(c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"),
               levels = levels(model_data$UI2)))


# back transform the predictors
nd2$NTL_DMSP <- BackTransformCentreredPredictor(
  transformedX = nd2$NTL_DMSP,
  originalX = predictsSites$NTL_DMSP)

# set richness and abundance to 0 - to be predicted
nd2$Species_richness <- 0
# 设置随机效应为NA（群体水平预测）
nd2$SS <- NA
nd2$SSB <- NA  

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2 == "Primary vegetation") & 
                  (nd2$NTL_DMSP == nd2$NTL_DMSP[which.min(abs(nd2$NTL_DMSP))]))
print(refRow)
# adjust plot 1: mean anomaly and abundance

# 计算各土地利用类型的分位数
QPV <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Primary vegetation"],
                probs = exclQuantiles)
QSV <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Urban_Low"],
                probs = exclQuantiles)
QAL <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Urban_Moderate"],
                probs = exclQuantiles)
QAH <- quantile(x = model_data$NTL_DMSP[model_data$UI2 == "Urban_High"],
                probs = exclQuantiles)

# 检查并打印分位数范围
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "
")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "
")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "
")
cat("Urban_High range:", QAH[1], "to", QAH[2], "
")

# 如果你有 PredictGLMERRandIter 函数，使用它

a.preds.tmean <- PredictGLMERRandIter(model = NTLModelRich, data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") 
                    & (nd2$NTL_DMSP < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation")
                    & (nd2$NTL_DMSP > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low")
                    & (nd2$NTL_DMSP < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low")
                    & (nd2$NTL_DMSP > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate")
                    & (nd2$NTL_DMSP < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate")
                    & (nd2$NTL_DMSP > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High")
                    & (nd2$NTL_DMSP < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High")
                    & (nd2$NTL_DMSP > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd2$UI2, !is.na(nd2$PredMedian)))

tapply(nd2$PredMedian, nd2$UI2, summary)
tapply(nd2$PredUpper, nd2$UI2, summary)
tapply(nd2$PredLower, nd2$UI2, summary)
# set factor levels
nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"))

nd2_sub <- subset(nd2, !is.na(NTL_DMSP) & !is.na(PredMedian) & NTL_DMSP >= 0 & NTL_DMSP <= 100)

# 检查每组数量
table(nd2_sub$UI2)
# plot 

p2 <- ggplot(data = nd2_sub, aes(x = NTL_DMSP, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 20,40,60), limits = c(0, 63)) +
  scale_y_continuous(breaks = c(-75,-50, -25, 0, 25, 50, 75,100), limits = c(-75, 75)) +
  ylab("Change in species richness (%)") +
  xlab("Nighttime Light Intensity") +
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(), 
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2)) +
  ggtitle("b")
print(p2)

library(cowplot)
# 去除主图图例
p1_nolegend <- p1 + theme(legend.position = "none")
p2_nolegend <- p2 + theme(legend.position = "none")

# 提取图例
leg <- get_legend(
  p1 + theme(legend.position = "right", 
             legend.justification = "center",
             legend.text = element_text(size = 7),
             legend.title = element_blank())
)

# 拼接主图
combined_plots <- plot_grid(p1_nolegend, p2_nolegend, nrow = 1, align = 'hv')

# 拼接主图和图例
final_plot <- plot_grid(combined_plots, leg, nrow = 1, rel_widths = c(1, 0.25))

# 显示
print(final_plot)


