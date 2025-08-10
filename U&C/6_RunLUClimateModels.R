##%######################################################%##
#                                                          #
####             Run the climate/LUI models             ####
#                                                          #
##%######################################################%##
# This script runs the mixed effects models looking at the effect of climate
# change and landuse/use intensity.
rm(list = ls())

# directories 
predictsDataDir <- "5_PREDICTSMatchPropNatHab/"
outDir <- "6_RunLUClimateModels/"

if(!dir.exists(outDir)) dir.create(outDir)

# sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)



# load libraries
library(devtools)
#install_github("timnewbold/StatisticalModels")
library(StatisticalModels)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions", force = T)
library(predictsFunctions)
source("Functions.R")
library(sjPlot)
library(cowplot)


###Create Models for all insects in predicts for standardised climate anomaly and Land interactions

# read in the predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesWithClimateAndNatHab.rds"))
predictsSites <- predictsSites@data

# remove the Argricture and Sencondary sites
predictsSites$UI2 <- factor(predictsSites$UI2)
predictsSites$UI2 <- relevel(predictsSites$UI2,ref="Primary vegetation")

# rescale the variable
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)
predictsSites$StdTmaxAnomalyRS <- StdCenterPredictor(predictsSites$StdTmaxAnomaly)

# organise NH data - combine primary and secondary vegetation at each scale
predictsSites$NH_1000 <- predictsSites$PV_1000 + predictsSites$SV_1000
predictsSites$NH_3000 <- predictsSites$PV_3000 + predictsSites$SV_3000
predictsSites$NH_5000 <- predictsSites$PV_5000 + predictsSites$SV_5000
predictsSites$NH_10000 <- predictsSites$PV_10000 + predictsSites$SV_10000


# rescale the variables
predictsSites$NH_1000.rs <- StdCenterPredictor(predictsSites$NH_1000)
predictsSites$NH_3000.rs <- StdCenterPredictor(predictsSites$NH_3000)
predictsSites$NH_5000.rs <- StdCenterPredictor(predictsSites$NH_5000)
predictsSites$NH_10000.rs <- StdCenterPredictor(predictsSites$NH_10000)

# rescaling abundance and log values
predictsSites <- RescaleAbundance(predictsSites)

# charlie added this line as later bits were throwing errors
predictsSites <- droplevels(predictsSites)

# some of the climate values are NA since they do not meet the thresholds
predictsSites <- predictsSites[!is.na(predictsSites$avg_temp), ]

# 图1：标准化变量的相关性分析
plot(predictsSites$StdTmeanAnomaly, predictsSites$StdTmaxAnomaly,
     xlab = "StdTmeanAnomaly",
     ylab = "StdTmaxAnomaly",
     main = "Correlation: StdTmeanAnomaly vs StdTmaxAnomaly")

# 计算标准化变量的相关系数,use = "complete.obs"去除缺失值的影响
std_cor <- cor(predictsSites$StdTmeanAnomaly, predictsSites$StdTmaxAnomaly, 
               use = "complete.obs")
cor(predictsSites$StdTmeanAnomaly, predictsSites$StdTmaxAnomaly, 
    use = "complete.obs")
# 在图上显示相关系数
legend("topleft", legend = paste("Correlation:", round(std_cor, 3)), 
       bty = "n")

# 图2：非标准化变量的相关性分析
par(new = FALSE)  # 确保开始新图
plot(predictsSites$TmeanAnomaly, predictsSites$TmaxAnomaly,
     xlab = "TmeanAnomaly",
     ylab = "TmaxAnomaly",
     main = "Correlation: TmeanAnomaly vs TmaxAnomaly")

# 计算非标准化变量的相关系数
non_std_cor <- cor(predictsSites$TmeanAnomaly, predictsSites$TmaxAnomaly, 
                   use = "complete.obs")

# 在图上显示相关系数
legend("topleft", legend = paste("Correlation:", round(non_std_cor, 3)), 
       bty = "n")

# 打印两个相关系数
cat("Standardized correlation:", std_cor, "\n")
cat("Non-standardized correlation:", non_std_cor, "\n")


##%######################################################%##
#                                                          #
####     箱线图看各个landuse的温度异常的分布情况        ####
#                                                          #
##%######################################################%##
# 设置基本图形参数
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# 设置颜色
colors <- c("#71A8CF", "#7FC8A9")

# 重新排序因子水平
predictsSites$UI2 <- factor(predictsSites$UI2, 
                            levels = c("Primary vegetation", "Urban_Low", 
                                       "Urban_Moderate", "Urban_High"))

# 图1：平均温度异常值分布
boxplot(StdTmeanAnomaly ~ UI2, data = predictsSites, 
        main = "Mean Temperature Anomaly",
        xlab = "",
        ylab = "Standardized Temperature Anomaly",
        col = colors[1],
        border = "gray30",
        outpch = 20,
        outcol = "gray40")

# 添加均值点
means <- tapply(predictsSites$StdTmeanAnomaly, predictsSites$UI2, mean, na.rm = TRUE)
points(1:length(means), means, col = "darkred", pch = 18, cex = 1.3)

# 添加整体平均线 - 使用更优雅的颜色
overall_mean <- mean(predictsSites$StdTmeanAnomaly, na.rm = TRUE)
abline(h = overall_mean, col = "red", lty = 2, lwd = 1.5)


# 图2：最高温度异常值分布
boxplot(StdTmaxAnomaly ~ UI2, data = predictsSites, 
        main = "Maximum Temperature Anomaly",
        xlab = "",
        ylab = "Standardized Max Temp Anomaly",
        col = colors[2],
        border = "gray30",
        outpch = 20,
        outcol = "gray40")

# 添加均值点
means_max <- tapply(predictsSites$StdTmaxAnomaly, predictsSites$UI2, mean, na.rm = TRUE)
points(1:length(means_max), means_max, col = "darkred", pch = 18, cex = 1.3)

# 添加整体平均线 - 使用更优雅的颜色
overall_mean_max <- mean(predictsSites$StdTmaxAnomaly, na.rm = TRUE)
abline(h = overall_mean_max, col = "red", lty = 2, lwd = 1.5)

# 重置绘图参数
par(mfrow = c(1, 1))

# take a look at possible correlations between variables
plot(predictsSites$avg_temp, predictsSites$TmeanAnomaly,
     xlab = "Average temperature", 
     ylab = "Anomaly (difference between present and baseline)")
cor(predictsSites$avg_temp, predictsSites$TmeanAnomaly) # -0.2309


plot(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly,
     xlab = "Average temperature", 
     ylab = "Standardised climate anomaly")
cor(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly) # 0.2453451

cor(predictsSites$TmeanAnomaly, predictsSites$StdTmeanAnomaly) # 0.2022808


# save the dataset
saveRDS(object = predictsSites,file = paste0(outDir,"PREDICTSSiteData.rds"))
#predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSiteData.rds"))

##%######################################################%##
#                                                          #
####        Running the model selection process         ####
#                                                          #
##%######################################################%##

# 1. Abundance, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ] # 1984
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ] # 1984
# 确保是日期格式
model_data$Sample_start_earliest <- as.Date(model_data$Sample_start_earliest)

# 提取年份
model_data$Year <- as.numeric(format(model_data$Sample_start_earliest, "%Y"))

# 检查结果
head(model_data[, c("Sample_start_earliest", "Year")])
MeanAnomalyModelAbund <- lmer(LogAbund ~ UI2 + poly(StdTmeanAnomalyRS,1) + 
                               UI2:poly(StdTmeanAnomalyRS,1) + 
                                (1|SS) + (1|SSB),
                             data = model_data)
r_squared <- r.squaredGLMM(MeanAnomalyModelAbund)
print(r_squared)
library(lmerTest)  # 提供F值和p值计算
library(car)       # 提供Anova()函数
# 1. 计算温度异常值与城市化强度(UI2)的交互作用
# 丰度模型
anova_abund_UI <- anova(MeanAnomalyModelAbund)
print(anova_abund_UI)  # 查看"UI2:poly(StdTmeanAnomalyRS, 1)"行的F值和p值

# 或使用car包的Anova函数获取III型方差分析结果
anova_abund_UI_type3 <- Anova(MeanAnomalyModelAbund, type = "III")
print(anova_abund_UI_type3)
# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "/MeanAnomalyModelAbund.rdata"))
# 2. Richness, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 2170
# 确保是日期格式
model_data$Sample_start_earliest <- as.Date(model_data$Sample_start_earliest)

# 提取年份
model_data$Year <- as.numeric(format(model_data$Sample_start_earliest, "%Y"))

# 检查结果
head(model_data[, c("Sample_start_earliest", "Year")])

MeanAnomalyModelRich <- lmer(Species_richness ~ UI2 + poly(StdTmeanAnomalyRS,1) + 
                               UI2:poly(StdTmeanAnomalyRS,1) + 
                               (1|SS) + (1|SSB) + (1|Year) + (1|n_months),
                             data = model_data)
r_squared <- r.squaredGLMM(MeanAnomalyModelRich)
print(r_squared)
# 物种丰富度模型
anova_rich_UI <- anova(MeanAnomalyModelRich)
print(anova_rich_UI)  # 查看"UI2:poly(StdTmeanAnomalyRS, 1)"行的F值和p值

save(MeanAnomalyModelRich, file = paste0(outDir, "/MeanAnomalyModelRich.rdata"))
# 3. Abundance, max anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ] # 1984
model_data <- model_data[!is.na(model_data$StdTmaxAnomalyRS), ] # 1984

# 确保是日期格式
model_data$Sample_start_earliest <- as.Date(model_data$Sample_start_earliest)

# 提取年份
model_data$Year <- as.numeric(format(model_data$Sample_start_earliest, "%Y"))

# 检查结果
head(model_data[, c("Sample_start_earliest", "Year")])


MaxAnomalyModelAbund <- lmer(LogAbund ~ UI2 + poly(StdTmaxAnomalyRS,1) + 
                               UI2:poly(StdTmaxAnomalyRS,1) + 
                               (1|SS) + (1|SSB) + (1|Year) + (1|n_months),
                               data = model_data)
r_squared <- r.squaredGLMM(MaxAnomalyModelAbund)
print(r_squared)  
 
# 4. Richness, max anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmaxAnomalyRS), ] # 2170

# 确保是日期格式
model_data$Sample_start_earliest <- as.Date(model_data$Sample_start_earliest)

# 提取年份
model_data$Year <- as.numeric(format(model_data$Sample_start_earliest, "%Y"))

# 检查结果
head(model_data[, c("Sample_start_earliest", "Year")])

MaxAnomalyModelRich <- lmer(Species_richness ~ UI2 + poly(StdTmaxAnomalyRS,1) + 
                              UI2:poly(StdTmaxAnomalyRS,1) + 
                              (1|SS) + (1|SSB) + (1|Year) + (1|n_months),
                            data = model_data)

r_squared <- r.squaredGLMM(MaxAnomalyModelRich)
print(r_squared)  

##%######################################################%##
#                                                          #
####                 Manuscript Figures                 ####
#                                                          #
##%######################################################%##
library(ggplot2)
library(dplyr)

# 设置分位数
exclQuantiles <- c(0.025, 0.975)

# 创建预测数据框
nd <- expand.grid(
  StdTmeanAnomalyRS = seq(from = min(model_data$StdTmeanAnomalyRS),
                          to = max(model_data$StdTmeanAnomalyRS),
                          length.out = 100),
  UI2 = factor(c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"),
               levels = levels(model_data$UI2)))

# 反向转换预测变量（如果你有这个函数的话）
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
   transformedX = nd$StdTmeanAnomalyRS,
   originalX = predictsSites$StdTmeanAnomaly)

# 设置丰度为0（将被预测）
nd$LogAbund <- 0

# 设置随机效应为NA（群体水平预测）
nd$SS <- NA
nd$SSB <- NA  
nd$Year <- NA
nd$n_months <- NA

# 参考行：主要植被且温度异常接近0
refRow <- which((nd$UI2 == "Primary vegetation") & 
                  (nd$StdTmeanAnomaly == nd$StdTmeanAnomaly[which.min(abs(nd$StdTmeanAnomaly))]))
print(refRow)
# 计算各土地利用类型的分位数
QPV <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Primary vegetation"],
                probs = exclQuantiles)
QSV <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Urban_Low"],
                probs = exclQuantiles)
QAL <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Urban_Moderate"],
                probs = exclQuantiles)
QAH <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Urban_High"],
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

a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund, data = nd)


# 反向转换丰度值（从对数尺度转换回原始尺度）
a.preds.tmean <- exp(a.preds.tmean) - 0.01

# 转换为相对于参考的值（相对于主要植被在零异常温度下的值）
a.preds.tmean <- sweep(x = a.preds.tmean, MARGIN = 2, 
                       STATS = a.preds.tmean[refRow, ], FUN = '/')

# 移除超出分位数范围的值（确保预测在观测数据范围内）
a.preds.tmean[which((nd$UI2 == "Primary vegetation") & 
                      (nd$StdTmeanAnomalyRS < QPV[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Primary vegetation") & 
                      (nd$StdTmeanAnomalyRS > QPV[2])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Low") & 
                      (nd$StdTmeanAnomalyRS < QSV[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Low") & 
                      (nd$StdTmeanAnomalyRS > QSV[2])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Moderate") & 
                      (nd$StdTmeanAnomalyRS < QAL[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_Moderate") & 
                      (nd$StdTmeanAnomalyRS > QAL[2])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_High") & 
                      (nd$StdTmeanAnomalyRS < QAH[1])), ] <- NA
a.preds.tmean[which((nd$UI2 == "Urban_High") & 
                      (nd$StdTmeanAnomalyRS > QAH[2])), ] <- NA

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
p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1,1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Mean Temperature Anomaly") +
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


## now the species richness plot
model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]
exclQuantiles <- c(0.025,0.975)
nd2 <- expand.grid(
  StdTmeanAnomalyRS = seq(from = min(model_data$StdTmeanAnomalyRS),
                          to = max(model_data$StdTmeanAnomalyRS),
                          length.out = 100),
  UI2 = factor(c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"),
               levels = levels(model_data$UI2)))


# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0
# 设置随机效应为NA（群体水平预测）
nd2$SS <- NA
nd2$SSB <- NA  
nd2$Year <- NA
nd2$n_months <- NA

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2 == "Primary vegetation") & 
                  (nd2$StdTmeanAnomaly == nd2$StdTmeanAnomaly[which.min(abs(nd2$StdTmeanAnomaly))]))
print(refRow)
# adjust plot 1: mean anomaly and abundance

# 计算各土地利用类型的分位数
QPV <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Primary vegetation"],
                probs = exclQuantiles)
QSV <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Urban_Low"],
                probs = exclQuantiles)
QAL <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Urban_Moderate"],
                probs = exclQuantiles)
QAH <- quantile(x = model_data$StdTmeanAnomalyRS[model_data$UI2 == "Urban_High"],
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

a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich, data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"))

# plot 
p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Mean Temperature Anomaly") +
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
# combine plots
cowplot::plot_grid(p1, p2)

### Extended Data 2 - maximum anomaly ###### Extended Data 2 - maximum anomaly ###
### Extended Data 2 - maximum anomaly ###### Extended Data 2 - maximum anomaly ###
### Extended Data 2 - maximum anomaly ###### Extended Data 2 - maximum anomaly ###
### Extended Data 2 - maximum anomaly ###### Extended Data 2 - maximum anomaly ###
# 设置分位数
exclQuantiles <- c(0.025, 0.975)
model_data <- predictsSites[!is.na(predictsSites$LogAbund), ] # 1984
model_data <- model_data[!is.na(model_data$StdTmaxAnomalyRS), ] # 1984

nd3 <- expand.grid(
  StdTmaxAnomalyRS = seq(from = min(model_data$StdTmaxAnomalyRS),
                         to = max(model_data$StdTmaxAnomalyRS),
                         length.out = 100),
  UI2 = factor(c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"),
               levels = levels(model_data$UI2)))

# 反向转换预测变量（如果你有这个函数的话）
nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# 设置丰度为0（将被预测）
nd3$LogAbund <- 0

# 设置随机效应为NA（群体水平预测）
nd3$SS <- NA
nd3$SSB <- NA  
nd3$Year <- NA
nd3$n_months <- NA

# 参考行：主要植被且温度异常接近0
refRow <- which((nd3$UI2 == "Primary vegetation") & 
                  (nd3$StdTmaxAnomaly == nd3$StdTmaxAnomaly[which.min(abs(nd3$StdTmaxAnomaly))]))
print(refRow)

# 计算各土地利用类型的分位数
QPV <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Primary vegetation"],
                probs = exclQuantiles)
QSV <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Urban_Low"],
                probs = exclQuantiles)
QAL <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Urban_Moderate"],
                probs = exclQuantiles)
QAH <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Urban_High"],
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
a.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelAbund, data = nd3)

# 反向转换丰度值（从对数尺度转换回原始尺度）
a.preds.tmean <- exp(a.preds.tmean) - 0.01

# 转换为相对于参考的值（相对于主要植被在零异常温度下的值）
a.preds.tmean <- sweep(x = a.preds.tmean, MARGIN = 2, 
                       STATS = a.preds.tmean[refRow, ], FUN = '/')

# 移除超出分位数范围的值（确保预测在观测数据范围内）
a.preds.tmean[which((nd3$UI2 == "Primary vegetation") & 
                      (nd3$StdTmaxAnomalyRS < QPV[1])), ] <- NA
a.preds.tmean[which((nd3$UI2 == "Primary vegetation") & 
                      (nd3$StdTmaxAnomalyRS > QPV[2])), ] <- NA
a.preds.tmean[which((nd3$UI2 == "Urban_Low") & 
                      (nd3$StdTmaxAnomalyRS < QSV[1])), ] <- NA
a.preds.tmean[which((nd3$UI2 == "Urban_Low") & 
                      (nd3$StdTmaxAnomalyRS > QSV[2])), ] <- NA
a.preds.tmean[which((nd3$UI2 == "Urban_Moderate") & 
                      (nd3$StdTmaxAnomalyRS < QAL[1])), ] <- NA
a.preds.tmean[which((nd3$UI2 == "Urban_Moderate") & 
                      (nd3$StdTmaxAnomalyRS > QAL[2])), ] <- NA
a.preds.tmean[which((nd3$UI2 == "Urban_High") & 
                      (nd3$StdTmaxAnomalyRS < QAH[1])), ] <- NA
a.preds.tmean[which((nd3$UI2 == "Urban_High") & 
                      (nd3$StdTmaxAnomalyRS > QAH[2])), ] <- NA

# 计算中位数和置信区间（转换为百分比变化）
nd3$PredMedian <- ((apply(X = a.preds.tmean, MARGIN = 1, 
                          FUN = median, na.rm = TRUE)) * 100) - 100
nd3$PredUpper <- ((apply(X = a.preds.tmean, MARGIN = 1, 
                         FUN = quantile, probs = 0.975, na.rm = TRUE)) * 100) - 100
nd3$PredLower <- ((apply(X = a.preds.tmean, MARGIN = 1, 
                         FUN = quantile, probs = 0.025, na.rm = TRUE)) * 100) - 100

# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd3$UI2, !is.na(nd3$PredMedian)))

# 设置因子水平
nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Urban_Low", 
                                      "Urban_Moderate", "Urban_High"))

# 加载ggplot2
library(ggplot2)

# 绘图
p3 <- ggplot(data = nd3, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1,1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Max Temperature Anomaly") +
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
print(p3)

###第四幅图
# 设置分位数
exclQuantiles <- c(0.025, 0.975)
model_data <- predictsSites[!is.na(predictsSites$StdTmaxAnomalyRS), ] # 2170

nd4 <- expand.grid(
  StdTmaxAnomalyRS = seq(from = min(model_data$StdTmaxAnomalyRS),
                         to = max(model_data$StdTmaxAnomalyRS),
                         length.out = 100),
  UI2 = factor(c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"),
               levels = levels(model_data$UI2)))

# 反向转换预测变量（如果你有这个函数的话）
nd4$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# 设置物种丰富度为0（将被预测）
nd4$Species_richness <- 0

# 设置随机效应为NA（群体水平预测）
nd4$SS <- NA
nd4$SSB <- NA  
nd4$Year <- NA
nd4$n_months <- NA

# 参考行：主要植被且温度异常接近0
refRow <- which((nd4$UI2 == "Primary vegetation") & 
                  (nd4$StdTmaxAnomaly == nd4$StdTmaxAnomaly[which.min(abs(nd4$StdTmaxAnomaly))]))
print(refRow)

# 计算各土地利用类型的分位数
QPV <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Primary vegetation"],
                probs = exclQuantiles)
QSV <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Urban_Low"],
                probs = exclQuantiles)
QAL <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Urban_Moderate"],
                probs = exclQuantiles)
QAH <- quantile(x = model_data$StdTmaxAnomalyRS[model_data$UI2 == "Urban_High"],
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
a.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelRich, data = nd4)

# 反向转换物种丰富度值（从对数尺度转换回原始尺度）
a.preds.tmax <- exp(a.preds.tmax) - 0.01

# 转换为相对于参考的值（相对于主要植被在零异常温度下的值）
a.preds.tmax <- sweep(x = a.preds.tmax, MARGIN = 2, 
                      STATS = a.preds.tmax[refRow, ], FUN = '/')

# 移除超出分位数范围的值（确保预测在观测数据范围内）
a.preds.tmax[which((nd4$UI2 == "Primary vegetation") & 
                     (nd4$StdTmaxAnomalyRS < QPV[1])), ] <- NA
a.preds.tmax[which((nd4$UI2 == "Primary vegetation") & 
                     (nd4$StdTmaxAnomalyRS > QPV[2])), ] <- NA
a.preds.tmax[which((nd4$UI2 == "Urban_Low") & 
                     (nd4$StdTmaxAnomalyRS < QSV[1])), ] <- NA
a.preds.tmax[which((nd4$UI2 == "Urban_Low") & 
                     (nd4$StdTmaxAnomalyRS > QSV[2])), ] <- NA
a.preds.tmax[which((nd4$UI2 == "Urban_Moderate") & 
                     (nd4$StdTmaxAnomalyRS < QAL[1])), ] <- NA
a.preds.tmax[which((nd4$UI2 == "Urban_Moderate") & 
                     (nd4$StdTmaxAnomalyRS > QAL[2])), ] <- NA
a.preds.tmax[which((nd4$UI2 == "Urban_High") & 
                     (nd4$StdTmaxAnomalyRS < QAH[1])), ] <- NA
a.preds.tmax[which((nd4$UI2 == "Urban_High") & 
                     (nd4$StdTmaxAnomalyRS > QAH[2])), ] <- NA

# 计算中位数和置信区间（转换为百分比变化）
nd4$PredMedian <- ((apply(X = a.preds.tmax, MARGIN = 1, 
                          FUN = median, na.rm = TRUE)) * 100) - 100
nd4$PredUpper <- ((apply(X = a.preds.tmax, MARGIN = 1, 
                         FUN = quantile, probs = 0.975, na.rm = TRUE)) * 100) - 100
nd4$PredLower <- ((apply(X = a.preds.tmax, MARGIN = 1, 
                         FUN = quantile, probs = 0.025, na.rm = TRUE)) * 100) - 100

# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd4$UI2, !is.na(nd4$PredMedian)))

# 设置因子水平
nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Urban_Low", 
                                      "Urban_Moderate", "Urban_High"))

# 加载ggplot2
library(ggplot2)

# 绘图
p4 <- ggplot(data = nd4, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1,1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Max Temperature Anomaly") +
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
  ggtitle("d")

# 显示图形
print(p4)
