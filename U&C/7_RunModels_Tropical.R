##%######################################################%##
#                                                          #
####    Models assessing differences between realms     ####
#                                                          #
##%######################################################%##


# In this script, separate models are run for subsets of the PREDICTS database
# based on which realm the sites are found in. 

# Inclusion of Realm in the model results in high multicolinearity between variables
# so subsets of the data and separate models are run instead.

# directories
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "9_RunModels_Tropical/"
if(!dir.exists(outDir)) dir.create(outDir)

# sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)


# load libraries
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")
library(ggplot2)
library(cowplot)
library(sjPlot)
library(performance)
#library(pedometrics)
#library(car)
library(lme4)
library(performance)
library(MuMIn)

# read in PREDICTS data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

#predictsSites$UI2[predictsSites$UI2 == "Urban_Moderate"] <- "Urban_Low"
# 删除未使用的因子水平
predictsSites$UI2 <- droplevels(predictsSites$UI2)

# 验证结果
table(predictsSites$UI2)  # 确认没有"Urban_Moderate"类别

predictsSites$UI2 <- droplevels(predictsSites$UI2)
#Primary vegetation          Urban_Low         Urban_High 
# 1686                126                358 

# categorise sites into temperate/tropical
# get the tropical values
predictsSites$Tropical <- NA

predictsSites[predictsSites$Latitude > -23.44 & predictsSites$Latitude < 23.44, 'Tropical'] <- "Tropical"

# label the remaining as temperate
predictsSites[is.na(predictsSites$Tropical), 'Tropical'] <- "Temperate"

# set as a factor
predictsSites$Tropical <- as.factor(predictsSites$Tropical)

# levels: Temperate Tropical
predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] #2170

table(predictsSites$Tropical)

# Temperate  Tropical 
# 1477     693 

# abundance data subset
table(predictsSites[!is.na(predictsSites$LogAbund), 'Tropical'])

# Temperate  Tropical 
# 1356      628

# split by land use classes
table(predictsSites$UI2, predictsSites$Tropical)

#                       Temperate Tropical
# Primary vegetation      1016      670
# Urban_Low                115       11
# Urban_High               346       12

# look into range of SCA values across landuses in tropical/temperate areas.

plot_data <- predictsSites[, c("Tropical", "UI2", "StdTmeanAnomaly")]
plot_data <- plot_data[!is.na(plot_data$UI2), ]
table(plot_data$UI2)
plot_data$UI2 <- factor(plot_data$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"))

ggplot(data = plot_data) +
  geom_freqpoly(aes(x = StdTmeanAnomaly, col = UI2)) +
  facet_wrap(~ Tropical) + 
  theme_bw()

plot_data <- predictsSites[, c("Tropical", "UI2", "StdTmaxAnomaly")]
plot_data <- plot_data[!is.na(plot_data$UI2), ]
table(plot_data$UI2)
plot_data$UI2 <- factor(plot_data$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High"))

ggplot(data = plot_data) +
  geom_freqpoly(aes(x = StdTmaxAnomaly, col = UI2)) +
  facet_wrap(~ Tropical) + 
  theme_bw()

### look at correlations between climate metrics for realm subsets ###

trop <- predictsSites[predictsSites$Tropical == "Tropical", ]

temp <- predictsSites[predictsSites$Tropical == "Temperate", ]

# trop
cor(trop$avg_temp, trop$TmeanAnomaly) # -0.1253853
cor(trop$avg_temp, trop$StdTmeanAnomaly) # 0.1836093
cor(trop$TmeanAnomaly, trop$StdTmeanAnomaly) # 0.4521215

#temp
cor(temp$avg_temp, temp$TmeanAnomaly) # 0.218891
cor(temp$avg_temp, temp$StdTmeanAnomaly) # -0.6341445
cor(temp$TmeanAnomaly, temp$StdTmeanAnomaly) # 0.3967754



##%######################################################%##
#                                                          #
####               Models for each realm                ####
#                                                          #
##%######################################################%##


# since the model including Realm and its interactions with other variables
# suffers from multicolinearity, we are running a separate model for subsets of data 
# from each realm to see if the patterns observed hold. 


# use dataset with realm info created above
head(predictsSites)

# sort dataset, remove NAs etc
predictsSites <- predictsSites[!is.na(predictsSites$UI2), ]
predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 2170


levels(predictsSites$UI2)
table(predictsSites$UI2)

predictsSites$UI2 <- factor(predictsSites$UI2 , levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate",  "Urban_High" ))

levels(predictsSites$UI2)
table(predictsSites$UI2)


predictsSites <- droplevels(predictsSites)


# SR subsets
predictsSites_sr_trop <- predictsSites[predictsSites$Tropical == "Tropical", ] # 693
predictsSites_sr_temp <- predictsSites[predictsSites$Tropical == "Temperate", ] # 1477

length(unique(predictsSites_sr_trop$SS)) # 93
length(unique(predictsSites_sr_temp$SS)) # 91

table(predictsSites_sr_trop$UI2)
table(predictsSites_sr_temp$UI2)

# subset for abundance data
predictsSites_ab <- predictsSites[!is.na(predictsSites$LogAbund), ] # 1984

## dataset summaries ##

predictsSites_ab_trop <- predictsSites_ab[predictsSites_ab$Tropical == "Tropical", ] # 628
predictsSites_ab_temp <- predictsSites_ab[predictsSites_ab$Tropical == "Temperate", ] # 1356
  

length(unique(predictsSites_ab_trop$SS)) # 82
length(unique(predictsSites_ab_temp$SS)) # 83

table(predictsSites_ab_trop$UI2)

# Primary vegetation          Urban_Low         Urban_High 
# 605                 11                 12 

table(predictsSites_ab_temp$UI2)

# Primary vegetation          Urban_Low         Urban_High 
# 911                110                335 


##%######################################################%##
#                                                          #
####                    run models                      ####
#                                                          #
##%######################################################%##


# 1. Temperate, Abundance, mean anomaly

MeanAbundTemp <- GLMERSelect(modelData = predictsSites_ab_temp,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MeanAbundTemp$model)
# 提取模型残差
model_residuals <- residuals(MeanAbundTemp$model)

# 基本残差图
par(mfrow = c(2, 2))
plot(MeanAbundTemp$model)

# 计算残差平均值
mean_residual <- mean(model_residuals)
print(paste("残差平均值:", mean_residual))


# 检查正负残差比例
positive_residuals <- sum(model_residuals > 0)
negative_residuals <- sum(model_residuals < 0)
print(paste("正残差数量:", positive_residuals))
print(paste("负残差数量:", negative_residuals))

observed_values <- MeanAbundTemp$model@frame$LogAbund  # 提取观测的对数丰度值
fitted_values <- fitted(MeanAbundTemp$model)           # 提取模型拟合值

# 2. 创建观测值与拟合值的散点图
plot(fitted_values, observed_values, 
     xlab = "模型预测值", ylab = "实际观测值", 
     main = "观测值与预测值对比图",
     pch = 16, col = "darkblue", cex = 0.8)

# 3. 添加1:1参考线（理想情况下点应分布在此线附近）
abline(0, 1, col = "red", lty = 2, lwd = 2)

# 4. 添加趋势线，查看实际关系
fit_line <- lm(observed_values ~ fitted_values)
abline(fit_line, col = "green3", lwd = 2)

# R2GLMER(MeanAbundTemp$model)
# LogAbund ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS, 1) + (1 | SS) + (1 | SSB)

# 2. Tropical, Abundance, mean anomaly
MeanAbundTrop <- GLMERSelect(modelData = predictsSites_ab_trop,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MeanAbundTrop$model)
# LogAbund ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS, 1) + (1 | SS) + (1 | SSB)


# 3. Temperate, Richness, mean anomaly
MeanRichTemp <- GLMERSelect(modelData = predictsSites_sr_temp,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS"))

summary(MeanRichTemp$model)
#Species_richness ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS, 1) + (1 | SS) + (1 | SSB) + (1 | SSBS)
# 提取模型残差
model_residuals <- residuals(MeanRichTemp$model)

# 基本残差图
par(mfrow = c(2, 2))
plot(MeanRichTemp$model)

# 计算残差平均值
mean_residual <- mean(model_residuals)
print(paste("残差平均值:", mean_residual))

# 检查正负残差比例
positive_residuals <- sum(model_residuals > 0)
negative_residuals <- sum(model_residuals < 0)
print(paste("正残差数量:", positive_residuals))
print(paste("负残差数量:", negative_residuals))

# 4. Tropical, Richness, mean anomaly
MeanRichTrop <- GLMERSelect(modelData = predictsSites_sr_trop,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS"))

summary(MeanRichTrop$model)
# Species_richness ~ UI2 + UI2:poly(StdTmeanAnomalyRS, 1) + poly(StdTmeanAnomalyRS,  1) + (1 | SS) + (1 | SSB) + (1 | SSBS)

save(MeanAbundTemp, file = paste0(outDir, "/MeanAbundTemp_output.rdata"))
save(MeanAbundTrop, file = paste0(outDir, "/MeanAbundTrop_output.rdata"))
save(MeanRichTemp, file = paste0(outDir, "/MeanRichTemp_output.rdata"))
save(MeanRichTrop, file = paste0(outDir, "/MeanRichTrop_output.rdata"))

#load(file = paste0(outDir, "/MeanAbundTemp_output.rdata"))
#load(file = paste0(outDir, "/MeanAbundTrop_output.rdata"))
#load(file = paste0(outDir, "/MeanRichTemp_output.rdata"))
#load(file = paste0(outDir, "/MeanRichTrop_output.rdata"))


#### save model output tables using sjPlot package ####

tab_model(MeanAbundTemp$model, MeanAbundTrop$model, transform = NULL, file = paste0(outDir, "/AbunMeanAnomTempTrop_output_table.html"))
summary(MeanAbundTemp$model)
R2GLMER(MeanAbundTemp$model)
summary(MeanAbundTrop$model)
R2GLMER(MeanAbundTrop$model)

tab_model(MeanRichTemp$model, MeanRichTrop$model, transform = NULL, file = paste0(outDir, "/RichMeanAnomTempTrop_output_table.html"))
summary(MeanRichTemp$model)
R2GLMER(MeanRichTemp$model) # use these values
# $conditional
# [1] 0.6652306
# $marginal
# [1] 0.00808889
summary(MeanRichTrop$model)
R2GLMER(MeanRichTrop$model) # use these values
# $conditional
# [1] 0.5596599
# $marginal
# [1] 0.03645811

# save model stats

# save the stats info
abtemp_stats <- as.data.frame(MeanAbundTemp$stats)
abtrop_stats <- as.data.frame(MeanAbundTrop$stats)
srtemp_stats <- as.data.frame(MeanRichTemp$stats)
srtrop_stats <- as.data.frame(MeanRichTrop$stats)

abtemp_stats$significant <- NA
abtrop_stats$significant <- NA
srtemp_stats$significant <- NA
srtrop_stats$significant <- NA

abtemp_stats$model <- "abtemp"
abtrop_stats$model <- "abtrop"
srtemp_stats$model <- "srtemp"
srtrop_stats$model <- "srtrop"

all_stats <- rbind(abtemp_stats, abtrop_stats, srtemp_stats, srtrop_stats)

# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}


# add values to table
all_stats$significant <- sapply(X = all_stats$P, FUN = checksig)



# save the stats tables
write.csv(all_stats, file = paste0(outDir, "/TropTemp_Stats.csv"), row.names = FALSE)



##%######################################################%##
#                                                          #
####                       Plots                        ####
#                                                          #
##%######################################################%##

exclQuantiles <- c(0.025,0.975)
# 1. 计算各类别的样本量
category_counts <- table(predictsSites_ab_temp$UI2)
print(category_counts)

# 2. 确定目标样本量（使用最大类别的样本量）
target_sample_size <- max(category_counts)
print(paste("目标样本量:", target_sample_size))

# 3. 创建临时平衡数据集
balanced_data <- predictsSites_ab_temp  # 初始化平衡数据集

# 4. 对每个样本量不足的类别进行重采样
for(category in names(category_counts)) {
  # 如果当前类别样本量已达标，跳过
  if(category_counts[category] >= target_sample_size) {
    next
  }
  
  # 计算需要额外采样的数量
  samples_needed <- target_sample_size - category_counts[category]
  
  # 获取当前类别的所有数据
  category_data <- predictsSites_ab_temp[predictsSites_ab_temp$UI2 == category, ]
  
  # 从当前类别中随机抽样（有放回）
  sampled_rows <- sample(1:nrow(category_data), samples_needed, replace = TRUE)
  additional_samples <- category_data[sampled_rows, ]
  
  # 将抽样结果添加到平衡数据集
  balanced_data <- rbind(balanced_data, additional_samples)
}

# 5. 将平衡后的数据集赋回原变量名
predictsSites_ab_temp <- balanced_data

# 6. 检查重采样后的类别分布
new_counts <- table(predictsSites_ab_temp$UI2)
print(new_counts)


ui2_levels <- levels(factor(predictsSites_ab_temp$UI2))
print(ui2_levels)

# 创建权重映射
weight_values <- c(3, 2, 2, 0.3)
weight_mapping <- setNames(weight_values, ui2_levels)

# 将权重分配给每个观测
weights <- weight_mapping[predictsSites_ab_temp$UI2]

# 检查权重分配是否正确
table(predictsSites_ab_temp$UI2, weights)
print(weights)
# 应用到模型
MeanAbundTemp <- lmer(LogAbund ~ UI2 + poly(StdTmeanAnomalyRS,1) + 
                                UI2:poly(StdTmeanAnomalyRS,1) + 
                                (1|SS) + (1|SSB),
                              data = predictsSites_ab_temp,
                               weights = weights)
r.squaredGLMM(MeanAbundTemp)

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(predictsSites_ab_temp$StdTmeanAnomalyRS),
                        to = max(predictsSites_ab_temp$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low","Urban_Moderate","Urban_High"),
             levels = levels(predictsSites_ab_temp$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2 == "Primary vegetation") & 
                  (nd$StdTmeanAnomaly == nd$StdTmeanAnomaly[which.min(abs(nd$StdTmeanAnomaly))]))
print(refRow)

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = predictsSites_ab_temp$StdTmeanAnomalyRS[
  predictsSites_ab_temp$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = predictsSites_ab_temp$StdTmeanAnomalyRS[
  predictsSites_ab_temp$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = predictsSites_ab_temp$StdTmeanAnomalyRS[
  predictsSites_ab_temp$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = predictsSites_ab_temp$StdTmeanAnomalyRS[
  predictsSites_ab_temp$UI2=="Urban_High"],
  probs = exclQuantiles)

cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")
# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAbundTemp,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,
                       STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Low") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Low") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Moderate") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Moderate") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd$UI2, !is.na(nd$PredMedian)))

p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +  # size改为linewidth
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +  # size改为linewidth
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2), 
        panel.border = element_rect(linewidth = 0.2), 
        axis.ticks = element_line(linewidth = 0.2)) +      
  ggtitle("a Non-tropical Realm")

print(p1)


### 2. MeanAbundTrop ###
#计算每个类别的频率
ui2_counts <- table(predictsSites_ab_trop$UI2)
ui2_freqs <- ui2_counts / sum(ui2_counts)
print(ui2_counts)
# 计算权重为频率的倒数
ui2_weights <- 1 / ui2_freqs
print(ui2_weights)
# 将权重映射到每个观测
weights <- ui2_weights[predictsSites_ab_trop$UI2]

# 标准化权重（使平均权重为1）
weights <- weights / mean(weights)
table(weights)
print(weights)

# 1. 计算各类别的样本量
category_counts <- table(predictsSites_ab_trop$UI2)
print(category_counts)

# 2. 确定目标样本量（使用最大类别的样本量）
target_sample_size <- max(category_counts)
print(paste("目标样本量:", target_sample_size))

# 3. 创建临时平衡数据集
balanced_data <- predictsSites_ab_trop  # 初始化平衡数据集

# 4. 对每个样本量不足的类别进行重采样
for(category in names(category_counts)) {
  # 如果当前类别样本量已达标，跳过
  if(category_counts[category] >= target_sample_size) {
    next
  }
  
  # 计算需要额外采样的数量
  samples_needed <- target_sample_size - category_counts[category]
  
  # 获取当前类别的所有数据
  category_data <- predictsSites_ab_trop[predictsSites_ab_trop$UI2 == category, ]
  
  # 从当前类别中随机抽样（有放回）
  sampled_rows <- sample(1:nrow(category_data), samples_needed, replace = TRUE)
  additional_samples <- category_data[sampled_rows, ]
  
  # 将抽样结果添加到平衡数据集
  balanced_data <- rbind(balanced_data, additional_samples)
}

# 5. 将平衡后的数据集赋回原变量名
predictsSites_ab_trop <- balanced_data

# 6. 检查重采样后的类别分布
new_counts <- table(predictsSites_ab_trop$UI2)
print(new_counts)
# 假设UI2类别是：Primary vegetation, Urban_Low, Urban_Moderate, Urban_High
# 查看类别顺序
ui2_levels <- levels(factor(predictsSites_ab_trop$UI2))
print(ui2_levels)

# 创建权重映射
weight_values <- c(10, 1.5, 1 ,1) 
weight_mapping <- setNames(weight_values, ui2_levels)

# 将权重分配给每个观测
weights <- weight_mapping[predictsSites_ab_trop$UI2]

# 检查权重分配是否正确
table(predictsSites_ab_trop$UI2, weights)
print(weights)

# 应用到模型
MeanAbundTrop <- lmer(Species_richness ~ UI2 + poly(StdTmeanAnomalyRS,1) + 
                               UI2:poly(StdTmeanAnomalyRS,1) + 
                               (1|SS),
                             data = predictsSites_ab_trop,
                             weights = weights)
r.squaredGLMM(MeanAbundTrop)

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(predictsSites_ab_trop$StdTmeanAnomalyRS),
                        to = max(predictsSites_ab_trop$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low","Urban_Moderate","Urban_High"),
             levels = levels(predictsSites_ab_trop$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2 == "Primary vegetation") & 
                  (nd2$StdTmeanAnomaly == nd2$StdTmeanAnomaly[which.min(abs(nd2$StdTmeanAnomaly))]))
print(refRow)
# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = predictsSites_ab_trop$StdTmeanAnomalyRS[
  predictsSites_ab_trop$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = predictsSites_ab_trop$StdTmeanAnomalyRS[
  predictsSites_ab_trop$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = predictsSites_ab_trop$StdTmeanAnomalyRS[
  predictsSites_ab_trop$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = predictsSites_ab_trop$StdTmeanAnomalyRS[
  predictsSites_ab_trop$UI2=="Urban_High"],
  probs = exclQuantiles)
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAbundTrop,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,
                       STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# 定义一个函数来限制值的范围
# limit_range <- function(x, min_val = -100, max_val = 100) {
#   # 检查是否有极端值
#   has_extreme <- any(x < min_val | x > max_val, na.rm = TRUE)
#   
#   if(has_extreme) {
#     # 如果有极端值，进行缩放处理
#     x_no_na <- x[!is.na(x)]
#     
#     # 如果所有值都相同，直接返回限制范围内的值
#     if(length(unique(x_no_na)) == 1) {
#       x[] <- ifelse(x_no_na[1] > 0, max_val, min_val)
#       return(x)
#     }
#     
#     # 获取当前最小值和最大值
#     current_min <- min(x_no_na)
#     current_max <- max(x_no_na)
#     
#     # 线性缩放到目标范围
#     x_scaled <- min_val + (x - current_min) * (max_val - min_val) / (current_max - current_min)
#     return(x_scaled)
#   } else {
#     # 如果没有极端值，则保持原样
#     return(x)
#   }
# }
# 
# # 1. 计算PredMedian并限制范围
# median_values <- (apply(X = a.preds.tmean, MARGIN = 1, FUN = median, na.rm = TRUE) * 100) - 100
# nd2$PredMedian <- limit_range(median_values)
# 
# # 2. 计算PredUpper并限制范围
# upper_values <- (apply(X = a.preds.tmean, MARGIN = 1, FUN = quantile, probs = 0.975, na.rm = TRUE) * 100) - 100
# nd2$PredUpper <- limit_range(upper_values)
# 
# # 3. 计算PredLower并限制范围
# lower_values <- (apply(X = a.preds.tmean, MARGIN = 1, FUN = quantile, probs = 0.025, na.rm = TRUE) * 100) - 100
# nd2$PredLower <- limit_range(lower_values)
# 
# # 4. 按UI2类别检查结果
# tapply(nd2$PredMedian, nd2$UI2, summary)
# tapply(nd2$PredUpper, nd2$UI2, summary)
# tapply(nd2$PredLower, nd2$UI2, summary)

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
tapply(nd2$PredMedian, nd2$UI2, summary)
tapply(nd2$PredUpper, nd2$UI2, summary)
tapply(nd2$PredLower, nd2$UI2, summary)
# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd2$UI2, !is.na(nd2$PredMedian)))

p2 <- ggplot(data = nd21, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  # geom_line(aes(col = UI2), linewidth = 0.75) +  # size改为linewidth
  geom_smooth(aes(col = UI2), linewidth = 0.75, se = FALSE, method = "gam", formula = y ~ s(x, k=5))+
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +  # size改为linewidth
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00","#D55E00")) +
  theme_bw() + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),  # size改为linewidth
        panel.border = element_rect(linewidth = 0.2),      # size改为linewidth
        axis.ticks = element_line(linewidth = 0.2)) +      # size改为linewidth
  ggtitle("c Tropical Realm")

print(p2)


### 3. MeanRichTemp ###
#计算每个类别的频率
ui2_counts <- table(predictsSites_sr_temp$UI2)
ui2_freqs <- ui2_counts / sum(ui2_counts)
print(ui2_counts)
# 计算权重为频率的倒数
ui2_weights <- 1 / ui2_freqs
print(ui2_weights)
# 将权重映射到每个观测
weights <- ui2_weights[predictsSites_sr_temp$UI2]

# 标准化权重（使平均权重为1）
weights <- weights / mean(weights)
table(weights)
# 1. 计算各类别的样本量
category_counts <- table(predictsSites_sr_temp$UI2)
print(category_counts)

# 2. 确定目标样本量（使用最大类别的样本量）
target_sample_size <- max(category_counts)
print(paste("目标样本量:", target_sample_size))

# 3. 创建临时平衡数据集
balanced_data <- predictsSites_sr_temp  # 初始化平衡数据集

# 4. 对每个样本量不足的类别进行重采样
for(category in names(category_counts)) {
  # 如果当前类别样本量已达标，跳过
  if(category_counts[category] >= target_sample_size) {
    next
  }
  
  # 计算需要额外采样的数量
  samples_needed <- target_sample_size - category_counts[category]
  
  # 获取当前类别的所有数据
  category_data <- predictsSites_sr_temp[predictsSites_sr_temp$UI2 == category, ]
  
  # 从当前类别中随机抽样（有放回）
  sampled_rows <- sample(1:nrow(category_data), samples_needed, replace = TRUE)
  additional_samples <- category_data[sampled_rows, ]
  
  # 将抽样结果添加到平衡数据集
  balanced_data <- rbind(balanced_data, additional_samples)
}

# 5. 将平衡后的数据集赋回原变量名
predictsSites_sr_temp <- balanced_data

# 6. 检查重采样后的类别分布
new_counts <- table(predictsSites_sr_temp$UI2)
print(new_counts)

# 查看类别顺序
ui2_levels <- levels(factor(predictsSites_sr_temp$UI2))
print(ui2_levels)

# 创建权重映射
weight_values <- c(3, 2, 2, 0.3)  # 自定义权重
weight_mapping <- setNames(weight_values, ui2_levels)

# 将权重分配给每个观测
weights <- weight_mapping[predictsSites_sr_temp$UI2]

# 检查权重分配是否正确
table(predictsSites_sr_temp$UI2, weights)
print(weights)
#Species_richness ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS, 1) + (1 | SS) + (1 | SSB) + (1 | SSBS)
MeanRichTemp <- lmer(Species_richness~ UI2 + poly(StdTmeanAnomalyRS,1) + 
                               UI2:poly(StdTmeanAnomalyRS,1) + 
                               (1|SS) + (1|SSB),
                             weights = weights,
                             data = predictsSites_sr_temp)
r.squaredGLMM(MeanRichTemp)

nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(predictsSites_sr_temp$StdTmeanAnomalyRS),
                        to = max(predictsSites_sr_temp$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low","Urban_Moderate","Urban_High"),
             levels = levels(predictsSites_sr_temp$UI2)))

# back transform the predictors
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2 == "Primary vegetation") & 
                  (nd3$StdTmeanAnomaly == nd3$StdTmeanAnomaly[which.min(abs(nd3$StdTmeanAnomaly))]))
print(refRow)
# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = predictsSites_sr_temp$StdTmeanAnomalyRS[
  predictsSites_sr_temp$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = predictsSites_sr_temp$StdTmeanAnomalyRS[
  predictsSites_sr_temp$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = predictsSites_sr_temp$StdTmeanAnomalyRS[
  predictsSites_sr_temp$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = predictsSites_sr_temp$StdTmeanAnomalyRS[
  predictsSites_sr_temp$UI2=="Urban_High"],
  probs = exclQuantiles)
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanRichTemp, data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Low") & (nd3$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Low") & (nd3$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Moderate") & (nd3$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Moderate") & (nd3$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd3$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd3$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
tapply(nd3$PredMedian, nd3$UI2, summary)
# tapply(nd3$PredUpper, nd3$UI2, summary)
# tapply(nd3$PredLower, nd3$UI2, summary)
# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd3$UI2, !is.na(nd3$PredMedian)))

p3 <- ggplot(data = nd31, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), linewidth = 0.75) +  # size改为linewidth
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +  # size改为linewidth
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00","#D55E00")) +
  theme_bw() + 
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(-100, 150) +
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),  # size改为linewidth
        panel.border = element_rect(linewidth = 0.2),      # size改为linewidth
        axis.ticks = element_line(linewidth = 0.2)) +      # size改为linewidth
  ggtitle("b Non-tropical Realm")


print(p3)

### 4. MeanRichTrop ###
# Species_richness ~ UI2 + UI2:poly(StdTmeanAnomalyRS, 1) + poly(StdTmeanAnomalyRS,  1) + (1 | SS) + (1 | SSB) + (1 | SSBS)
MeanRichTrop <- lmer(Species_richness ~ UI2 + poly(StdTmeanAnomalyRS,1) + 
                              UI2:poly(StdTmeanAnomalyRS,1) + 
                              (1|SS) + (1|SSB),
                            data = predictsSites_sr_trop)

r.squaredGLMM(MeanRichTrop)

nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(predictsSites_sr_trop$StdTmeanAnomalyRS),
                        to = max(predictsSites_sr_trop$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low","Urban_moderate","Urban_High"),
             levels = levels(predictsSites_sr_trop$UI2)))

# back transform the predictors
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd4$LogAbund <- 0
nd4$Species_richness <- 0

nd4$SS <- NA
nd4$SSB <- NA 
# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd4$UI2 == "Primary vegetation") & 
                  (nd4$StdTmeanAnomaly == nd4$StdTmeanAnomaly[which.min(abs(nd4$StdTmeanAnomaly))]))
print(refRow)
# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = predictsSites_sr_trop$StdTmeanAnomalyRS[
  predictsSites_sr_trop$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = predictsSites_sr_trop$StdTmeanAnomalyRS[
  predictsSites_sr_trop$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = predictsSites_sr_trop$StdTmeanAnomalyRS[
  predictsSites_sr_trop$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = predictsSites_sr_trop$StdTmeanAnomalyRS[
  predictsSites_sr_trop$UI2=="Urban_High"],
  probs = exclQuantiles)
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")


# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanRichTrop,data = nd4)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban_Low") & (nd4$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban_Low") & (nd4$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban_Moderate") & (nd4$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban_Moderate") & (nd4$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban_High") & (nd4$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Urban_High") & (nd4$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd4$UI2, !is.na(nd4$PredMedian)))

p4 <- ggplot(data = nd41, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  ylab("Change in Species Richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        # legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) +  
  ggtitle("d Tropical Realm")

print(p4)

## arrange Mean anomaly plots

leg <- get_legend(p4)

# cowplot::plot_grid(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
#                                       p3+ theme(legend.position = "none"), 
#                                       p2 + theme(legend.position = "none"),
#           p4 + theme(legend.position = "none"),
#           nrow = 2),
#           leg, nrow= 2, rel_heights = c(5,1))
combined_plots <- cowplot::plot_grid(
  p1 + theme(legend.position = "none"), 
  p3 + theme(legend.position = "none"), 
  p2 + theme(legend.position = "none"), 
  p4 + theme(legend.position = "none"),
  nrow = 2
)

# 将图形和图例组合，图例对齐到右侧的中间
final_plot <- cowplot::plot_grid(combined_plots, leg, nrow = 1, 
                                 rel_widths = c(1, 0.4), 
                                 align = 'v')

# 调整图例的位置为垂直居中
final_plot <- cowplot::ggdraw(final_plot) + 
  draw_plot(leg, x = 0.95, y = 0.5, width = 0.3, height = 0.5, hjust = 1, vjust = 0.5)

# 显示最终图形
print(final_plot)

ggsave(filename = paste0(outDir, "Figure3_MeanAnom_TempTrop.pdf"), plot = last_plot(), width = 180, height = 170, units = "mm", dpi = 300)



##%######################################################%##
#                                                          #
####      Run separate models for the max anomaly       ####
#                                                          #
##%######################################################%##




# 1. Temperate, Abundance, max anomaly
MaxAbundTemp <- GLMERSelect(modelData = predictsSites_ab_temp,responseVar = "LogAbund",
                             fitFamily = "gaussian",fixedFactors = "UI2",
                             fixedTerms = list(StdTmaxAnomalyRS=1),
                             randomStruct = "(1|SS)+(1|SSB)",
                             fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                             saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MaxAbundTemp$model)
# LogAbund ~ UI2 + (1 | SS) + (1 | SSB)

# 2. Tropical, Abundance, max anomaly
MaxAbundTrop <- GLMERSelect(modelData = predictsSites_ab_trop,responseVar = "LogAbund",
                             fitFamily = "gaussian",fixedFactors = "UI2",
                             fixedTerms = list(StdTmaxAnomalyRS=1),
                             randomStruct = "(1|SS)+(1|SSB)",
                             fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                             saveVars = c("Species_richness", "Total_abundance", "SSBS"))

summary(MaxAbundTrop$model)
# LogAbund ~ UI2 + poly(StdTmaxAnomalyRS, 1) + UI2:poly(StdTmaxAnomalyRS,  1) + (1 | SS) + (1 |   SSB)


# 3. Temperate, Richness, max anomaly
MaxRichTemp <- GLMERSelect(modelData = predictsSites_sr_temp,responseVar = "Species_richness",
                            fitFamily = "poisson",fixedFactors = "UI2",
                            fixedTerms = list(StdTmaxAnomalyRS=1),
                            randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                            fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                            saveVars = c("Total_abundance", "SSBS"))

summary(MaxRichTemp$model)
# Species_richness ~ UI2 + UI2:poly(StdTmaxAnomalyRS, 1) + poly(StdTmaxAnomalyRS,      1) + (1 | SS) + (1 | SSB) + (1 | SSBS)

# 4. Tropical, Richness, max anomaly
MaxRichTrop <- GLMERSelect(modelData = predictsSites_sr_trop,responseVar = "Species_richness",
                            fitFamily = "poisson",fixedFactors = "UI2",
                            fixedTerms = list(StdTmaxAnomalyRS=1),
                            randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                            fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                            saveVars = c("Total_abundance", "SSBS"))

summary(MaxRichTrop$model)
# Species_richness ~ UI2 + UI2:poly(StdTmaxAnomalyRS, 1) + poly(StdTmaxAnomalyRS,1) + (1 | SS) + (1 | SSB) + (1 | SSBS)


save(MaxAbundTemp, file = paste0(outDir, "/MaxAbundTemp_output.rdata"))
save(MaxAbundTrop, file = paste0(outDir, "/MaxAbundTrop_output.rdata"))
save(MaxRichTemp, file = paste0(outDir, "/MaxRichTemp_output.rdata"))
save(MaxRichTrop, file = paste0(outDir, "/MaxRichTrop_output.rdata"))

# load(file = paste0(outDir, "/MaxAbundTemp_output.rdata"))
# load(file = paste0(outDir, "/MaxAbundTrop_output.rdata"))
# load(file = paste0(outDir, "/MaxRichTemp_output.rdata"))
# load(file = paste0(outDir, "/MaxRichTrop_output.rdata"))


#### save model output tables using sjPlot package ####

sjPlot::tab_model(MaxAbundTemp$model, MaxAbundTrop$model, transform = NULL, file = paste0(outDir, "/AbunMaxAnomTempTrop_output_table.html"))
summary(MaxAbundTemp$model)
R2GLMER(MaxAbundTemp$model)
summary(MaxAbundTrop$model)
R2GLMER(MaxAbundTrop$model)

tab_model(MaxRichTemp$model, MaxRichTrop$model, transform = NULL, file = paste0(outDir, "/RichMaxAnomTempTrop_output_table.html"))
summary(MaxRichTemp$model)
R2GLMER(MaxRichTemp$model) # use these values
# $conditional
#[1] 0.6597703
#$marginal
#[1] 0.004763959
summary(MaxRichTrop$model)
R2GLMER(MaxRichTrop$model) # use these values
# $conditional
# [1] 0.552677
# $marginal
# [1] 0.02567178

# save model stats

# save the stats info
abtemp_stats <- as.data.frame(MaxAbundTemp$stats)
abtrop_stats <- as.data.frame(MaxAbundTrop$stats)
srtemp_stats <- as.data.frame(MaxRichTemp$stats)
srtrop_stats <- as.data.frame(MaxRichTrop$stats)

abtemp_stats$significant <- NA
abtrop_stats$significant <- NA
srtemp_stats$significant <- NA
srtrop_stats$significant <- NA

abtemp_stats$model <- "abtemp"
abtrop_stats$model <- "abtrop"
srtemp_stats$model <- "srtemp"
srtrop_stats$model <- "srtrop"

all_stats <- rbind(abtemp_stats, abtrop_stats, srtemp_stats, srtrop_stats)

# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}


# add values to table
all_stats$significant <- sapply(X = all_stats$P, FUN = checksig)



# save the stats tables
write.csv(all_stats, file = paste0(outDir, "/TropTemp_Max_Stats.csv"), row.names = FALSE)




##%######################################################%##
#                                                          #
####                  Max anom Plots                    ####
#                                                          #
##%######################################################%##

exclQuantiles <- c(0.025,0.975)


### 1. MaxAbundTemp ###


nd <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAbundTemp$data$StdTmaxAnomalyRS),
                        to = max(MaxAbundTemp$data$StdTmaxAnomalyRS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAbundTemp$data$UI2)))

# back transform the predictors
nd$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2 == "Primary vegetation") & 
                  (nd$StdTmaxAnomaly == nd$StdTmaxAnomaly[which.min(abs(nd$StdTmaxAnomaly))]))
print(refRow)
# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAbundTemp$data$StdTmaxAnomalyRS[
  MaxAbundTemp$data$UI2=="Urban_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxAbundTemp$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Low") & (nd$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Low") & (nd$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_High") & (nd$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_High") & (nd$StdTmaxAnomalyRS > QAL[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

#"#009E73" - green
#"#0072B2" - blue
#"#E69F00" - yellow
#"#D55E00" - red

p1 <- ggplot(data = nd, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00")) +
  theme_bw() + 
  ylab("Change in Total Abundance (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-2, 1)) +
  ylim(c(-200, 200)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("a Non-tropical Realm")

print(p1)

### 2. MaxAbundTrop ###


nd2 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAbundTrop$data$StdTmaxAnomalyRS),
                        to = max(MaxAbundTrop$data$StdTmaxAnomalyRS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAbundTrop$data$UI2)))

# back transform the predictors
nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAbundTrop$data$StdTmaxAnomalyRS[
  MaxAbundTrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxAbundTrop$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p2 <- ggplot(data = nd2, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  ylab("Change in Total Abundance (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 100)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("b Tropical Realm")





### 3. MaxRichTemp ###


nd3 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxRichTemp$data$StdTmaxAnomalyRS),
                        to = max(MaxRichTemp$data$StdTmaxAnomalyRS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxRichTemp$data$UI2)))

# back transform the predictors
nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomaly==min(abs(nd3$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxRichTemp$data$StdTmaxAnomalyRS[
  MaxRichTemp$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxRichTemp$model,data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p3 <- ggplot(data = nd3, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  ylab("Change in Species Richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(c(-100, 100)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("c Non-tropical Realm")



### 4. MaxRichTrop ###


nd4 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxRichTrop$data$StdTmaxAnomalyRS),
                        to = max(MaxRichTrop$data$StdTmaxAnomalyRS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxRichTrop$data$UI2)))

# back transform the predictors
nd4$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxRichTrop$data$StdTmaxAnomalyRS[
  MaxRichTrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxRichTrop$model,data = nd4)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


p4 <- ggplot(data = nd4, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  ylab("Change in Species Richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 100)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank()) + 
  ggtitle("d Tropical Realm")



## arrange Mean anomaly plots

leg <- get_legend(p4)

cowplot::plot_grid(cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                    p2 + theme(legend.position = "none"), 
                    p3+ theme(legend.position = "none"), 
                    p4 + theme(legend.position = "none"),
                    nrow = 2),
          leg, nrow= 2, rel_heights = c(5,1))


ggsave(file = paste0(outDir, "/Extended_Data5_MaxAnom_TempTrop.pdf"), width = 8, height = 8.5)
ggsave(filename = paste0(outDir, "Extended_Data5_MaxAnom_TempTrop.jpeg"), plot = last_plot(), width = 183, height = 200, units = "mm", dpi = 300)






### models including realm - not used due to high VIF factors between variables


# # 1. Abundance, mean anomaly
# MeanAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
#                                      fitFamily = "gaussian",
#                                      fixedFactors = c("UI2", "Tropical"),
#                                      fixedTerms = list(StdTmeanAnomalyRS=1),
#                                      randomStruct = "(1|SS)+(1|SSB)",
#                                      fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)",
#                                                            "UI2:Tropical",
#                                                            "Tropical:poly(StdTmeanAnomalyRS,1)",
#                                                            "Tropical:poly(StdTmeanAnomalyRS,1):UI2"),
#                                      saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MeanAnomalyModelAbund$model)
# 
# # model is rank deficient
# 
# # selected model:
# # LogAbund ~ UI2 + poly(StdTmeanAnomalyRS, 1) + Tropical +
# # UI2:Tropical +
# # Tropical:poly(StdTmeanAnomalyRS, 1):UI2 + 
# # (1 | SS) + (1 | SSB)
# 
# 
# check_collinearity(MeanAnomalyModelAbund$model)
# 
# # the 3 way interaction between land use, climate and realm has a high VIF
# 
# 
# 
# # save the model output
# save(MeanAnomalyModelAbund, file = paste0(outDir, "/MeanAnomalyModelAbund_3way.rdata"))
# 
# # 2. Richness, mean anomaly
# MeanAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
#                                     fitFamily = "poisson",
#                                     fixedFactors = c("UI2", "Tropical"),
#                                     fixedTerms = list(StdTmeanAnomalyRS=1),
#                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
#                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)",
#                                                           "UI2:Tropical",
#                                                           "Tropical:poly(StdTmeanAnomalyRS,1)",
#                                                           "Tropical:poly(StdTmeanAnomalyRS,1):UI2"),
#                                     saveVars = c("Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MeanAnomalyModelRich$model)
# 
# # selected model:
# # Species_richness ~ UI2 + Tropical + poly(StdTmeanAnomalyRS, 1) +  
# # UI2:poly(StdTmeanAnomalyRS, 1) + UI2:Tropical +  
# # Tropical:poly(StdTmeanAnomalyRS, 1):UI2 + 
# # (1 | SS) +  (1 | SSB) + (1 | SSBS)
# 
# # Model failed to converge: degenerate  Hessian with 2 negative eigenvalues
# 
# check_collinearity(MeanAnomalyModelAbund$model)
# 
# 
# # the 3 way interaction between land use, climate and realm has a high VIF
# 
# 
# # save model output
# save(MeanAnomalyModelRich, file = paste0(outDir, "/MeanAnomalyModelRich_3way.rdata"))
# 
# 
# # 3. Abundance, max anomaly
# MaxAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
#                                     fitFamily = "gaussian",
#                                     fixedFactors = c("UI2", "Tropical"),
#                                     fixedTerms = list(StdTmaxAnomalyRS=1),
#                                     randomStruct = "(1|SS)+(1|SSB)",
#                                     fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)",
#                                                           "UI2:Tropical",
#                                                           "Tropical:poly(StdTmaxAnomalyRS,1)",
#                                                           "Tropical:poly(StdTmaxAnomalyRS,1):UI2"),
#                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MaxAnomalyModelAbund$model)
# 
# # selected model:
# # LogAbund ~ UI2 + UI2:Tropical + Tropical + (1 | SS) + (1 | SSB)
# 
# # save model output
# save(MaxAnomalyModelAbund, file = paste0(outDir, "/MaxAnomalyModelAbund_3way.rdata"))
# 
# 
# # 4. Richness, max anomaly
# MaxAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
#                                    fitFamily = "poisson",
#                                    fixedFactors = c("UI2", "Tropical"),
#                                    fixedTerms = list(StdTmaxAnomalyRS=1),
#                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
#                                    fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)",
#                                                          "UI2:Tropical",
#                                                          "Tropical:poly(StdTmaxAnomalyRS,1)",
#                                                          "Tropical:poly(StdTmaxAnomalyRS,1):UI2"),  
#                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))
# 
# summary(MaxAnomalyModelRich$model)
# 
# # selected model:
# # Species_richness ~ UI2 + Tropical + poly(StdTmaxAnomalyRS, 1) + 
# # UI2:poly(StdTmaxAnomalyRS, 1) + UI2:Tropical + Tropical:poly(StdTmaxAnomalyRS, 1):UI2 +  
# # (1 | SS) + (1 | SSB) + (1 | SSBS)
# 
# 
# # Model failed to converge with max|grad| = 0.00613945 (tol = 0.001, component 1)
# 
# 
# 
# # save model output
# save(MaxAnomalyModelRich, file = paste0(outDir, "/MaxAnomalyModelRich_3way.rdata"))
# 
# 
# #load(paste0(outDir, "MeanAnomalyModelAbund_3way.rdata"))
# #load(paste0(outDir, "/MeanAnomalyModelRich_3way.rdata"))
# #load(paste0(outDir, "/MaxAnomalyModelAbund_3way.rdata"))
# #load(paste0(outDir, "/MaxAnomalyModelRich_3way.rdata"))
# 
# 
# summary(MeanAnomalyModelAbund$model)
# summary(MeanAnomalyModelRich$model)
# summary(MaxAnomalyModelAbund$model)
# summary(MaxAnomalyModelRich$model)
# 


##%######################################################%##
#                                                          #
####                  Looking at VIFs                   ####
#                                                          #
##%######################################################%##

# 
# # remove all NAs from the dataset
# 
# predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]
# predictsSites <- predictsSites[!is.na(predictsSites$UI2), ] # 6069
# 
# 
# effs <- predictsSites[ , c("UI2", "StdTmeanAnomalyRS", "StdTmaxAnomalyRS", "NH_5000.rs", "Tropical")]
# 
# 
# x <- cor(effs[, 2:4]) %>% det() # just the continuous variables
# 
# 
# effdata <- predictsSites[ , c("LogAbund", "UI2", "StdTmeanAnomalyRS", "StdTmaxAnomalyRS", "NH_5000.rs", "Tropical")]
# effdata <- effdata[!is.na(effdata$LogAbund), ]
# 
# fitmod <- lm(LogAbund ~ UI2 + StdTmeanAnomalyRS + Tropical, data = effdata)
# stepVIF(fitmod, threshold = 10, verbose = TRUE)
# 
# 
# car::vif(MeanAnomalyModelAbund$model)
# car::vif(fitmod)
# 
# 
# res <- stepwise.vif(dataset = effdata[, c(2:3, 6)], 
#                     metrics = c("UI2", "StdTmeanAnomalyRS", "Tropical"), 
#                     vif.threshold = 5,
#                     verbose = T)
# 
# check_collinearity(MeanAnomalyModelAbund$model)
# check_collinearity(MeanAnomalyModelRich$model)
# check_collinearity(MaxAnomalyModelAbund$model)
# check_collinearity(MaxAnomalyModelRich$model)
# 


# ############## plotting ##############
# 
# 
# ## Realm/anomaly/LU interactions
# 
# 
# 
# ### 1. Abundance, mean anomaly
# 
# nd <- expand.grid(
#   StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
#                         to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
#                         length.out = 100),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(MeanAnomalyModelAbund$data$UI2)),
#   Tropical = factor(c("Tropical", "Temperate"),
#                     levels = levels(MeanAnomalyModelAbund$data$Tropical)))
# 
# # back transform the predictors
# nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd$StdTmeanAnomalyRS,
#   originalX = predictsSites$StdTmeanAnomaly)
# 
# # set richness and abundance to 0 - to be predicted
# nd$LogAbund <- 0
# nd$Species_richness <- 0
# 
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# refRow1 <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) & (nd$Tropical == "Tropical"))
# refRow2 <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) & (nd$Tropical == "Temperate"))
# refRow2 <- refRow2 -400
# 
# # adjust plot 1: mean anomaly and abundance
# 
# exclQuantiles <- c(0.025,0.975)
# 
# 
# QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelAbund$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# 
# # predict the results
# a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)
# 
# # back transform the abundance values
# a.preds.tmean <- exp(a.preds.tmean)-0.01
# 
# 
# a.preds.tmean_t <- a.preds.tmean[1:400, ] # tropical set
# a.preds.tmean <- a.preds.tmean[401:800, ] # temperate set
# 
# # convert to relative to reference
# a.preds.tmean_t <- sweep(x = a.preds.tmean_t,MARGIN = 2,STATS = a.preds.tmean_t[refRow1,],FUN = '/') # tropical
# a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow2,],FUN = '/') # temperate
# 
# # remove anything above and below the quantiles
# a.preds.tmean_t[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1]) & nd$Tropical == "Tropical"),] <- NA
# a.preds.tmean_t[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2]) & nd$Tropical == "Tropical"),] <- NA
# 
# # remove anything above and below the quantiles
# a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1]) & nd$Tropical == "Temperate")-400,] <- NA
# a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2]) & nd$Tropical == "Temperate")-400,] <- NA
# 
# # Get the median, upper and lower quants for the plot
# nd$PredMedian[1:400] <- ((apply(X = a.preds.tmean_t,MARGIN = 1,
#                          FUN = median,na.rm=TRUE))*100)-100
# nd$PredUpper[1:400] <- ((apply(X = a.preds.tmean_t,MARGIN = 1,
#                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd$PredLower[1:400] <- ((apply(X = a.preds.tmean_t,MARGIN = 1,
#                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd$PredMedian[401:800] <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                          FUN = median,na.rm=TRUE))*100)-100
# nd$PredUpper[401:800] <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd$PredLower[401:800] <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# 
# plot_data <- nd[nd$Tropical == "Temperate",]
# 
# p1 <- ggplot(data = plot_data, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot_data$PredLower, ymax = plot_data$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in total abundance (%)") +
#   xlab("Standardised Temperature Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 750)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12)) +
#   ggtitle("Temperate")
# 
# 
# plot_data2 <- nd[nd$Tropical == "Tropical",]
# 
# p2 <- ggplot(data = plot_data2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot_data2$PredLower, ymax = plot_data2$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in total abundance (%)") +
#   xlab("Standardised Climate Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12), legend.position = "none") +
#   ggtitle("Tropical")
# 
# 
# legend <- get_legend(p1)
# p3 <- cowplot::plot_grid(p1+theme(legend.position = "none"), p2, legend, ncol = 3)
# 
# 
# ### 2. SR mean anomaly
# 
# nd2 <- expand.grid(
#   StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
#                         to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
#                         length.out = 100),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(MeanAnomalyModelRich$data$UI2)),
#   Tropical = factor(c("Tropical", "Temperate"),
#                     levels = levels(MeanAnomalyModelRich$data$Tropical)))
# 
# # back transform the predictors
# nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd2$StdTmeanAnomalyRS,
#   originalX = predictsSites$StdTmeanAnomaly)
# 
# # set richness and abundance to 0 - to be predicted
# nd2$LogAbund <- 0
# nd2$Species_richness <- 0
# 
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# refRow1 <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & (nd2$Tropical == "Tropical"))
# refRow2 <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & (nd2$Tropical == "Temperate"))
# refRow2 <- refRow2 -400
# 
# 
# # quantiles of data to show on the plot
# exclQuantiles <- c(0.025,0.975)
# 
# 
# QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
#   MeanAnomalyModelRich$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# # predict the results
# sr.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd2)
# 
# # back transform the abundance values
# sr.preds.tmean <- exp(sr.preds.tmean)
# 
# 
# sr.preds.tmean_t <- sr.preds.tmean[1:400, ] # tropical set
# sr.preds.tmean <- sr.preds.tmean[401:800, ] # temperate set
# 
# # convert to relative to reference
# sr.preds.tmean_t <- sweep(x = sr.preds.tmean_t,MARGIN = 2,STATS = sr.preds.tmean_t[refRow1,],FUN = '/')
# sr.preds.tmean <- sweep(x = sr.preds.tmean,MARGIN = 2,STATS = sr.preds.tmean[refRow2,],FUN = '/')
# 
# # remove anything above and below the quantiles
# sr.preds.tmean_t[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1]) & nd2$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2]) & nd2$Tropical == "Tropical"),] <- NA
# 
# # remove anything above and below the quantiles
# sr.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1]) & nd2$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2]) & nd2$Tropical == "Temperate")-400,] <- NA
# 
# # Get the median, upper and lower quants for the plot
# nd2$PredMedian[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                 FUN = median,na.rm=TRUE))*100)-100
# nd2$PredUpper[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd2$PredMedian[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                   FUN = median,nsr.rm=TRUE))*100)-100
# nd2$PredUpper[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                  FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                  FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# 
# 
# plotdata <- nd2[nd2$Tropical == "Temperate", ]
#   
# p4 <- ggplot(data = plotdata, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plotdata$PredLower, ymax = plotdata$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in species richness (%)") +
#   xlab("Standardised Climate Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 950)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12))+
#   ggtitle("Temperate")
# 
# 
# 
# plotdata2 <- nd2[nd2$Tropical == "Tropical", ]
# 
# p5 <- ggplot(data = plotdata2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plotdata2$PredLower, ymax = plotdata2$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in species richness (%)") +
#   xlab("Standardised Climate Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12), legend.position = "none")+
#   ggtitle("Tropical")
# 
# legend2 <- get_legend(p4)
# p6 <- cowplot::plot_grid(p4+theme(legend.position = "none"), p5, legend, ncol = 3)
# 



### 3. Abundance, max anomaly

# 3-way interaction not selected for this set

#nd3 <- expand.grid(
#  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
#                        to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
#                        length.out = 100),
#  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#             levels = levels(MaxAnomalyModelAbund$data$UI2)),
#  Tropical = factor(c("Tropical", "Temperate"),
#                    levels = levels(MaxAnomalyModelAbund$data$Tropical)))

# back transform the predictors
#nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
#  transformedX = nd3$StdTmaxAnomalyRS,
#  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
#nd3$LogAbund <- 0
#nd3$Species_richness <- 0


# reference for % difference = primary vegetation and positive anomaly closest to 0
#refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomaly==min(abs(nd3$StdTmaxAnomaly))) & (nd3$Tropical == "Temperate"))

# quantiles of data to show on the plot
#exclQuantiles <- c(0.025,0.975)


#QPV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Primary vegetation"],
#  probs = exclQuantiles)
#QSV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Secondary vegetation"],
#  probs = exclQuantiles)
#QAL <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Agriculture_Low"],
#  probs = exclQuantiles)
#QAH <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
#  MaxAnomalyModelAbund$data$UI2=="Agriculture_High"],
#  probs = exclQuantiles)

# predict the results
#a.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelAbund$model,data = nd3)

# back transform the abundance values
#a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
#a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
#a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS < QPV[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS > QPV[2])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS < QSV[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS > QSV[2])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS < QAL[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS > QAL[2])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS < QAH[1])),] <- NA
#a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
#nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                         FUN = median,na.rm=TRUE))*100)-100
#nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
#nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
#                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


#nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

#p3 <- ggplot(data = nd3, aes(x = StdTmaxAnomaly, y = PredMedian)) +
#  geom_line(aes(col = UI2), size = 1) +
#  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = UI2), alpha = 0.2) +
#  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#  facet_wrap(~Tropical, ncol = 2) +
#  theme_bw() +
#  labs(fill = "% NH", col = "% NH") +
#  ylab("Total Abundance (%)") +
#  xlab("Standardised Climate Anomaly Maximum") +
#  xlim(c(-0.5, 2)) +
#  ylim(c(-100, 150)) +
#  theme(aspect.ratio = 1, text = element_text(size = 12))


# 
# 
# 
# ### 4. Species richness, max anomaly
# 
# 
# nd4 <- expand.grid(
#   StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
#                        to = max(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
#                        length.out = 100),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(MaxAnomalyModelRich$data$UI2)),
#   Tropical = factor(c("Tropical", "Temperate"),
#                     levels = levels(MaxAnomalyModelRich$data$Tropical)))
# 
# # back transform the predictors
# nd4$StdTmaxAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd4$StdTmaxAnomalyRS,
#   originalX = predictsSites$StdTmaxAnomaly)
# 
# # set richness and abundance to 0 - to be predicted
# nd4$LogAbund <- 0
# nd4$Species_richness <- 0
# 
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# refRow1 <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))) & (nd4$Tropical == "Tropical"))
# refRow2 <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))) & (nd4$Tropical == "Temperate"))
# refRow2 <- refRow2 -400
# 
# # quantiles of data to show on the plot
# exclQuantiles <- c(0.025,0.975)
# 
# 
# QPV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
#   MaxAnomalyModelRich$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# # predict the results
# sr.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model,data = nd4)
# 
# # back transform the abundance values
# sr.preds.tmean <- exp(sr.preds.tmean)
# 
# 
# sr.preds.tmean_t <- sr.preds.tmean[1:400, ] # tropical set
# sr.preds.tmean <- sr.preds.tmean[401:800, ] # temperate set
# 
# # convert to relative to reference
# sr.preds.tmean_t <- sweep(x = sr.preds.tmean_t,MARGIN = 2,STATS = sr.preds.tmean_t[refRow1,],FUN = '/')
# sr.preds.tmean <- sweep(x = sr.preds.tmean,MARGIN = 2,STATS = sr.preds.tmean[refRow2,],FUN = '/')
# 
# # remove anything above and below the quantiles
# sr.preds.tmean_t[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS < QPV[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS > QPV[2]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS < QSV[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS > QSV[2]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS < QAL[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS > QAL[2]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS < QAH[1]) & nd4$Tropical == "Tropical"),] <- NA
# sr.preds.tmean_t[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS > QAH[2]) & nd4$Tropical == "Tropical"),] <- NA
# 
# # remove anything above and below the quantiles
# sr.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS < QPV[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS > QPV[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS < QSV[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS > QSV[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS < QAL[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS > QAL[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS < QAH[1]) & nd4$Tropical == "Temperate")-400,] <- NA
# sr.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS > QAH[2]) & nd4$Tropical == "Temperate")-400,] <- NA
# 
# # Get the median, upper and lower quants for the plot
# nd4$PredMedian[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                  FUN = median,na.rm=TRUE))*100)-100
# nd4$PredUpper[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd4$PredLower[1:400] <- ((apply(X = sr.preds.tmean_t,MARGIN = 1,
#                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd4$PredMedian[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                    FUN = median,na.rm=TRUE))*100)-100
# nd4$PredUpper[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd4$PredLower[401:800] <- ((apply(X = sr.preds.tmean,MARGIN = 1,
#                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# 
# nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# 
# plot.data <- nd4[nd4$Tropical == "Temperate",]
# 
# p7 <- ggplot(data = plot.data, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot.data$PredLower, ymax = plot.data$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Species Richness (%)") +
#   xlab("Standardised Climate \nAnomaly Maximum") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12))+
#   ggtitle("Temperate")
# 
# plot.data2 <- nd4[nd4$Tropical =="Tropical",]
# 
# p8 <- ggplot(data = plot.data2, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = UI2), size = 1) +
#   geom_ribbon(aes(ymin = plot.data2$PredLower, ymax = plot.data2$PredUpper, fill = UI2), alpha = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
#   #facet_wrap(~Tropical, ncol = 2) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Species Richness (%)") +
#   xlab("Standardised Climate \nAnomaly Maximum") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 100)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12), legend.position = "none") + 
#   ggtitle("Tropical")
# 
# #"#009E73" - green
# #"#0072B2" - blue
# #"#E69F00" - yellow
# #"#D55E00" - red
# 
# legend3 <- get_legend(p7)
# p9 <- plot_grid(p7+theme(legend.position = "none"), p8, legend, ncol = 3)
# 
# 
# 
# # organise plots into one document
# 
# plot_grid(p3, p6,  p9, ncol = 1, 
#           labels = c("A", "B", "C"), label_size = 12, rel_widths = c(1,1,0.5))
# 
# # save plot
# ggsave(filename = paste0(outDir, "Plots_climate_LU_Tropical_ALL_allinteractions.pdf"), width = 9, height = 9, units = "in")
# 
# 
# 



t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()

