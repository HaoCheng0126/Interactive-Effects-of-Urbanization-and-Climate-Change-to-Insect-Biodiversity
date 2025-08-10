##%######################################################%##
#                                                          #
####             Testing buffer zone sizes              ####
####             for natural habitat metric             ####
#                                                          #
##%######################################################%##

# In this script, we test the difference between model results when different 
# buffer sizes are used to determine amount of nearby natural habitat. 


rm(list = ls())


# directories
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "14_Additional_Tests/"
if(!dir.exists(outDir)) dir.create(outDir)


# sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)


library(StatisticalModels)
library(ggplot2)
library(cowplot)
source('Functions.R')

###Create Models for all insects in predicts for Standardised climate anomaly and Land interactions

predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))


#### Abundance models ####

modelData <- predictsSites[!is.na(predictsSites$LogAbund), ] 
modelData <- modelData[!is.na(modelData$StdTmeanAnomalyRS), ] # 1984
nrow(modelData)
AbundMeanAnomalyModel_1000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_1000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))
AbundMeanAnomalyModel_3000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_3000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))
AbundMeanAnomalyModel_5000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_5000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))
AbundMeanAnomalyModel_10000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_10000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))


### models for the stats


MeanAbun_1000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, GS_1000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(GS_1000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))


test <- GLMERSelect(modelData = modelData,
            responseVar = "LogAbund",
            fitFamily = "gaussian",fixedFactors = "UI2",
            fixedTerms = list(StdTmeanAnomalyRS=1, GS_1000.rs=1),
            randomStruct = "(1|SS)+(1|SSB)",
            fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(GS_1000.rs,1)"))

test$model


save(MeanAbun_1000, file = paste0(outDir, "MeanAbun_buffer_1000.rdata"))

MeanAbun_3000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, GS_3000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(GS_3000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))

save(MeanAbun_3000, file = paste0(outDir, "MeanAbun_buffer_3000.rdata"))

MeanAbun_5000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, GS_5000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(GS_5000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))

save(MeanAbun_5000, file = paste0(outDir, "MeanAbun_buffer_5000.rdata"))


MeanAbun_10000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, GS_10000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(GS_10000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))

save(MeanAbun_10000, file = paste0(outDir, "MeanAbun_buffer_10000.rdata"))


# take a look at the results

MeanAbun_1000$stats
MeanAbun_3000$stats
MeanAbun_5000$stats
MeanAbun_10000$stats

#### Species richness models ####

modelData <- predictsSites # 2170

RichMeanAnomalyModel_1000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_1000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

RichMeanAnomalyModel_3000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_3000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

RichMeanAnomalyModel_5000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_5000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

RichMeanAnomalyModel_10000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * GS_10000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))


# take a look at the results
fixef(RichMeanAnomalyModel_1000$model)
fixef(RichMeanAnomalyModel_3000$model)
fixef(RichMeanAnomalyModel_5000$model)
fixef(RichMeanAnomalyModel_10000$model)




#### plots of abundance results ####

### 1000 buffer

AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_1000

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low","Urban_Moderate","Urban_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  GS_1000.rs=c(-0.96,0.01493712,1.00,2.0))
# GS_1000.rs=AbundMeanAnomalyModel1$data$GS_1000.rs)
# GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd$GS_1000 <- round(BackTransformCentreredPredictor(
  transformedX = nd$GS_1000.rs,originalX = predictsSites$GS_1000)*100,0)

nd$GS_1000[nd$GS_1000 == 99] <- 100

# 定义区间和对应的标签
breaks <- c(-Inf, 25, 55, 83, 120)  # 划分的区间
labels <- c(25, 50, 75, 100)         # 区间对应的标签

# 使用 cut 函数将 GS_1000 划分为区间，并转换为数字型
nd$GS_1000 <- cut(nd$GS_1000, 
                  breaks = breaks, 
                  labels = labels, 
                  right = TRUE,  # 右开区间
                  include.lowest = TRUE)  # 包含下限

# 将切割后的类别转化为数值型（如果需要）
nd$GS_1000 <- as.numeric(as.character(nd$GS_1000))

# set values for richness and Abundance
nd$LogAbund <- 0
nd$Species_richness <- 0

# set the reference row
refRow <- which((nd$UI2=="Primary vegetation") & 
                  nd$StdTmeanAnomaly== min(nd$StdTmeanAnomaly[nd$StdTmeanAnomaly > 0]) &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd$GS_1000==100))
refRow <- which((nd$UI2 == "Primary vegetation") & 
                  nd$StdTmeanAnomaly == nd$StdTmeanAnomaly[which.min(abs(nd$StdTmeanAnomaly))]&
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd$GS_1000==100))[1]

print(refRow)
# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_High"],
  probs = exclQuantiles)
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")
# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Low") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Low") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Moderate") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_Moderate") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Urban_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs
# 检查是否有有效的预测值
cat("Number of valid predictions:
")
print(table(nd$UI2, !is.na(nd$PredMedian)))
tapply(nd$PredMedian, nd$UI2, summary)
tapply(nd$PredUpper, nd$UI2, summary)
tapply(nd$PredLower, nd$UI2, summary)
# set factor levels
nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High" ))

nd$GS_1000 <- factor(nd$GS_1000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd <- nd[nd$UI2 %in% c("Urban_Low", "Urban_Moderate","Urban_High"), ]

# plot
p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = GS_1000), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = GS_1000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 3, labeller = as_labeller(c('Urban_Low' = "a Urban_Low", 'Urban_Moderate' = "b Urban_Moderate",'Urban_High' = "c Urban_High"))) + 
  theme_bw() + 
  labs(fill = "% GS", col = "% GS") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(c(-100, 110)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))

print(p1)

### 3000 buffer

AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_3000

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low", "Urban_Moderate", "Urban_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  GS_3000.rs=c(-0.99,0.01493712,1.02,2.03))
# GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd2$GS_3000 <- round(BackTransformCentreredPredictor(
  transformedX = nd2$GS_3000.rs,originalX = predictsSites$GS_3000)*100,0)

nd2$GS_3000[nd2$GS_3000 == 99] <- 100

# 定义区间和对应的标签
breaks <- c(-Inf, 25, 55, 83, 120)  # 划分的区间
labels <- c(25, 50, 75, 100)         # 区间对应的标签

# 使用 cut 函数将 GS_1000 划分为区间，并转换为数字型
nd2$GS_3000 <- cut(nd2$GS_3000, 
                  breaks = breaks, 
                  labels = labels, 
                  right = TRUE,  # 右开区间
                  include.lowest = TRUE)  # 包含下限

# 将切割后的类别转化为数值型（如果需要）
nd2$GS_3000 <- as.numeric(as.character(nd2$GS_3000))

# set values for richness and Abundance
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# set the reference row
refRow <- which((nd2$UI2=="Primary vegetation") & 
                  nd2$StdTmeanAnomaly== nd2$StdTmeanAnomaly[which.min(abs(nd2$StdTmeanAnomaly))]&
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd2$GS_3000==100))
print(refRow)
# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_High"],
  probs = exclQuantiles)
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")
# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd2, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Low") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_Moderate") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Urban_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
print(table(nd2$UI2, !is.na(nd2$PredMedian)))
tapply(nd2$PredMedian, nd2$UI2, summary)
tapply(nd2$PredUpper, nd2$UI2, summary)
tapply(nd2$PredLower, nd2$UI2, summary)
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High" ))

nd2$GS_3000 <- factor(nd2$GS_3000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd2 <- nd2[nd2$UI2 %in% c("Urban_Low", "Urban_Moderate", "Urban_High"), ]

# plot
p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = GS_3000), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = GS_3000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 3, labeller = as_labeller(c('Urban_Low' = "a Urban_Low", 'Urban_Moderate' = "b Urban_Moderate", 'Urban_High' = "c Urban_High"))) + 
  theme_bw() + 
  labs(fill = "% GS", col = "% GS") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(c(-100, 110)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))

print(p2)



### 5000 buffer

AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_5000

nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low","Urban_Moderate","Urban_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  GS_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))
# GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd3$GS_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd3$GS_5000.rs,originalX = predictsSites$GS_5000)*100,0)

nd3$GS_5000[nd3$GS_5000 == 99] <- 100

# 定义区间和对应的标签
breaks <- c(-Inf, 25, 55, 83, 120)  # 划分的区间
labels <- c(25, 50, 75, 100)         # 区间对应的标签

# 使用 cut 函数将 GS_1000 划分为区间，并转换为数字型
nd3$GS_5000 <- cut(nd3$GS_5000, 
                   breaks = breaks, 
                   labels = labels, 
                   right = TRUE,  # 右开区间
                   include.lowest = TRUE)  # 包含下限

# 将切割后的类别转化为数值型（如果需要）
nd3$GS_5000 <- as.numeric(as.character(nd3$GS_5000))

# set values for richness and Abundance
nd3$LogAbund <- 0
nd3$Species_richness <- 0


# set the reference row
refRow <- which((nd3$UI2=="Primary vegetation") & 
                  nd3$StdTmeanAnomaly== nd3$StdTmeanAnomaly[which.min(abs(nd3$StdTmeanAnomaly))] &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd3$GS_5000==100))
print(refRow)
# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Urban_High"],
  probs = exclQuantiles)
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")
# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd3, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Low") & (nd3$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Low") & (nd3$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Moderate") & (nd3$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_Moderate") & (nd3$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_High") & (nd3$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Urban_High") & (nd3$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs
print(table(nd3$UI2, !is.na(nd3$PredMedian)))
tapply(nd3$PredMedian, nd3$UI2, summary)
tapply(nd3$PredUpper, nd3$UI2, summary)
tapply(nd3$PredLower, nd3$UI2, summary)
# set factor levels
nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High" ))

nd3$GS_5000 <- factor(nd3$GS_5000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd3 <- nd3[nd3$UI2 %in% c("Urban_Low","Urban_Moderate", "Urban_High"), ]

# plot
p3 <- ggplot(data = nd3, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = GS_5000), size = 1) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = GS_5000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 3, labeller = as_labeller(
    c('Urban_Low' = "a. Urban Low",
      'Urban_Moderate' = "b. Urban Moderate",
      'Urban_High' = "c. Urban High")),
    strip.position = "top") + 
  theme_bw() + 
  labs(fill = "% GS", col = "% GS") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(c(-100, 110)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))

print(p3)
# 计算不同城市化水平和绿色空间覆盖率下的温度-丰度关系斜率
# 计算不同城市化水平和绿色空间覆盖率下的温度-丰度关系斜率
slope_results <- data.frame()

for(urban in c("Urban_Low", "Urban_Moderate", "Urban_High")) {
  for(gs in c(25, 50, 75, 100)) {
    # 提取并拟合数据
    subset_data <- nd3[nd3$UI2 == urban & nd3$GS_5000 == gs & !is.na(nd3$PredMedian), ]
    if(nrow(subset_data) > 5) {
      model <- lm(PredMedian ~ StdTmeanAnomaly, data = subset_data)
      slope <- coef(model)[2]
      slope_results <- rbind(slope_results, data.frame(
        Urban = urban,
        GS = gs,
        Slope = slope
      ))
    }
  }
}

# 计算每增加25%绿色空间的效应变化
change_results <- data.frame()

for(urban in c("Urban_Low", "Urban_Moderate", "Urban_High")) {
  urban_slopes <- slope_results[slope_results$Urban == urban, ]
  
  if(nrow(urban_slopes) > 1) {
    for(i in 1:(nrow(urban_slopes)-1)) {
      if(urban_slopes$GS[i+1] - urban_slopes$GS[i] == 25) {
        # 计算斜率变化及百分比
        slope_change <- urban_slopes$Slope[i+1] - urban_slopes$Slope[i]
        if(urban_slopes$Slope[i] != 0) {
          percent_change <- (slope_change / abs(urban_slopes$Slope[i])) * 100
        } else {
          percent_change <- NA
        }
        
        change_results <- rbind(change_results, data.frame(
          Urban = urban,
          GS_from = urban_slopes$GS[i],
          GS_to = urban_slopes$GS[i+1],
          Percent_change = percent_change
        ))
      }
    }
  }
}

# 计算平均效应变化及标准误
avg_effects <- aggregate(
  Percent_change ~ Urban, 
  data = change_results, 
  FUN = function(x) {
    c(mean = mean(x, na.rm = TRUE),
      se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  }
)

# 输出结果
cat("\n绿色空间覆盖率每增加25%对昆虫丰度温度敏感性的影响:\n")
for(i in 1:nrow(avg_effects)) {
  cat(sprintf("%s: %.1f%% (±%.1f SE)\n", 
              avg_effects$Urban[i],
              avg_effects$Percent_change[i, "mean"], 
              avg_effects$Percent_change[i, "se"]))
}

# 输出每个具体变化区间的效应
cat("\n每个25%绿色空间变化区间的具体影响:\n")
for(urban in c("Urban_Low", "Urban_Moderate", "Urban_High")) {
  urban_changes <- change_results[change_results$Urban == urban, ]
  if(nrow(urban_changes) > 0) {
    cat("\n", urban, ":\n")
    for(i in 1:nrow(urban_changes)) {
      cat(sprintf("  %d%% to %d%%: %.1f%% 变化\n", 
                  urban_changes$GS_from[i],
                  urban_changes$GS_to[i],
                  urban_changes$Percent_change[i]))
    }
  }
}

### 5000m buffer

RichMeanAnomalyModel1 <- RichMeanAnomalyModel_5000

### mean anomaly species richness 
nd31 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Urban_Low","Urban_Moderate","Urban_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  GS_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))
#GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd31$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd31$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd31$GS_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd31$GS_5000.rs,originalX = predictsSites$GS_5000)*100,0)

nd31$GS_5000[nd31$GS_5000 == 99] <- 100

# 定义区间和对应的标签
breaks <- c(-Inf, 25, 55, 83, 120)  # 划分的区间
labels <- c(25, 50, 75, 100)         # 区间对应的标签

# 使用 cut 函数将 GS_1000 划分为区间，并转换为数字型
nd31$GS_5000 <- cut(nd31$GS_5000, 
                   breaks = breaks, 
                   labels = labels, 
                   right = TRUE,  # 右开区间
                   include.lowest = TRUE)  # 包含下限

# 将切割后的类别转化为数值型（如果需要）
nd31$GS_5000 <- as.numeric(as.character(nd31$GS_5000))
# set values for richness and Abundance
nd31$LogAbund <- 0
nd31$Species_richness <- 0

# set the reference row
refRow <- which((nd31$UI2=="Primary vegetation") & 
                  nd31$StdTmeanAnomaly== nd31$StdTmeanAnomaly[which.min(abs(nd31$StdTmeanAnomaly))] &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd31$GS_5000==100))
print(refRow)
# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Urban_Low"],
  probs = exclQuantiles)
QAL <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Urban_Moderate"],
  probs = exclQuantiles)
QAH <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Urban_High"],
  probs = exclQuantiles)
cat("Primary vegetation range:", QPV[1], "to", QPV[2], "")
cat("Urban_Low range:", QSV[1], "to", QSV[2], "")
cat("Urban_Moderate range:", QAL[1], "to", QAL[2], "")
cat("Urban_High range:", QAH[1], "to", QAH[2], "")
# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd31)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd31$UI2=="Primary vegetation") & (nd31$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd31$UI2=="Primary vegetation") & (nd31$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd31$UI2=="Urban_Low") & (nd31$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd31$UI2=="Urban_Low") & (nd31$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd31$UI2=="Urban_Moderate") & (nd31$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd31$UI2=="Urban_Moderate") & (nd31$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd31$UI2=="Urban_High") & (nd31$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd31$UI2=="Urban_High") & (nd31$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd31$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd31$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd31$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs
print(table(nd31$UI2, !is.na(nd31$PredMedian)))
tapply(nd31$PredMedian, nd31$UI2, summary)
tapply(nd31$PredUpper, nd31$UI2, summary)
tapply(nd31$PredLower, nd31$UI2, summary)
# set factor levels
nd31$UI2 <- factor(nd31$UI2, levels = c("Primary vegetation", "Urban_Low", "Urban_Moderate", "Urban_High" ))

nd31$GS_5000 <- factor(nd31$GS_5000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd31 <- nd31[nd31$UI2 %in% c("Urban_Low","Urban_Moderate", "Urban_High"), ]

# plot
p4 <- ggplot(data = nd31, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = GS_5000), size = 1) +
  geom_ribbon(aes(ymin = nd31$PredLower, ymax = nd31$PredUpper, fill = GS_5000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 3, labeller = as_labeller(
    c('Urban_Low' = "d. Urban Low",
      'Urban_Moderate' = "e. Urban Moderate",
      'Urban_High' = "f. Urban High")),
    strip.position = "top") +  # 将标签移动到顶部
  theme_bw() +
  labs(fill = "% GS", col = "% GS") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(c(-100, 110)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))
print(p4)
# 使用 gridExtra 合并
grid.arrange(p3, p4, nrow = 2) 
### 10000 buffer


AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_10000

nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  GS_10000.rs=c(-1.08,-0.01,1.045653,2.08))
# GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd4$GS_10000 <- round(BackTransformCentreredPredictor(
  transformedX = nd4$GS_10000.rs,originalX = predictsSites$GS_10000)*100,0)

nd4$GS_10000[nd4$GS_10000 == 99] <- 100

# set values for richness and Abundance
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# set the reference row
refRow <- which((nd4$UI2=="Primary vegetation") & 
                  nd4$StdTmeanAnomaly== min(nd4$StdTmeanAnomaly[nd4$StdTmeanAnomaly > 0]) &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd4$GS_10000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd4, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd4$GS_10000 <- factor(nd4$GS_10000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd4 <- nd4[nd4$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = GS_10000), size = 1) +
  geom_ribbon(aes(ymin = nd4$PredLower, ymax = nd4$PredUpper, fill = GS_10000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "g              Agriculture_Low", 'Agriculture_High' = "h              Agriculture_High"))) + 
  theme_bw() + 
  labs(fill = "% GS", col = "% GS") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 1)) +
  ylim(c(-100, 110)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))



cowplot::plot_grid(p1, p2, p3, p4, nrow = 4, labels = c("1000m", "3000m", "5000m", "10000m"))

# save
ggsave(filename = paste0(outDir, "Abun_buffers_plots.pdf"), height = 16, width = 9)



#### Plots for Richness results ####

### 1000m buffer

RichMeanAnomalyModel1 <- RichMeanAnomalyModel_1000

### mean anomaly species richness 
nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  GS_1000.rs=c(-0.96,0.01493712,1.00,2.0))
#GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd$GS_1000 <- round(BackTransformCentreredPredictor(
  transformedX = nd$GS_1000.rs,originalX = predictsSites$GS_1000)*100,0)

nd$GS_1000[nd$GS_1000 == 99] <- 100

# set values for richness and Abundance
nd$LogAbund <- 0
nd$Species_richness <- 0

# set the reference row
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) &
                  (nd$GS_1000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd$GS_1000 <- factor(nd$GS_1000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd <- nd[nd$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = GS_1000), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = GS_1000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "a              Agriculture_Low", 'Agriculture_High' = "b              Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% GS", col = "% GS") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))



### 3000m buffer

RichMeanAnomalyModel1 <- RichMeanAnomalyModel_3000

### mean anomaly species richness 
nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  GS_3000.rs=c(-0.99,0.01493712,1.02,2.03))
#GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd2$GS_3000 <- round(BackTransformCentreredPredictor(
  transformedX = nd2$GS_3000.rs,originalX = predictsSites$GS_3000)*100,0)

nd2$GS_3000[nd2$GS_3000 == 99] <- 100

# set values for richness and Abundance
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# set the reference row
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) &
                  (nd2$GS_3000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd2)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd2$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd2$GS_3000 <- factor(nd2$GS_3000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd2 <- nd2[nd2$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = GS_3000), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = GS_3000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "c              Agriculture_Low", 'Agriculture_High' = "d              Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% GS", col = "% GS") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))





### 10000m buffer


RichMeanAnomalyModel1 <- RichMeanAnomalyModel_10000

### mean anomaly species richness 
nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  GS_10000.rs=c(-1.08,-0.01,1.045653,2.08))
#GS_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform GS data range
nd4$GS_10000 <- round(BackTransformCentreredPredictor(
  transformedX = nd4$GS_10000.rs,originalX = predictsSites$GS_10000)*100,0)

nd4$GS_10000[nd4$GS_10000 == 99] <- 100

# set values for richness and Abundance
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# set the reference row
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomaly==min(abs(nd4$StdTmeanAnomaly))) &
                  (nd4$GS_10000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd4)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd4$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd4$GS_10000 <- factor(nd4$GS_10000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd4 <- nd4[nd4$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = GS_10000), size = 1) +
  geom_ribbon(aes(ymin = nd4$PredLower, ymax = nd4$PredUpper, fill = GS_10000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "g              Agriculture_Low", 'Agriculture_High' = "h              Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% GS", col = "% GS") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))




cowplot::plot_grid(p1, p2, p3, p4, nrow = 4, labels = c("1000m", "3000m", "5000m", "10000m"))

# save
ggsave(filename = paste0(outDir, "Rich_buffers_plots.pdf"), height = 16, width = 9)


t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()

