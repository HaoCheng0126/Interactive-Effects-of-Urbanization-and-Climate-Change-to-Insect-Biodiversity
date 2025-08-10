##%######################################################%##  
#  
#  
####  
# 气候数据组织与绘图  
####  
#  
#  
##%######################################################%##  
# 本脚本用于探索气候数据，生成温度异常值及相关地图  
# 本脚本也可用于使用不同温度阈值生成异常值数据集  
# 论文中使用的阈值是10摄氏度  

# 设置目录  
dataDir <- "/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/0_data/"  
outDir <- "/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/3_PrepareClimateIndexMaps/"  
#outDir <- "10_SCA_Baseline_testing/"  
#outDir <- "14_Additional_Tests/"  
if(!dir.exists(outDir)) dir.create(outDir)  
# sink(paste0(outDir,"log.txt"))  
t.start <- Sys.time()  # 记录开始时间  
print(t.start)  

# 加载所需的R包  
library(raster)      # 用于栅格数据处理  
library(sp)          # 空间数据处理  
library(dismo)       # 物种分布建模  
library(sf)
library(rgdal)       # 地理空间数据抽象库  
library(RColorBrewer) # 颜色配置  
library(ncdf4)       # 处理NetCDF文件  
library(rasterVis)   # 栅格数据可视化  
library(gridExtra)   # 组合图表  
library(cowplot)     # 组合多个ggplot图  
library(ggplot2)     # 绘图系统  
library(viridis)     # 颜色映射  
library(snow)        # 并行计算  

# 从CRU加载平均温度数据  
tmp <- stack(paste0(dataDir,"cru_ts4.09.1901.2024.tmp.dat.nc"),varname="tmp")  

# 获取1901至1930年的值作为基线  
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]  

# 用于测试不同基线长度：  
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]  
tmp1901_1910 <- tmp[[names(tmp)[1:120]]]  
tmp1901_1920 <- tmp[[names(tmp)[1:240]]]  

# 当代温度范围(2005年中期)  
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]  

# 更近期的数据，用于Extended Data Figure 7  
#tmp2016_18<- tmp[[names(tmp)[1381:1416]]]  

# 设置昆虫活跃月份的温度阈值  
thresh <- 10 # 可选 6, 8, 10 摄氏度  

# 确定使用哪个栅格作为现在  
pre_ras <- tmp2004_6  
#pre_ras <- tmp2016_18  

# 获取栅格值中非NA位置的向量  
ras <- pre_ras[[1]]  
vals <- values(ras)  

# 基于非NA单元格创建点集  
wgs84 <- crs(tmp)  
pnts <- rasterToPoints(ras, spatial = T)  
SP <- SpatialPoints(pnts, proj4string=wgs84)  

# 设置并行处理环境  
nCores <- parallel::detectCores()  
st1 <- Sys.time()  
cl <- snow::makeCluster(nCores-1)  

# 导出变量到并行集群  
snow::clusterExport(  
  cl = cl,  
  list = c('pre_ras', 'values', 'names', 'length', 'mean', 'sd',  
           'tmp', 'SP','rasterize','crop','trim', 'grep', 'sapply', 'strsplit',  
           'cellStats', 'thresh', 'tmp1901_1930', 'tmp1901_1905', 'tmp1901_1910', 'tmp1901_1920'),envir = environment())  

# 使用并行计算处理温度变量  
temperatureVars <- data.frame(t(parSapply(  
  cl = cl,X = (1:length(SP)),FUN = function(i){  
    # 对每个空间点进行处理  
    temperatureVars <-NULL  
    
    # 使用掩码提高处理速度  
    mask <- trim(rasterize(SP[i, ], pre_ras[[1]]))  
    mapCrop <- crop(pre_ras, mask)  
    if(!length(names(mapCrop)[values(mapCrop) >= thresh]) == 0 &  
       length(values(mapCrop)[!is.na(values(mapCrop))]) > 0 ){  
      # 首先识别该单元格的昆虫活跃月份  
      # 获取5年中每个月的平均温度  
      vals <- NULL  
      # 对于每个月，获取5年的平均温度  
      for(j in 1:12){  
        if(j < 10){ mon <- paste0(0, j) }else {mon <- j}  
        monthmean <- values(mean(mapCrop[[grep(mon, sapply(strsplit(names(mapCrop), "[.]"), "[[", 2))  
        ]]))  
        vals <- rbind(vals, c(mon, monthmean))  
      }  
      vals <- as.data.frame(vals)  
      vals$V2 <- as.numeric(as.character(vals$V2))  
      # 哪些月份的5年平均值大于等于阈值  
      vals <- vals[vals$V2 >= thresh, ]  
      
      # 有时平均值不超过阈值，即使个别月份超过  
      if(nrow(vals) == 0){  
        avg_temp = NA  
        n_months = NA  
        Anom <- NA  
        StdAnom <- NA  
        return(c(avg_temp = avg_temp, n_months = n_months, Anom = Anom, StdAnom = StdAnom))  
      }else{  
        # 哪些是合适的月份  
        months <- vals$V1  
        # 有多少月份达到或超过阈值？  
        n_months <- length(months)  
        # 计算"当代"平均值和标准差  
        avg_temp <- mean(vals$V2)  
        
        ### 现在计算活跃月份的基线平均值和标准差 ###  
        # 获取该网格单元所有年份的值  
        baseline <- crop(tmp1901_1930, mask)  
        #baseline <- crop(tmp1901_1905, mask)  
        #baseline <- crop(tmp1901_1910, mask)  
        #baseline <- crop(tmp1901_1920, mask)  
        
        # 将基线仅限于所需的月份  
        baseline <-  
          baseline[[names(baseline)[sapply(strsplit(names(baseline), "[.]"), "[[", 2) %in% months]]]  
        
        # 获取平均值和标准差  
        mean_baseline <- mean(values(baseline))  
        sd_mean_baseline <- sd(values(baseline))  
        
        # 现在计算异常值和标准化异常值  
        Anom <- avg_temp - mean_baseline  
        StdAnom <-  
          Anom/sd_mean_baseline  
        
        return(c(avg_temp = avg_temp, n_months = n_months, Anom = Anom, StdAnom = StdAnom))  
      }}else{ # 0/NA检查后  
        avg_temp = NA  
        n_months = NA  
        Anom <- NA  
        StdAnom <- NA  
        return(c(avg_temp = avg_temp, n_months = n_months, Anom = Anom, StdAnom = StdAnom))  
      }  
  } # 函数结束  
)))  

# 停止并行集群  
snow::stopCluster(cl)  
st2 <- Sys.time()  
print(st2 - st1) # 时间差5.196104小时  

# 保存结果  
save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2004_06_thresh_", thresh, ".rdata"))  

#### 查看不同指标之间的相关性 ####  
temperatureVars <- as.data.frame(temperatureVars) # 67420行  
# 移除NA值  
temp_data <- temperatureVars[!is.na(temperatureVars$avg_temp), ] # 58319行  

# 计算平均温度与异常值的相关性  
cor(temp_data$avg_temp, temp_data$Anom) # -0.2040909, 2004-6版本, 阈值10  


# Enhanced correlation scatter plot
ggplot(data = temp_data, aes(x = avg_temp, y = Anom)) + 
  # Add semi-transparent points with improved color
  geom_point(size = 0.8, alpha = 0.6, color = "#1F77B4") + 
  # Fix the warning by using linewidth instead of size
  geom_smooth(method = "lm", linewidth = 1.5, color = "#D62728", fill = "#D62728", alpha = 0.15) +
  # Use a cleaner theme with customizations
  theme_minimal() + 
  # Add grid lines only for major breaks
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 12, color = "#555555", hjust = 0),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Improve labels with better formatting
  labs(
    x = "Average temperature of the location (2004-2006, active months, °C)",
    y = "Temperature anomaly (difference between
present and baseline temperatures, °C)",
    title = "Global Temperature Anomaly vs. Average Temperature",
    subtitle = paste0("Pearson's correlation coefficient = ", 
                      round(cor(temp_data$avg_temp, temp_data$Anom, use = "complete.obs"), digits = 2))
  )

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_Anom.pdf"))


# 移除少数异常值，一些无穷大值  
nrow(temp_data[temp_data$StdAnom >12, ]) # 14行  
temp_data <- temp_data[temp_data$StdAnom < 12, ] # 57998行  

# 计算平均温度与标准化异常值的相关性  
cor(temp_data$avg_temp, temp_data$StdAnom) # -0.155, 2004-06, 阈值10  

# 绘制相关性散点图  
# Enhanced standardized anomaly scatter plot with consistent styling
ggplot(data = temp_data, aes(x = avg_temp, y = StdAnom)) + 
  # Add semi-transparent points with improved color
  geom_point(size = 0.8, alpha = 0.6, color = "#1F77B4") + 
  # Fix the warning by using linewidth instead of size
  geom_smooth(method = "lm", linewidth = 1.5, color = "#D62728", fill = "#D62728", alpha = 0.15) +
  # Use a cleaner theme with customizations
  theme_minimal() + 
  # Add grid lines only for major breaks
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 12, color = "#555555", hjust = 0),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Improve labels with better formatting
  labs(
    x = "Average temperature of the location (2004-2006, active months, °C)",
    y = "Standardised climate anomaly",
    title = "Global Temperature Standardised Anomaly vs. Average Temperature",
    subtitle = paste0("Pearson's correlation coefficient = ", 
                      round(cor(temp_data$avg_temp, temp_data$StdAnom, use = "complete.obs"), digits = 2))
  )
ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_StdAnom.pdf"))

# 移除大于3的标准化异常值  
# Enhanced standardized anomaly scatter plot with outliers removed
ggplot(data = temp_data[temp_data$StdAnom <= 3, ], aes(x = avg_temp, y = StdAnom)) + 
  # Add semi-transparent points with improved color
  geom_point(size = 0.8, alpha = 0.6, color = "#1F77B4") + 
  # Fix the warning by using linewidth instead of size
  geom_smooth(method = "lm", linewidth = 1.5, color = "#D62728", fill = "#D62728", alpha = 0.15) +
  # Use a cleaner theme with customizations
  theme_minimal() + 
  # Add grid lines only for major breaks
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 12, color = "#555555", hjust = 0),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Improve labels with better formatting in English
  labs(
    x = "Average temperature of the location (2004-2006, active months, °C)",
    y = "Standardised climate anomaly",
    title = "Global Temperature Standardised Anomaly vs. Average Temperature",
    subtitle = paste0("Outliers removed (StdAnom ≤ 3), Pearson's correlation coefficient = ", 
                      round(cor(temp_data[temp_data$StdAnom <= 3, 'avg_temp'], 
                                temp_data[temp_data$StdAnom <= 3, 'StdAnom'], 
                                use = "complete.obs"), digits = 2))
  )

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_StdAnom_outliersrem.pdf"))  

# 计算异常值与标准化异常值的相关性  
cor(temp_data$Anom, temp_data$StdAnom) # 0.3247709, 2004-06, 阈值10  

# 现在查看不同区域  
# 添加经纬度信息  
SP_df <- as.data.frame(SP)  
SP_df <- cbind(SP_df, temperatureVars)  

# 按区域分类  
SP_df$Tropical <- NA  
SP_df[SP_df$y > -23.44 & SP_df$y < 23.44, 'Tropical'] <- "Tropical"  # 热带区域  
SP_df[is.na(SP_df$Tropical), 'Tropical'] <- "Temperate"              # 温带区域  
SP_df <- SP_df[!is.na(SP_df$avg_temp), ]  

table(SP_df$Tropical)  
# Temperate  Tropical 
# 39745     18267 

# 计算不同区域的相关性  
cor_temp <- round(cor(SP_df[SP_df$Tropical == "Temperate", "avg_temp"], SP_df[SP_df$Tropical == "Temperate", "Anom"]), digits = 2) # 0.02  
cor_trop <- round(cor(SP_df[SP_df$Tropical == "Tropical", "avg_temp"], SP_df[SP_df$Tropical == "Tropical", "Anom"]), digits = 2) # 0.14  

# 在标签中添加相关系数  
SP_df$Tropical <- sub("Temperate", paste0("Temperate, cor = ", cor_temp), SP_df$Tropical)  
SP_df$Tropical <- sub("Tropical", paste0("Tropical, cor = ", cor_trop), SP_df$Tropical)  

# 按区域绘制相关性散点图  
# Enhanced faceted scatter plot by tropical regions
ggplot(data = SP_df, aes(x = avg_temp, y = Anom)) + 
  # Add semi-transparent points with improved color
  geom_point(size = 0.8, alpha = 0.6, color = "#1F77B4") + 
  # Fix the warning by using linewidth instead of size
  geom_smooth(method = "lm", linewidth = 1.5, color = "#D62728", fill = "#D62728", alpha = 0.15) +
  # Maintain the faceting by tropical regions
  facet_wrap(~ Tropical) +
  # Use a cleaner theme with customizations
  theme_minimal() + 
  # Add grid lines only for major breaks
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Improve labels with better formatting
  labs(
    x = "Average temperature of the location (2016-2018, active months, °C)",
    y = "Temperature anomaly (difference between
present and baseline temperatures, °C)"
  )

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_Anom_REALM.pdf"))  

# 重置区域分类  
SP_df$Tropical <- NA  
SP_df[SP_df$y > -23.44 & SP_df$y < 23.44, 'Tropical'] <- "Tropical"  
SP_df[is.na(SP_df$Tropical), 'Tropical'] <- "Temperate"  

# 移除异常值  
SP_df <- SP_df[SP_df$StdAnom < 12, ] # 58313行  

# 计算各区域标准化异常值的相关性  
cor_temp <- round(cor(SP_df[SP_df$Tropical == "Temperate", "avg_temp"], SP_df[SP_df$Tropical == "Temperate", "StdAnom"]), digits = 2) # -0.55  
cor_trop <- round(cor(SP_df[SP_df$Tropical == "Tropical", "avg_temp"], SP_df[SP_df$Tropical == "Tropical", "StdAnom"]), digits = 2) # 0.08  

# 在标签中添加相关系数  
SP_df$Tropical2 <- sub("Temperate", paste0("Temperate, cor = ", cor_temp), SP_df$Tropical)  
SP_df$Tropical2 <- sub("Tropical", paste0("Tropical, cor = ", cor_trop), SP_df$Tropical)  

# 按区域绘制标准化异常值相关性散点图  
# Enhanced faceted scatter plot of standardized anomalies by tropical regions
ggplot(data = SP_df, aes(x = avg_temp, y = StdAnom)) + 
  # Add semi-transparent points with improved color
  geom_point(size = 0.8, alpha = 0.6, color = "#1F77B4") + 
  # Fix the warning by using linewidth instead of size
  geom_smooth(method = "lm", linewidth = 1.5, color = "#D62728", fill = "#D62728", alpha = 0.15) +
  # Maintain the faceting by tropical regions with free scales
  facet_wrap(~ Tropical2, scales = "free") +
  # Use a cleaner theme with customizations
  theme_minimal() + 
  # Add grid lines only for major breaks and better facet styling
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Improve labels with better formatting
  labs(
    x = "Average temperature of the location (2016-2018, active months, °C)",
    y = "Standardised climate anomaly"
  )

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_StdAnom_REALM_1618.pdf"))  

cor_temp <- round(cor(SP_df[SP_df$Tropical == "Temperate", "Anom"], SP_df[SP_df$Tropical == "Temperate", "StdAnom"]), digits = 2) # 0.37, 2004-6  

##%######################################################%##  
#  
#  
####  
# 论文图表  
####  
#  
#  
##%######################################################%##  

# 将数据转换为数据框  
temperatureVars2 <- as.data.frame(temperatureVars)  

# 将数据添加到点的经纬度  
SP_df <- as.data.frame(SP)  
SP_df <- cbind(SP_df, temperatureVars2)  

# 快速查看全球数据图  
avg_temp_ras <- rasterFromXYZ(SP_df[ , 1:3])  
Anom_ras <- rasterFromXYZ(SP_df[ , c(1,2,5)])  
StdAnom_ras <- rasterFromXYZ(SP_df[ , c(1,2,6)])  
n_months <- rasterFromXYZ(SP_df[ , c(1,2,4)])  
plot(avg_temp_ras)  
plot(Anom_ras)  
plot(StdAnom_ras)  
plot(n_months)  

### 首先，绝对变化 ###  
# 将栅格转换为数据框  
plot_data <- SP_df[, c(1,2,5)]  

# 组织断点、颜色和标签  
brks <- c(-0.6,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,3,5)  
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:8],  
          (brewer.pal(n = 8,name = "Purples"))[4:6],  
          (brewer.pal(n = 8,name = "Oranges"))[3:5])  
labs <- c("-0.6 : -0.2","-0.2 : -0.1","-0.1 : 0",  
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 3", ">3")  

# 将值分配到分类区间  
plot_data$bins <- cut(plot_data$Anom,  
                      breaks = brks,  
                      labels = labs,  
                      include.lowest = TRUE)  

# 获取用于轮廓的世界地图  
world_map <- map_data("world")  

# 创建绝对温度变化地图  
# Enhanced absolute temperature change map
p1 <- ggplot(plot_data[!is.na(plot_data$Anom),]) + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "lightgrey", fill = "white",  size = 0.1) +
  geom_tile(aes(x = x, y = y, fill = bins), na.rm = TRUE) +
  scale_fill_manual(values = cols) + 
  xlab("") +
  ylab("") +
  labs(fill = "Absolute\nTemperature\nChange") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.3), 
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        #legend.key.width = unit(3, "cm"),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7), 
        legend.key.size = unit(0.2,"cm"), legend.direction = "vertical",
        title = element_text(size = 8, face = "bold")) +
  guides(fill = guide_legend(title.position = "top") ) + 
  ggtitle("a")

print(p1)
# 计算每行的平均气候值  
rows <- init(Anom_ras, v='row')  
ravg <- zonal(Anom_ras, rows, fun = 'mean', na.rm = T)  
ravg[is.nan(ravg)] <- NA  
ravg <- as.data.frame(ravg)  

# 绘制边缘分布图  
# Improved marginal distribution plot
p2 <- ggplot(data = ravg) +
  # Use proper aesthetic mapping without direct $ references
  geom_line(aes(x = zone, y = mean), col = "#8379BD", linewidth = 0.8) +
  geom_ribbon(aes(x = zone, ymin = min(mean, na.rm = TRUE), ymax = mean), 
              fill = "#473C8B", alpha = 0.7) +
  # Set scales properly without direct $ references
  scale_x_reverse(limits = c(300, 1), expand = c(0, 0)) +
  # Calculate limits using summarize functions instead of $ notation
  scale_y_continuous(
    limits = function(x) {
      c(min(ravg$mean, na.rm = TRUE), max(ravg$mean, na.rm = TRUE))
    }, 
    expand = c(0, 0)
  ) +
  # Use minimal theme base
  theme_minimal() +
  # Custom theme elements
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    # Add clean white background for consistency
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Maintain flipped coordinates
  coord_flip()

print(p2)
#### 现在标准化异常值 ####  
# 将栅格转换为数据框  
plot_data2 <- SP_df[, c(1,2,6)]  

# 组织断点、颜色和标签  
brks2 <- c(-0.65,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,3, 5)  
cols2 <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:8],  
           (brewer.pal(n = 8,name = "Purples"))[4:6],  
           (brewer.pal(n = 8,name = "Oranges"))[3:5])  
labs2 <- c("-0.6 : -0.2","-0.2 : -0.1","-0.1 : 0",  
           "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 3", "> 3")  

# 将值分配到分类区间  
plot_data2$bins <- cut(plot_data2$StdAnom,  
                       breaks = brks2,  
                       labels = labs2,  
                       include.lowest = TRUE)  
plot_data2 <- plot_data2[!is.na(plot_data2$bins), ]  

# 创建标准化温度异常地图  
p3 <- ggplot(plot_data2[!is.na(plot_data2$StdAnom),]) +  
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "lightgrey", fill = "white", size = 0.1) +  
  geom_raster(aes(x = x, y = y, fill = bins), na.rm = TRUE) +  
  scale_fill_manual(values = cols2) +  
  xlab("") +  
  ylab("") +  
  labs(fill = "Standardised\nTemperature\nAnomaly") +  
  theme_bw() +  
  theme(legend.position = c(0.15, 0.3),  
        panel.border = element_blank(),  
        panel.grid = element_blank(),  
        axis.text = element_blank(),  
        axis.ticks = element_blank(),  
        legend.text = element_text(size = 6),  
        legend.title = element_text(size = 7),  
        legend.key.size = unit(0.2,"cm"), legend.direction = "vertical",  
        title = element_text(size = 8, face = "bold")) +  
  guides(fill = guide_legend(title.position = "top") ) +  
  ggtitle("b")  

print(p3)
# 计算每行的平均气候值  
rows2 <- init(StdAnom_ras, v='row')  
ravg2 <- zonal(StdAnom_ras, rows2, fun = 'mean', na.rm = T)  
ravg2[is.nan(ravg2)] <- NA  
ravg2 <- as.data.frame(ravg2)  

# 移除无穷大值  
ravg2 <- ravg2[!ravg2$mean == Inf, ]  

# 绘制边缘分布图  
p4 <- ggplot(data = ravg2) +  
  geom_line(aes(x = zone, y = mean), col = c("#8379BD")) +  
  geom_ribbon(aes(ymin = min(ravg2$mean, na.rm = T), ymax = mean, x = zone), fill = c("#473C8B"), alpha = 0.7) +  
  theme_bw() +  
  scale_x_reverse(limits = c(300, 1), expand = c(0,0)) +  
  scale_y_continuous(limits = c(min(ravg2$mean, na.rm = T), max(ravg2$mean, na.rm = T)), expand = c(0,0)) +  
  theme(panel.border = element_blank(),  
        panel.grid = element_blank(),  
        axis.text = element_blank(),  
        axis.title = element_blank(),  
        axis.ticks = element_blank()  
  ) +  
  coord_flip()  

# 将图和图例组织成一个对象  
final_plot <- cowplot::plot_grid(  
  cowplot::plot_grid(  
    p1, p2, nrow = 1, align = "hv", rel_widths = c(3,1)),  
  cowplot::plot_grid(  
    p3, p4, nrow = 1, align = "hv", rel_widths = c(3,1)),  
  nrow = 2  
)  
print(final_plot)
# 保存为PDF  
ggsave(filename = paste0(outDir, "Extended_Data1_maps_thresh_", thresh, ".pdf"), plot = last_plot(), width = 183, height = 200, units = "mm", dpi = 300)  
ggsave(filename = paste0(outDir, "Extended_Data1_maps_thresh_", thresh, ".jpeg"), plot = last_plot(), width = 183, height = 200, units = "mm", dpi = 300)  

#### 图 - 活跃月份数量地图 ####  
map_data <- SP_df[, c(1,2,4)]  

# 创建活跃月份地图  
ggplot(map_data[!is.na(map_data$n_months),]) +  
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "lightgrey", fill = "white", size = 0.1) +  
  geom_raster(aes(x = x, y = y, fill = n_months), na.rm = TRUE) +  
  scale_fill_gradient(low = c("#79CDCD"), high = c("#00688B"), na.value = NA, limits = c(1, 12)) +  
  xlab("") +  
  ylab("") +  
  labs(fill = "Number of months\nabove 10 degrees C") +  
  theme_bw() +  
  theme(legend.position = 'bottom',  
        panel.border = element_blank(),  
        panel.grid = element_blank(),  
        axis.text = element_blank(),  
        axis.ticks = element_blank(),  
        legend.text = element_text(size = 8),  
        legend.title = element_text(size = 10)) +  
  guides(colour = guide_colourbar(show.limits = TRUE))  
ggsave(filename = paste0(outDir, "Nmonths_plot_thresh_", thresh, ".pdf"), width = 4, height = 3)  

# 记录结束时间  
t.end <- Sys.time()  
print(round(t.end - t.start,0))  
sink()  
```