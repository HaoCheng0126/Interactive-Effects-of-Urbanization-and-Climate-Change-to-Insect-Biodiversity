##%######################################################%##
#                                                          #
####         Organise PREDICTS data for insects         ####
#                                                          #
##%######################################################%##

# This script takes the complete PREDICTS database, selects those entries for
# insects, and organises the data for analysis.

# directories
dataDir <- "/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/0_data/"
outDir <- "/Users/chenghao/Desktop/Final Project-R CODE/CH_Urbanization&CliamteChange/1_PreparePREDICTSData/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load required libraries
library(predictsFunctions)
library(ggplot2)

sessionInfo()

# Set the path to your local copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)

# Para contar los valores únicos en la columna "study"
num_unique_studies <- length(unique(predicts$Study_name))

# Para mostrar el resultado
print(paste("Hay", num_unique_studies, "estudios únicos en la base de datos PREDICTS"))

# Select only data for insects
predicts <- predicts[(predicts$Class=="Insecta"),]

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
# Correcting 870378 values for sensitivity to sampling effort (most of these are 0s, 19184 non-zero)

table(predicts$Diversity_metric)


# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those sites that are the problem
# problem reported to NHM PREDICTS team.
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]


# Merge sites that have the same coordinates (e.g. multiple traps on a single transect)
predicts <- MergeSites(diversity = predicts)


# remove rows where land use or use intensity info is missing
predicts.complete <- droplevels(predicts[(
  predicts$Predominant_land_use!="Cannot decide"),])
predicts.complete <- droplevels(predicts.complete[(
  predicts.complete$Use_intensity!="Cannot decide"),])

nrow(predicts.complete)
# 851631 records

# get counts of n species in major groups
species <- unique(predicts.complete[,c('Order','Taxon_name_entered')])

order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))


# Calculate site metrics of diversity
sites <- SiteMetrics(diversity = predicts,
                     extra.cols = c("Predominant_land_use",
                                    "SSB","SSBS", "Biome"))
# Computing site metrics for 899323 measurements
# The data contain 244 sources, 343 studies and 8685 sites
# Computing site-level values
# Computing total abundance
# Computing species richness
# Computing Simpson's diversity
# Computing Chao
# Computing Rarefied Species Richness

# First, we will rearrange the land-use classification a bit
sites$LandUse <- paste(sites$Predominant_land_use)

# Drop classification where land use could not be identified
sites$LandUse[(sites$LandUse=="Cannot decide")] <- NA

# Now make the variable a factor, and set the reference level to primary vegetation
sites$LandUse <- factor(sites$LandUse)
sites$LandUse <- relevel(sites$LandUse,ref="Primary vegetation")

sites$Use_intensity[sites$Use_intensity=="Cannot decide"] <- NA

# combine LU and UI 
sites$UI <- paste0(sites$LandUse,'_',sites$Use_intensity)
sites$UI[grep("NA",sites$UI)] <- NA

# recode according to land use and use intensity combinations
sites$UI2 <- dplyr::recode(sites$UI,
                           'Primary vegetation_Minimal use' = 'Primary vegetation',
                           'Cropland_Light use' = 'Agriculture',
                           'Secondary vegetation (indeterminate age)_Minimal use' = 'Secondary vegetation',
                           'Urban_Light use' = 'Urban_Moderate',
                           'Secondary vegetation (indeterminate age)_Light use' = 'Secondary vegetation',
                           'Cropland_Intense use' = 'Agriculture',
                           'Cropland_Minimal use' = 'Agriculture',
                           'Pasture_Light use' = 'Agriculture',
                           'Pasture_Minimal use' = 'Agriculture',
                           'Intermediate secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Secondary vegetation (indeterminate age)_Intense use' = 'Secondary vegetation',
                           'Pasture_Intense use' = 'Agriculture',
                           'Urban_Minimal use' = 'Urban_Low',
                           'Primary vegetation_Light use' = 'Primary vegetation',
                           'Young secondary vegetation_Light use' = 'Secondary vegetation',
                           'Urban_Intense use' = 'Urban_High',
                           'Primary vegetation_Intense use' = 'Primary vegetation',
                           'Young secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Minimal use' = 'Agriculture',
                           'Plantation forest_Intense use' = 'Agriculture',
                           'Young secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Light use' = 'Agriculture',
                           'Mature secondary vegetation_Light use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Light use' = 'Secondary vegetation')

# 
sites$Use_intensity[((sites$LandUse=="Mature secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Intermediate secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Young secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"

# remove the Agriculture sites and sites that are NA in UI2
sites <- sites[!sites$UI2 == "Agriculture", ]
sites <- sites[!sites$UI2 == "Secondary vegetation", ]
sites <- sites[!is.na(sites$UI2), ]


sites <- droplevels(sites)

# transform abundance values 
sites$LogAbund <- log(sites$Total_abundance+1)


# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]


# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds"))



##%######################################################%##
#                                                          #
#### redo dataset summaries after removing other sites  ####
#                                                          #
##%######################################################%##


predicts2 <- predicts[predicts$SSBS %in% sites$SSBS, ]


table(predicts2$Diversity_metric)


nrow(predicts2)
# 300038 records

# get counts of n species in major groups
species <- unique(predicts2[,c('Order','Taxon_name_entered')])

order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

print(order.counts)
# Archaeognatha     Blattodea    Coleoptera    Dermaptera       Diptera    Embioptera 
# 5            11            45          5285            14          1134             2 
# Ephemeroptera     Hemiptera   Hymenoptera      Isoptera   Lepidoptera      Mantodea     Mecoptera 
# 3          1080          3994            93          3441            29             1 
# Neuroptera       Odonata    Orthoptera      Phasmida  Phthiraptera      Psocodea  Siphonaptera 
# 39            93           210             2             2            36             2 
# Thysanoptera   Trichoptera     Zoraptera     Zygentoma   Megaloptera    Plecoptera Raphidioptera 
# 24             5            NA             2             1             6             1 

##%######################################################%##
#                                                          #
####            basic map of PREDICTS sites             ####
#                                                          #
##%######################################################%##


# plot the raster in ggplot
map.world <- map_data('world')

# map of sites
# p1 <- ggplot() +
#   geom_map(data = map.world, map = map.world,
#            aes(map_id = region),
#            fill = "grey", colour = "grey", linewidth = 0.2) +
#   geom_point(data = sites, aes(x = Longitude, y = Latitude), 
#              col = "#1E90FF", fill = "#104E8B", shape = 21) +
#   theme(axis.title = element_blank(), 
#         plot.background = element_blank(), 
#         panel.background = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank()) +
#   ggtitle("a.")

# 设置热带区域（南北回归线）
tropic_of_cancer <- 23.43656
tropic_of_capricorn <- -23.43656

p1 <- ggplot() +
  # 底图：世界地图
  geom_map(data = map.world, map = map.world,
           aes(map_id = region),
           fill = "grey95", colour = "grey70", linewidth = 0.2) +
  
  # 添加经纬度网格线（浅灰色）
  geom_hline(yintercept = seq(-90, 90, by = 30), color = "grey80", linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = seq(-180, 180, by = 30), color = "grey80", linewidth = 0.2, linetype = "dashed") +
  
  # 添加赤道线（稍暗一些）
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3, linetype = "dashed") +
  
  # 添加南北回归线（标记热带区域的边界）
  geom_hline(yintercept = tropic_of_cancer, color = "#FF7F50", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = tropic_of_capricorn, color = "#FF7F50", linewidth = 0.5, linetype = "dashed") +
  
  # 热带区域着色（添加透明的带状区域）
  annotate("rect", xmin = -180, xmax = 180, 
           ymin = tropic_of_capricorn, ymax = tropic_of_cancer,
           fill = "#FF7F5030", alpha = 0.1) +
  
  # 采样点（使用渐变颜色区分热带/非热带）
  geom_point(data = sites, 
             aes(x = Longitude, y = Latitude,
                 fill = abs(Latitude) <= 23.43656), # 基于纬度确定是否在热带
             color = "white", shape = 21, size = 2.5, stroke = 0.3) +
  
  # 设置点的填充颜色
  scale_fill_manual(values = c("FALSE" = "#1E90FF", "TRUE" = "#FF8C00"),
                    name = "Region",
                    labels = c("FALSE" = "Non-Tropical", "TRUE" = "Tropical")) +
  
  # 美化主题
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 9, color = "grey30"),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  # 限制地图范围（可选）
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-60, 90))

print(p1)


# save plot
ggsave(filename = paste0(outDir, "/PREDICTS_points_map.pdf"), height = 4, width = 8)


# 创建分面板图，每个面板显示一种土地利用类型
create_map_panel <- function(land_use_type) {
  subset_data <- sites[sites$UI2 == land_use_type, ]
  
  p <- ggplot() +
    # 底图：世界地图
    geom_map(data = map.world, map = map.world,
             aes(map_id = region),
             fill = "grey95", colour = "grey70", linewidth = 0.2) +
    # 添加经纬度网格线
    geom_hline(yintercept = seq(-90, 90, by = 30), color = "grey80", linewidth = 0.2, linetype = "dashed") +
    geom_vline(xintercept = seq(-180, 180, by = 30), color = "grey80", linewidth = 0.2, linetype = "dashed") +
    # 添加赤道线
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3, linetype = "dashed") +
    # 添加南北回归线
    geom_hline(yintercept = tropic_of_cancer, color = "#FF7F50", linewidth = 0.5, linetype = "dashed") +
    geom_hline(yintercept = tropic_of_capricorn, color = "#FF7F50", linewidth = 0.5, linetype = "dashed") +
    # 热带区域着色
    annotate("rect", xmin = -180, xmax = 180,
             ymin = tropic_of_capricorn, ymax = tropic_of_cancer,
             fill = "#FF7F5030", alpha = 0.1) +
    # 采样点（使用渐变颜色区分热带/非热带）
    geom_point(data = subset_data,
               aes(x = Longitude, y = Latitude,
                   fill = abs(Latitude) <= 23.43656), # 基于纬度确定是否在热带
               color = "white", shape = 21, size = 2, stroke = 0.3) +
    # 设置点的填充颜色
    scale_fill_manual(values = c("FALSE" = "#1E90FF", "TRUE" = "#FF8C00"),
                      name = "Region",
                      labels = c("FALSE" = "Non-Tropical", "TRUE" = "Tropical")) +
    # 美化主题
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 8, color = "grey30"),
      axis.title = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    ) +
    # 限制地图范围
    coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-60, 90)) +
    # 添加标题
    ggtitle(land_use_type) +
    # 添加站点数量信息
    labs(subtitle = paste("n =", nrow(subset_data), "sites"))
  
  return(p)
}

# 创建四个面板图
p2 <- create_map_panel("Primary vegetation")
p3 <- create_map_panel("Urban_Low")
p4 <- create_map_panel("Urban_Moderate") 
p5 <- create_map_panel("Urban_High")

# 创建共享的图例
legend_plot <- ggplot() +
  geom_point(data = data.frame(x = 1:2, y = 1:2, Region = c(TRUE, FALSE)),
             aes(x = x, y = y, fill = Region),
             color = "white", shape = 21, size = 3, stroke = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#FF8C00", "FALSE" = "#1E90FF"),
                    name = "Region",
                    labels = c("TRUE" = "Tropical", "FALSE" = "Non-Tropical")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# 提取图例
legend <- cowplot::get_legend(legend_plot)

# 合并所有面板
combined_plot <- gridExtra::grid.arrange(
  p2, p3, p4, p5, legend,
  layout_matrix = rbind(
    c(1, 2),
    c(3, 4),
    c(5, 5)
  ),
  heights = c(1, 1, 0.2)
)


# 保存最终图片
ggsave(filename = paste0(outDir, "/land_use_spatial_distribution.pdf"), width = 4, height = 8)


# 计算热带与非热带数据总量的简单代码
calculate_tropical_nontropical <- function(data) {
  # 定义热带区域边界
  tropic_of_cancer <- 23.43656
  tropic_of_capricorn <- -23.43656
  
  # 计算热带数据数量
  tropical_count <- sum(abs(data$Latitude) <= tropic_of_cancer)
  
  # 计算非热带数据数量
  non_tropical_count <- nrow(data) - tropical_count
  
  # 计算百分比
  total_count <- nrow(data)
  tropical_percent <- round(tropical_count / total_count * 100, 1)
  non_tropical_percent <- round(non_tropical_count / total_count * 100, 1)
  
  # 创建结果数据框
  result <- data.frame(
    Region = c("Tropical", "Non-tropical", "Total"),
    Count = c(tropical_count, non_tropical_count, total_count),
    Percentage = c(tropical_percent, non_tropical_percent, 100)
  )
  
  # 打印结果
  cat("热带与非热带数据数量:\n")
  cat("------------------------------\n")
  cat(paste("热带数据点: ", tropical_count, " (", tropical_percent, "%)\n", sep=""))
  cat(paste("非热带数据点: ", non_tropical_count, " (", non_tropical_percent, "%)\n", sep=""))
  cat(paste("总数据点: ", total_count, "\n", sep=""))
  
  return(result)
}

# 计算热带与非热带数据数量
tropical_stats <- calculate_tropical_nontropical(sites)


### Basic summaries ###

# nstudies/ nsites - all
length(unique(sites$SS)) # 183
length(unique(sites$SSBS)) # 2171

# nstudies/nsites - abun
length(unique(sites[!is.na(sites$LogAbund) , 'SS'])) # 243
length(unique(sites[!is.na(sites$LogAbund) , 'SSBS'])) # 3565

# reviewer request, UI2 by Biome

table(sites$Biome, sites$UI2)

#                                                         Agriculture_High Agriculture_Low Primary vegetation Secondary vegetation
# Boreal Forests/Taiga                                                    2               6                172                   13
# Temperate Conifer Forests                                               4             103                 10                   88
# Temperate Broadleaf & Mixed Forests                                  1308             671                315                  787
# Montane Grasslands & Shrublands                                         2             200                247                   33
# Temperate Grasslands, Savannas & Shrublands                            15              11                 15                   27
# Mediterranean Forests, Woodlands & Scrub                               21              32                 96                   58
# Deserts & Xeric Shrublands                                              0              30                 16                   16
# Tropical & Subtropical Grasslands, Savannas & Shrublands               56              47                175                   78
# Tropical & Subtropical Coniferous Forests                               2              26                 32                   43
# Flooded Grasslands & Savannas                                           6               6                  0                    0
# Tropical & Subtropical Dry Broadleaf Forests                           62              66                 13                   50
# Tropical & Subtropical Moist Broadleaf Forests                        293             119                420                  283
# Mangroves                                                               8               0                  5                    7

table(sites[!is.na(sites$LogAbund), 'Biome'], sites[!is.na(sites$LogAbund), 'UI2'])

#                                                          Agriculture_High Agriculture_Low Primary vegetation Secondary vegetation
# Boreal Forests/Taiga                                                    2               6                172                   13
# Temperate Conifer Forests                                               4             103                 10                   88
# Temperate Broadleaf & Mixed Forests                                  1280             669                289                  729
# Montane Grasslands & Shrublands                                         2             200                247                   33
# Temperate Grasslands, Savannas & Shrublands                            15              11                 15                   27
# Mediterranean Forests, Woodlands & Scrub                               21              26                 92                   58
# Deserts & Xeric Shrublands                                              0              30                 16                   16
# Tropical & Subtropical Grasslands, Savannas & Shrublands               56              47                165                   70
# Tropical & Subtropical Coniferous Forests                               2              26                 32                   43
# Flooded Grasslands & Savannas                                           0               0                  0                    0
# Tropical & Subtropical Dry Broadleaf Forests                           62              66                 13                   50
# Tropical & Subtropical Moist Broadleaf Forests                        265             110                354                  204
# Mangroves                                                               8               0                  5                    7

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
