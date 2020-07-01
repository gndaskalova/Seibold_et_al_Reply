# Figures for the reply to Seibold et al.
# Gergana Daskalova
# gndaskalova@gmail.com

# Note - the model scripts have to be run first to generate the objects needed here
# or the model objects need to be loaded

# Libraries ----
library(tidyverse)
library(lme4)
library(MCMCglmm)
library(ggeffects)
library(lmerTest)
library(gridExtra)
library(ggrepel)
library(readr)
library(ggforce)
library(ggalt)
library(ggthemes)

# Theme ----
theme_inv <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1,
                                    hjust = -0.1, face = "bold", family = "Arial"),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 1), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

# Load model output table
modeltable <- read_delim("figures/modeltable_parsed.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

modeltable$response <- factor(modeltable$response,
                              levels = c("biomass", "abundance", "species richness"),
                              labels = c("Biomass", "Abundance", "Richness"))
modeltable$test <- modeltable$pvalues < 0.05
modeltable$modeltest <- paste0(modeltable$modeltype, modeltable$test)
modeltable_f <- modeltable %>% filter(habitat == "forest" & plots == 30)
modeltable_f2 <- modeltable %>% filter(habitat == "forest" & plots == 140)
modeltable_g <- modeltable %>% filter(habitat == "grassland")

# Figure 1 ----
(eff_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(data = modeltable_f, aes(x = response, ymin = coef - se,
                                           ymax = coef + se, group = factor(modeltype),
                                           linetype = test),
                  position = position_dodge(1), colour = "#FE8324", size = 0.8) +
    geom_point(data = modeltable_f, aes(x = response, y = coef, group = factor(modeltype),
                                        shape = modeltest, fill = test),
               position = position_dodge(1), 
               size = 5.8, colour = "#FE8324") +
    geom_point(data = modeltable_f, aes(x = response, y = coef, group = factor(modeltype),
                                        shape = modeltest, fill = test),
               position = position_dodge(1), 
               size = 5, colour = "#FE8324") +
    scale_shape_manual(values = c(19, 24, 17,
                                  22, 15)) +
    scale_fill_manual(values = c("white", "#FE8324")) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    guides(fill = FALSE, linetype = FALSE, shape = FALSE) +
    theme_inv() +
    scale_y_continuous(limits = c(-0.08, 0.015)) +
    labs(x = NULL, y = "Coefficient\n", title = "a   Forest (30 plots)\n"))

(eff_plot2 <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(data = modeltable_f2, aes(x = response, ymin = coef - se,
                                            ymax = coef + se, group = factor(modeltype),
                                            linetype = test),
                  position = position_dodge(1), colour = "#FE8324", size = 0.8) +
    geom_point(data = modeltable_f2, aes(x = response, y = coef, group = factor(modeltype),
                                         shape = modeltest, fill = test),
               position = position_dodge(1), size = 5.8, colour = "#FE8324") +
    geom_point(data = modeltable_f2, aes(x = response, y = coef, group = factor(modeltype),
                                         shape = modeltest, fill = test),
               position = position_dodge(1), size = 5, colour = "#FE8324") +
    scale_shape_manual(values = c(21, 19, 24,
                                  22, 15)) +
    scale_fill_manual(values = c("white", "#FE8324")) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    guides(fill = FALSE, linetype = FALSE, shape = FALSE) +
    theme_inv() +
    labs(x = NULL, y = "Coefficient\n", title = "b   Forest (140 plots)\n"))

(eff_plot3 <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(data = modeltable_g, aes(x = response, ymin = coef - se,
                                           ymax = coef + se, group = factor(modeltype),
                                           linetype = test),
                  position = position_dodge(1), colour = "#5597AF", size = 0.8) +
    geom_point(data = modeltable_g, aes(x = response, y = coef, group = factor(modeltype),
                                        shape = modeltest, fill = test),
               position = position_dodge(1), size = 5.8,
               colour = "#5597AF") +
    geom_point(data = modeltable_g, aes(x = response, y = coef, group = factor(modeltype),
                                        shape = modeltest, fill = test),
               position = position_dodge(1), size = 5,
               colour = "#5597AF") +
    scale_shape_manual(values = c(21, 24, 24,
                                  22, 15)) +
    scale_fill_manual(values = c("white", "#5597AF")) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    guides(fill = FALSE, linetype = FALSE, shape = FALSE) +
    theme_inv() +
    labs(x = NULL, y = "Coefficient\n", title = "c   Grassland (150 plots)\n"))

figure1 <- grid.arrange(eff_plot, eff_plot2, eff_plot3, ncol = 3)

# Make legend
triangle <- data.frame(x = c(10, 11, 12, 10), y = c(0.5, 2.5, 0.5, 0.5))

(legend <- ggplot() +
    geom_circle(aes(x0 = 1, y0 = 1.5, r = 1), color = "grey40", fill = "grey40") +
    geom_rect(aes(xmin = 22, xmax = 24, 
                  ymin = 0.5, ymax = 2.5), color = "grey40", fill = "grey40") +
    geom_polygon(data = triangle, aes(x = x, y = y), color = "grey40", fill = "grey40") +
    theme_void() +
    scale_y_continuous(limits = c(0, 3)) +
    scale_x_continuous(limits = c(0, 33.5)) +
    annotate("text", x = 6, y = 1.5, label = "No random\nyear effects", size = 5) +
    annotate("text", x = 16, y = 1.5, label = "Random year\nintercept", size = 5) +
    annotate("text", x = 30, y = 1.5, label = "Random year\nintercept and slope", size = 5) +
    coord_equal())

figure1_legend <- grid.arrange(legend, figure1, heights = c(0.09, 0.91))

ggsave(figure1_legend, filename = "figures/figure1_2ndDec.png", height = 6, width = 15)
ggsave(figure1_legend, filename = "figures/figure1_2ndDec.pdf",
       device = cairo_pdf, height = 6, width = 15)

# Figure 2 ----
# Load models and duration of time series
# Note that the model output files are too large to be stored on GitHub

setwd("~/Downloads/Seibold June 2020 models")
load("biomass_terr_inv_m24June.RData")
load("biomass_fresh_inv_m24June.RData")
load("abundance_fresh_inv_m24June.RData")
load("abundance_terr_inv_m24June.RData")
load("richness_terr_inv_m24June.RData")
load("richness_fresh_inv_m24June.RData")
load("duration_abundance2.RData")
load("duration_abundance.RData")
load("duration_biomass.RData")
load("duration_biomass2.RData")
load("biomass_fresh_inv_m_vanklink24June.RData")
load("abundance_fresh_inv_m_vanklink24June.RData")
load("abundance_terr_inv_m_vanklink24June.RData")
load("biomass_terr_inv_m_vanklink24June.RData")

setwd("~/LandUseHub/data/input")
load("rarefied_mediansOct2017.Rdata")

setwd("~/ArmageddonHub")
# ** Biomass ----
# Extract slopes and 95% CIs for each time series
biomass_fresh_eff <- as.data.frame(biomass_fresh_inv$Sol)[76:80]
biomass_fresh_eff <- biomass_fresh_eff %>% gather(rarefyID, slope, 1:5) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.022560))
biomass_fresh_eff$taxa <- "Freshwater"
biomass_fresh_intervals <- as.data.frame(HPDinterval(biomass_fresh_inv$Sol + biomass_fresh_inv$Sol[, 'year2']))
biomass_fresh_intervals <- biomass_fresh_intervals[76:80,]
biomass_fresh_intervals$rarefyID <- rownames(biomass_fresh_intervals)
biomass_fresh_eff <- left_join(biomass_fresh_eff, biomass_fresh_intervals)

biomass_terr_eff <- as.data.frame(biomass_terr_inv$Sol)[46:66]
biomass_terr_eff <- biomass_terr_eff %>% gather(rarefyID, slope, 1:21) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + 0.02336))
biomass_terr_eff$taxa <- "Terrestrial"
biomass_terr_intervals <- as.data.frame(HPDinterval(biomass_terr_inv$Sol + biomass_terr_inv$Sol[, 'year2']))
biomass_terr_intervals <- biomass_terr_intervals[46:66,]
biomass_terr_intervals$rarefyID <- rownames(biomass_terr_intervals)
biomass_terr_eff <- left_join(biomass_terr_eff, biomass_terr_intervals)

# Combine effect sizes from freshwater and terrestrial in one object
biomass_effs <- bind_rows(biomass_fresh_eff, biomass_terr_eff)
biomass_effs$rarefyID <- gsub("year2.rarefyID.", "", biomass_effs$rarefyID)

# Add duration
biomass_duration <- bind_rows(duration_biomass, duration_biomass2)
biomass_effs <- left_join(biomass_effs, biomass_duration, by = "rarefyID")

# Create a column for significant or not
biomass_effs <- biomass_effs %>%  mutate(category = with(., case_when(
  (lower < 0 & upper < 0) ~ 'significant',
  (lower > 0 & upper > 0) ~ 'significant',
  (lower < 0 & upper > 0) ~ 'non-significant')))

seibold_effs3 <- data.frame(9, -0.03625)
colnames(seibold_effs3) <- c("duration", "slope")
seibold_effs3$label <- "Seibold et al. 2019 (forest)"

seibold_effs4 <- data.frame(10, -8.950e-02)
colnames(seibold_effs4) <- c("duration", "slope")
seibold_effs4$label <- "Seibold et al. 2019 (grassland)"

# Add the rest of the effect sizes from published studies
published_effect_sizes <- read_delim("data/published_effect_sizes.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

biomass_effs_published <- published_effect_sizes %>% filter(Taxon %in% c("flying insects", "moths"))
biomass_effs_published$label <- NA
biomass_effs_published[1,11] <- "Hallmann et al. 2017"
biomass_effs_published[2,11] <- "Macgregor et al. 2019"

# Add van Klink et al. effect sizes
# Extract slopes and 95% CIs for each time series
biomass_fresh_eff_vanklink <- as.data.frame(biomass_fresh_inv2$Sol)[339:377]
biomass_fresh_eff_vanklink <- biomass_fresh_eff_vanklink %>% gather(rarefyID, slope, 1:39) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.023896))
biomass_fresh_eff_vanklink$taxa <- "Freshwater"
biomass_fresh_intervals_vanklink <- as.data.frame(HPDinterval(biomass_fresh_inv2$Sol + biomass_fresh_inv2$Sol[, 'year2']))
biomass_fresh_intervals_vanklink <- biomass_fresh_intervals_vanklink[339:377,]
biomass_fresh_intervals_vanklink$rarefyID <- rownames(biomass_fresh_intervals_vanklink)
biomass_fresh_eff_vanklink <- left_join(biomass_fresh_eff_vanklink, biomass_fresh_intervals_vanklink)

biomass_terr_eff_vanklink <- as.data.frame(biomass_terr_inv2$Sol)[135:152]
biomass_terr_eff_vanklink <- biomass_terr_eff_vanklink %>% gather(rarefyID, slope, 1:18) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.04130))
biomass_terr_eff_vanklink$taxa <- "Terrestrial"
biomass_terr_intervals_vanklink <- as.data.frame(HPDinterval(biomass_terr_inv2$Sol + biomass_terr_inv2$Sol[, 'year2']))
biomass_terr_intervals_vanklink <- biomass_terr_intervals_vanklink[135:152,]
biomass_terr_intervals_vanklink$rarefyID <- rownames(biomass_terr_intervals_vanklink)
biomass_terr_eff_vanklink <- left_join(biomass_terr_eff_vanklink, biomass_terr_intervals_vanklink)

# Combine effect sizes from freshwater and terrestrial in one object
biomass_effs_vanklink <- bind_rows(biomass_fresh_eff_vanklink, biomass_terr_eff_vanklink)
biomass_effs_vanklink$rarefyID <- gsub("year2.StudyDur_ID.", "", biomass_effs_vanklink$rarefyID)

# Add duration
vanKlink_S1 <- readRDS("data/aax9931-vanKlink-SM-Data-S1.rds")
vanKlink_S1$StudyDur_ID <- paste0(vanKlink_S1$Datasource_ID, vanKlink_S1$Duration)
duration_vk <- vanKlink_S1 %>% dplyr::select(StudyDur_ID, Duration) %>% distinct()
colnames(duration_vk)[1] <- "rarefyID"
str(duration_vk)
duration_vk$rarefyID <- as.character(duration_vk$rarefyID)
biomass_effs_vanklink <- left_join(biomass_effs_vanklink, duration_vk, by = "rarefyID")

# Create a column for significant or not
biomass_effs_vanklink <- biomass_effs_vanklink %>%  mutate(category = with(., case_when(
  (lower < 0 & upper < 0) ~ 'significant',
  (lower > 0 & upper > 0) ~ 'significant',
  (lower < 0 & upper > 0) ~ 'non-significant')))

colnames(biomass_effs)
colnames(biomass_effs_vanklink)

colnames(biomass_effs_vanklink)[6] <- "duration"

biomass_effs$dataset <- "BioTIME"
biomass_effs_vanklink$dataset <- "van Klink et al."

biomass_effs_all <- bind_rows(biomass_effs, biomass_effs_vanklink)
biomass_effs_all$cat2 <- paste0(biomass_effs_all$category, biomass_effs_all$dataset)

(figure2a <- ggplot() +
    geom_point(data = biomass_effs_all, aes(x = duration, y = slope,
                                        colour = taxa, shape = cat2), size = 3, alpha = 0.5) +
    geom_point(data = biomass_effs_all[biomass_effs_all$category == "non-significant",], aes(x = duration, y = slope,
                                                                                     colour = taxa, shape = cat2), 
               size = 3.2, alpha = 0.5) +
    geom_point(data = biomass_effs_all[biomass_effs_all$category == "non-significant",], aes(x = duration, y = slope,
                                                                                     colour = taxa, shape = cat2), 
               size = 3.4, alpha = 0.5) +
    geom_point(data = seibold_effs3, aes(x = duration, y = slope), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = seibold_effs3,
                     aes(x = duration, y = slope,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     nudge_x = 40,
                     nudge_y = -0.06,
                     min.segment.length = 1, inherit.aes = FALSE) +
    geom_point(data = seibold_effs4, aes(x = duration, y = slope), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = seibold_effs4,
                     aes(x = duration, y = slope,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     nudge_x = 30,
                     nudge_y = -0.05,
                     min.segment.length = 1, inherit.aes = FALSE) +
    geom_point(data = biomass_effs_published, aes(x = Duration, y = `Effect size per year (log scale)`), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = biomass_effs_published,
                     aes(x = Duration, y = `Effect size per year (log scale)`,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     nudge_x = 40,
                     nudge_y = 0.04,
                     min.segment.length = 1, inherit.aes = FALSE) +
    theme_inv() +
    guides(shape = FALSE) +
    scale_shape_manual(values = c(21, 2, 19, 17)) +
    scale_x_continuous(limits = c(0, 101)) +
    scale_colour_manual(values = c("#59BAC0", "#d8b70a")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_y_continuous(limits = c(-0.2, 0.2)) +
    guides(colour = FALSE) +
    labs(x = "\nDuration (years)", y = "log(Biomass) / year\n", title = "a"))

ggsave(figure2a, filename = "figures/figure2a_June.png", height = 5, width = 5)

# ** Abundance ----

# Extract slopes and upper and lower 95% CI

# Freshwater
abundance_fresh_eff <- as.data.frame(abundance_fresh_inv$Sol)[144:156]
abundance_fresh_eff <- abundance_fresh_eff %>% gather(rarefyID, slope, 1:13) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + 0.008035))
abundance_fresh_eff$taxa <- "Freshwater"
abundance_fresh_intervals <- as.data.frame(HPDinterval(abundance_fresh_inv$Sol + abundance_fresh_inv$Sol[, 'year2']))
abundance_fresh_intervals <- abundance_fresh_intervals[144:156,]
abundance_fresh_intervals$rarefyID <- rownames(abundance_fresh_intervals)
abundance_fresh_eff <- left_join(abundance_fresh_eff, abundance_fresh_intervals)
abundance_fresh_eff$rarefyID <- gsub("year2.rarefyID.", "", abundance_fresh_eff$rarefyID)

# Terrestrial
abundance_terr_eff <- as.data.frame(abundance_terr_inv$Sol)[286:356]
abundance_terr_eff <- abundance_terr_eff %>% gather(rarefyID, slope, 1:70) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.01011))
abundance_terr_eff$taxa <- "Terrestrial"
abundance_terr_intervals <- as.data.frame(HPDinterval(abundance_terr_inv$Sol + abundance_terr_inv$Sol[, 'year2']))
abundance_terr_intervals <- abundance_terr_intervals[286:356,]
abundance_terr_intervals$rarefyID <- rownames(abundance_terr_intervals)
abundance_terr_eff <- left_join(abundance_terr_eff, abundance_terr_intervals)
abundance_terr_eff$rarefyID <- gsub("year2.rarefyID.", "", abundance_terr_eff$rarefyID)

# Add duration
abundance_terr_eff <- left_join(abundance_terr_eff, duration_abundance2, by = "rarefyID")
abundance_fresh_eff <- left_join(abundance_fresh_eff, duration_abundance, by = "rarefyID")
abundance_effs <- bind_rows(abundance_terr_eff, abundance_fresh_eff)

# Create a column for significant or not
abundance_effs <- abundance_effs %>%  mutate(category = with(., case_when(
  (lower < 0 & upper < 0) ~ 'significant',
  (lower > 0 & upper > 0) ~ 'significant',
  (lower < 0 & upper > 0) ~ 'non-significant')))

# Add effect sizes from published studies
seibold_effs5 <- data.frame(9, -0.019977)
colnames(seibold_effs5) <- c("duration", "slope")
seibold_effs5$label <- "Seibold et al. 2019 (forest)"

seibold_effs6 <- data.frame(10, -0.1616735)
colnames(seibold_effs6) <- c("duration", "slope")
seibold_effs6$label <- "Seibold et al. 2019 (grassland)"

abundance_effs_published <- published_effect_sizes %>% filter(Taxon %in% c("walking sticks"))
abundance_effs_published$label <- "Lister and Garcia 2018"

# Add van Klink et al. effect sizes
# Extract slopes and 95% CIs for each time series
abundance_fresh_eff_vanklink <- as.data.frame(abundance_fresh_inv2$Sol)[619:697]
abundance_fresh_eff_vanklink <- abundance_fresh_eff_vanklink %>% gather(rarefyID, slope, 1:79) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + 0.001444))
abundance_fresh_eff_vanklink$taxa <- "Freshwater"
abundance_fresh_intervals_vanklink <- as.data.frame(HPDinterval(abundance_fresh_inv2$Sol + abundance_fresh_inv2$Sol[, 'year2']))
abundance_fresh_intervals_vanklink <- abundance_fresh_intervals_vanklink[619:697,]
abundance_fresh_intervals_vanklink$rarefyID <- rownames(abundance_fresh_intervals_vanklink)
abundance_fresh_eff_vanklink <- left_join(abundance_fresh_eff_vanklink, abundance_fresh_intervals_vanklink)

abundance_terr_eff_vanklink <- as.data.frame(abundance_terr_inv2$Sol)[675:895]
abundance_terr_eff_vanklink <- abundance_terr_eff_vanklink %>% gather(rarefyID, slope, 1:221) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.013943))
abundance_terr_eff_vanklink$taxa <- "Terrestrial"
abundance_terr_intervals_vanklink <- as.data.frame(HPDinterval(abundance_terr_inv2$Sol + abundance_terr_inv2$Sol[, 'year2']))
abundance_terr_intervals_vanklink <- abundance_terr_intervals_vanklink[675:895,]
abundance_terr_intervals_vanklink$rarefyID <- rownames(abundance_terr_intervals_vanklink)
abundance_terr_eff_vanklink <- left_join(abundance_terr_eff_vanklink, abundance_terr_intervals_vanklink)

# Combine effect sizes from freshwater and terrestrial in one object
abundance_effs_vanklink <- bind_rows(abundance_fresh_eff_vanklink, abundance_terr_eff_vanklink)
abundance_effs_vanklink$rarefyID <- gsub("year2.StudyDur_ID.", "", abundance_effs_vanklink$rarefyID)

# Add duration
abundance_effs_vanklink <- left_join(abundance_effs_vanklink, duration_vk, by = "rarefyID")

# Create a column for significant or not
abundance_effs_vanklink <- abundance_effs_vanklink %>%  mutate(category = with(., case_when(
  (lower < 0 & upper < 0) ~ 'significant',
  (lower > 0 & upper > 0) ~ 'significant',
  (lower < 0 & upper > 0) ~ 'non-significant')))

colnames(abundance_effs)
colnames(abundance_effs_vanklink)

colnames(abundance_effs_vanklink)[6] <- "duration"

abundance_effs$dataset <- "BioTIME"
abundance_effs_vanklink$dataset <- "van Klink et al."

abundance_effs_all <- bind_rows(abundance_effs, abundance_effs_vanklink)
abundance_effs_all$cat2 <- paste0(abundance_effs_all$category, abundance_effs_all$dataset)

(figure2b <- ggplot() +
    geom_point(data = abundance_effs_all, aes(x = duration, y = slope,
                                          colour = taxa, shape = cat2), size = 3, alpha = 0.5) +
    geom_point(data = abundance_effs_all[abundance_effs_all$category == "non-significant",], aes(x = duration, y = slope,
                                                                                        colour = taxa, shape = cat2), size = 3.2, alpha = 0.5) +
    geom_point(data = abundance_effs_all[abundance_effs_all$category == "non-significant",], aes(x = duration, y = slope,
                                                                                        colour = taxa, shape = cat2), size = 3.4, alpha = 0.5) +
    geom_point(data = seibold_effs5, aes(x = duration, y = slope), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = seibold_effs5,
                     aes(x = duration, y = slope,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     nudge_x = 40,
                     nudge_y = -0.06,
                     min.segment.length = 1, inherit.aes = FALSE) +
    geom_point(data = seibold_effs6, aes(x = duration, y = slope), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = seibold_effs6,
                     aes(x = duration, y = slope,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     nudge_x = 40,
                     nudge_y = -0.05,
                     min.segment.length = 1, inherit.aes = FALSE) +
    geom_point(data = abundance_effs_published, aes(x = Duration, y = `Effect size per year (log scale)`), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = abundance_effs_published,
                     aes(x = Duration, y = `Effect size per year (log scale)`,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     # nudge_x = 30,
                     #nudge_y = -0.15,
                     min.segment.length = 1, inherit.aes = FALSE) +
    theme_inv() +
    scale_shape_manual(values = c(21, 2, 19, 17)) +
    scale_x_continuous(limits = c(0, 101)) +
    scale_colour_manual(values = c("#59BAC0", "#d8b70a")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_y_continuous(limits = c(-0.2, 0.2)) +
    guides(colour = FALSE, shape = FALSE) +
    labs(x = "\nDuration (years)", y = "log(Abundance) / year\n", title = "b"))

# Note that the message about the one point excluded is the Lister and Garcia
# effect size which is too much of an outlier and doesn't fit on the graph

ggsave(figure2b, filename = "figures/figure2b_June.png", height = 5, width = 5)

# ** Species richness ----

# Extract slopes and 95% CIs for each time series

# Freshwater
richness_fresh_eff <- as.data.frame(richness_fresh_inv$Sol)[146:160]
richness_fresh_eff <- richness_fresh_eff %>% gather(rarefyID, slope, 1:15) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.008924))
richness_fresh_eff$taxa <- "Freshwater"
richness_fresh_intervals <- as.data.frame(HPDinterval(richness_fresh_inv$Sol + richness_fresh_inv$Sol[, 'year2']))
richness_fresh_intervals <- richness_fresh_intervals[146:160,]
richness_fresh_intervals$rarefyID <- rownames(richness_fresh_intervals)
richness_fresh_eff <- left_join(richness_fresh_eff, richness_fresh_intervals)
richness_fresh_eff$rarefyID <- gsub("year2.rarefyID.", "", richness_fresh_eff$rarefyID)

# Terrestrial
richness_terr_eff <- as.data.frame(richness_terr_inv$Sol)[298:377]
richness_terr_eff <- richness_terr_eff %>% gather(rarefyID, slope, 1:80) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + 0.008119))
richness_terr_eff$taxa <- "Terrestrial"
richness_terr_intervals <- as.data.frame(HPDinterval(richness_terr_inv$Sol + richness_terr_inv$Sol[, 'year2']))
richness_terr_intervals <- richness_terr_intervals[298:377,]
richness_terr_intervals$rarefyID <- rownames(richness_terr_intervals)
richness_terr_eff <- left_join(richness_terr_eff, richness_terr_intervals)
richness_terr_eff$rarefyID <- gsub("year2.rarefyID.", "", richness_terr_eff$rarefyID)

richness_effs <- bind_rows(richness_fresh_eff, richness_terr_eff)

# Add duration
duration <- rarefied_medians %>% dplyr::select(rarefyID, duration) %>% distinct()
richness_effs <- left_join(richness_effs, duration, by = "rarefyID")

# Create a column for significant or not
richness_effs <- richness_effs %>%  mutate(category = with(., case_when(
  (lower < 0 & upper < 0) ~ 'significant',
  (lower > 0 & upper > 0) ~ 'significant',
  (lower < 0 & upper > 0) ~ 'non-significant')))

# Add effect sizes from published studies
seibold_effs <- data.frame(9, -0.040391)
colnames(seibold_effs) <- c("duration", "slope")
seibold_effs$label <- "Seibold et al. 2019 (forest)"

seibold_effs2 <- data.frame(10, -0.034460)
colnames(seibold_effs2) <- c("duration", "slope")
seibold_effs2$label <- "Seibold et al. 2019 (grassland)"

(figure2c <- ggplot() +
    geom_point(data = richness_effs, aes(x = duration, y = slope,
                                         colour = taxa, shape = category), size = 3, alpha = 0.5) +
    geom_point(data = richness_effs[richness_effs$category == "non-significant",], aes(x = duration, y = slope,
                                                                                       colour = taxa, shape = category), size = 3.2, alpha = 0.5) +
    geom_point(data = richness_effs[richness_effs$category == "non-significant",], aes(x = duration, y = slope,
                                                                                       colour = taxa, shape = category), size = 3.4, alpha = 0.5) +
    geom_point(data = seibold_effs, aes(x = duration, y = slope), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = seibold_effs,
                     aes(x = duration, y = slope,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     nudge_x = 40,
                     nudge_y = -0.04,
                     min.segment.length = 1, inherit.aes = FALSE) +
    geom_point(data = seibold_effs2, aes(x = duration, y = slope), 
               colour = "red", size = 3, alpha = 1) +
    geom_label_repel(data = seibold_effs2,
                     aes(x = duration, y = slope,
                         label = label),
                     box.padding = 1, size = 3, 
                     segment.size = 0.25,
                     nudge_x = 60, 
                     #nudge_y = 1,
                     min.segment.length = 1, inherit.aes = FALSE) +
    theme_inv() +
    scale_shape_manual(values = c(21, 19)) +
    guides(shape = FALSE) +
    scale_x_continuous(limits = c(0, 101)) +
    scale_colour_manual(values = c("#59BAC0", "#d8b70a")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_y_continuous(limits = c(-0.1, 0.1)) +
    labs(x = "\nDuration (years)", y = "log(Species richness) / year\n", title = "c"))

ggsave(figure2c, filename = "figures/figure2c.png", height = 5, width = 5)

figure2_panel <- grid.arrange(figure2a, figure2b, figure2c, ncol = 3)
ggsave(figure2_panel, filename = "figures/Figure2.png", height = 5, width = 15)
ggsave(figure2_panel, filename = "figures/Figure2.pdf", height = 5, width = 15, device = cairo_pdf)

# Figure S1 ----
# Model predictions figure

# Biomass figure
biomass_forest1_predictions <- ggpredict(biomass_forest1, terms = "yearcenter")
biomass_forest1_predictions$model <- "Model 1 (30 plots) *"
biomass_forest2_predictions <- ggpredict(biomass_forest2, terms = "yearcenter")
biomass_forest2_predictions$model <- "Model 2 (30 plots)"
biomass_forest3_predictions <- ggpredict(biomass_forest3, terms = "yearcenter")
biomass_forest3_predictions$model <- "Model 3 (30 plots)"

# Merge
biomass_predictions_f <- bind_rows(biomass_forest1_predictions,
                                   biomass_forest2_predictions,
                                   biomass_forest3_predictions)

(biomass_forest_graph <-  ggplot() +
    geom_boxplot(data = forest, aes(x = CollectionYear, y = log(biomass), 
                                    group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#FE8324", alpha = 0.7, width = 0.3) +
    geom_line(data = biomass_predictions_f, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = biomass_predictions_f, aes(x = x + 2012, ymin = conf.low, 
                                                  ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    scale_colour_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    theme_inv() +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    labs(x = NULL, y = "log(biomass)\n", title = "a"))

biomass_forest_all1_predictions <- ggpredict(biomass_forest_all1, terms = "yearcenter")
biomass_forest_all1_predictions$model <- "Model 1 (140 plots) *"
biomass_forest_all2_predictions <- ggpredict(biomass_forest_all2, terms = "yearcenter")
biomass_forest_all2_predictions$model <- "Model 2 (140 plots)"
biomass_forest_all3_predictions <- ggpredict(biomass_forest_all3, terms = "yearcenter")
biomass_forest_all3_predictions$model <- "Model 3 (140 plots)"

# Merge
biomass_predictions_all_f <- bind_rows(biomass_forest_all1_predictions,
                                       biomass_forest_all2_predictions,
                                       biomass_forest_all3_predictions)

(biomass_forest_graph2 <-  ggplot() +
    geom_boxplot(data = forest_all, aes(x = CollectionYear, y = log(biomass), 
                                        group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#FE8324", alpha = 0.7, width = 0.3) +
    geom_line(data = biomass_predictions_all_f, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = biomass_predictions_all_f, aes(x = x + 2012, ymin = conf.low, 
                                                      ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    scale_colour_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    theme_inv() +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    labs(x = NULL, y = "log(biomass)\n", title = "d"))

biomass_grassland1_predictions <- ggpredict(biomass_grass1, terms = "yearcenter")
biomass_grassland1_predictions$model <- "Model 1 (150 plots) *"
biomass_grassland2_predictions <- ggpredict(biomass_grass2, terms = "yearcenter")
biomass_grassland2_predictions$model <- "Model 2 (150 plots)"
biomass_grassland3_predictions <- ggpredict(biomass_grass3, terms = "yearcenter")
biomass_grassland3_predictions$model <- "Model 3 (150 plots)"

# Merge
biomass_predictions_g <- bind_rows(biomass_grassland1_predictions,
                                   biomass_grassland2_predictions,
                                   biomass_grassland3_predictions)

(biomass_grassland_graph <-  ggplot() +
    geom_boxplot(data = grass, aes(x = CollectionYear, y = log(biomass), 
                                   group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#5597AF", alpha = 0.7, width = 0.3) +
    geom_line(data = biomass_predictions_g, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = biomass_predictions_g, aes(x = x + 2012, ymin = conf.low, 
                                                  ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#5597AF", "#3b6271", "#151c1f")) +
    scale_colour_manual(values = c("#5597AF", "#3b6271", "#151c1f")) +
    theme_inv() +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    scale_y_continuous(limits = c(4, 10)) +
    labs(x = NULL, y = "log(biomass)\n", title = "h"))

# Abundance figure ----
abundance_forest1_predictions <- ggpredict(abundance_forest1, terms = "yearcenter")
abundance_forest1_predictions$model <- "Model 1 (30 plots) *"
abundance_forest2_predictions <- ggpredict(abundance_forest2, terms = "yearcenter")
abundance_forest2_predictions$model <- "Model 2 (30 plots)"
abundance_forest3_predictions <- ggpredict(abundance_forest3, terms = "yearcenter")
abundance_forest3_predictions$model <- "Model 3 (30 plots)"

# Merge
abundance_predictions_f <- bind_rows(abundance_forest1_predictions,
                                     abundance_forest2_predictions,
                                     abundance_forest3_predictions)

(abundance_forest_graph <-  ggplot() +
    geom_boxplot(data = forest, aes(x = CollectionYear, y = abundance_identified, 
                                    group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#FE8324", alpha = 0.7, width = 0.3) +
    geom_line(data = abundance_predictions_f, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = abundance_predictions_f, aes(x = x + 2012, ymin = conf.low, 
                                                    ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    scale_colour_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    theme_inv() +
    scale_y_continuous(limits = c(0, 800)) +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    labs(x = NULL, y = "Abundance\n", title = "b"))

abundance_forest_all1_predictions <- ggpredict(abundance_forest_all1, terms = "yearcenter")
abundance_forest_all1_predictions$model <- "Model 1 (140 plots)"
abundance_forest_all2_predictions <- ggpredict(abundance_forest_all2, terms = "yearcenter")
abundance_forest_all2_predictions$model <- "Model 2 (140 plots)"
abundance_forest_all3_predictions <- ggpredict(abundance_forest_all3, terms = "yearcenter")
abundance_forest_all3_predictions$model <- "Model 3 (140 plots)"

# Merge
abundance_predictions_all_f <- bind_rows(abundance_forest_all1_predictions,
                                         abundance_forest_all2_predictions,
                                         abundance_forest_all3_predictions)

(abundance_forest_graph2 <-  ggplot() +
    geom_boxplot(data = forest, aes(x = CollectionYear, y = abundance_identified, 
                                    group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#FE8324", alpha = 0.7, width = 0.3) +
    geom_line(data = abundance_predictions_all_f, aes(x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = abundance_predictions_all_f, aes(x = x + 2012, ymin = conf.low, 
                                                        ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    scale_colour_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    theme_inv() +
    scale_y_continuous(limits = c(0, 800)) +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    labs(x = NULL, y = "Abundance\n", title = "e"))

abundance_grassland1_predictions <- ggpredict(abundance_grass1, terms = "yearcenter")
abundance_grassland1_predictions$model <- "Model 1 (150 plots) *"
abundance_grassland2_predictions <- ggpredict(abundance_grass2, terms = "yearcenter")
abundance_grassland2_predictions$model <- "Model 2 (150 plots) *"
abundance_grassland3_predictions <- ggpredict(abundance_grass3, terms = "yearcenter")
abundance_grassland3_predictions$model <- "Model 3 (150 plots) *"

# Merge
abundance_predictions_g <- bind_rows(abundance_grassland1_predictions,
                                     abundance_grassland2_predictions,
                                     abundance_grassland3_predictions)

(abundance_grassland_graph <-  ggplot() +
    geom_boxplot(data = grass, aes(x = CollectionYear, y = abundance_identified, 
                                   group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#5597AF", alpha = 0.7, width = 0.3) +
    geom_line(data = abundance_predictions_g, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = abundance_predictions_g, aes(x = x + 2012, ymin = conf.low, 
                                                    ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#5597AF", "#3b6271", "#151c1f")) +
    scale_colour_manual(values = c("#5597AF", "#3b6271", "#151c1f")) +
    theme_inv() +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    scale_y_continuous(limits = c(0, 800)) +
    labs(x = NULL,y = "Abundance\n", title = "i"))

# Richness figure ----
species_forest1_predictions <- ggpredict(species_forest1, terms = "yearcenter")
species_forest1_predictions$model <- "Model 1 (30 plots) *"
species_forest2_predictions <- ggpredict(species_forest2, terms = "yearcenter")
species_forest2_predictions$model <- "Model 2 (30 plots) *"
species_forest3_predictions <- ggpredict(species_forest3, terms = "yearcenter")
species_forest3_predictions$model <- "Model 3 (30 plots) *"

# Merge
species_predictions_f <- bind_rows(species_forest1_predictions,
                                   species_forest2_predictions,
                                   species_forest3_predictions)

(species_forest_graph <-  ggplot() +
    geom_boxplot(data = forest, aes(x = CollectionYear, y = species, 
                                    group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#FE8324", alpha = 0.7, width = 0.3) +
    geom_line(data = species_predictions_f, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = species_predictions_f, aes(x = x + 2012, ymin = conf.low, 
                                                  ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    scale_colour_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    theme_inv() +
    scale_y_continuous(limits = c(0, 150)) +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    labs(x = NULL, y = "Species number\n", title = "c"))

species_forest_all1_predictions <- ggpredict(species_forest_all1, terms = "yearcenter")
species_forest_all1_predictions$model <- "Model 1 (140 plots) *"
species_forest_all2_predictions <- ggpredict(species_forest_all2, terms = "yearcenter")
species_forest_all2_predictions$model <- "Model 2 (140 plots)"
species_forest_all3_predictions <- ggpredict(species_forest_all3, terms = "yearcenter")
species_forest_all3_predictions$model <- "Model 3 (140 plots)"

# Merge
species_predictions_all_f <- bind_rows(species_forest_all1_predictions,
                                       species_forest_all2_predictions,
                                       species_forest_all3_predictions)

(species_forest_graph2 <-  ggplot() +
    geom_boxplot(data = forest, aes(x = CollectionYear, y = species, 
                                    group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#FE8324", alpha = 0.7, width = 0.3) +
    geom_line(data = species_predictions_all_f, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = species_predictions_all_f, aes(x = x + 2012, ymin = conf.low, 
                                                      ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    scale_colour_manual(values = c("#FE8324", "#a3561e", "#2b1a0e")) +
    theme_inv() +
    scale_y_continuous(limits = c(0, 150)) +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    labs(x = NULL, y = "Species number\n", title = "e"))

species_grassland1_predictions <- ggpredict(species_grass1, terms = "yearcenter")
species_grassland1_predictions$model <- "Model 1 (150 plots) *"
species_grassland2_predictions <- ggpredict(species_grass2, terms = "yearcenter")
species_grassland2_predictions$model <- "Model 2 (150 plots)"
species_grassland3_predictions <- ggpredict(species_grass3, terms = "yearcenter")
species_grassland3_predictions$model <- "Model 3 (150 plots)"

# Merge
species_predictions_g <- bind_rows(species_grassland1_predictions,
                                   species_grassland2_predictions,
                                   species_grassland3_predictions)

(species_grassland_graph <-  ggplot() +
    geom_boxplot(data = grass, aes(x = CollectionYear, y = species, 
                                   group = factor(CollectionYear)), 
                 outlier.shape = " ", fill = "#5597AF", alpha = 0.7, width = 0.3) +
    geom_line(data = species_predictions_g, aes( x = x + 2012, y = predicted, colour = model),
              size = 1) +
    geom_ribbon(data = species_predictions_g, aes(x = x + 2012, ymin = conf.low, 
                                                  ymax = conf.high, fill = model),
                alpha = 0.2) +
    scale_fill_manual(values = c("#5597AF", "#3b6271", "#151c1f")) +
    scale_colour_manual(values = c("#5597AF", "#3b6271", "#151c1f")) +
    theme_inv() +
    scale_x_continuous(limits = c(2007.8, 2017.2),
                       breaks = c(2008, 2011, 2014, 2017)) +
    scale_y_continuous(limits = c(0, 70)) +
    labs(x = NULL, y = "Species number\n", title = "j"))

figure1SI <- grid.arrange(biomass_forest_graph, 
                          abundance_forest_graph, 
                          species_forest_graph,
                          biomass_forest_graph2,
                          abundance_forest_graph2,
                          species_forest_graph2,
                          biomass_grassland_graph,
                          abundance_grassland_graph,
                          species_grassland_graph,
                          ncol = 3)

(legendSI <- ggplot() +
    geom_rect(aes(xmin = 1, xmax = 2, 
                  ymin = 1, ymax = 2), color = "#FE8324", fill = "#FE8324") +
    geom_rect(aes(xmin = 8, xmax = 9, 
                  ymin = 1, ymax = 2), color = "#5597AF", fill = "#5597AF") +
    theme_void() +
    scale_y_continuous(limits = c(0, 3)) +
    scale_x_continuous(limits = c(0, 15)) +
    annotate("text", x = 4, y = 1.5, label = "Forest", size = 6) +
    annotate("text", x = 12, y = 1.5, label = "Grassland", size = 6) +
    coord_equal())

figure1SI_legend <- grid.arrange(legendSI, figure1SI, heights = c(0.08, 0.92))

ggsave(figure1SI_legend, filename = "figures/Figure1_SI.pdf",
       height = 15, width = 19, device = cairo_pdf)

ggsave(figure1SI_legend, filename = "figures/Figure1_SI.png",
       height = 15, width = 19)

# Map
# BioTIME
inv <- rarefied_medians %>% filter(TAXA %in% c("Terrestrial invertebrates",
                                               "Freshwater invertebrates")) %>% 
  filter(length(unique(YEAR)) > 4)

inv_simple <- inv %>% dplyr::select(rarefyID_x, rarefyID_y, TAXA) %>% distinct()
inv_simple$TAXA <- factor(inv_simple$TAXA, levels = c("Freshwater invertebrates",
                                                         "Terrestrial invertebrates"),
                             labels = c("Freshwater", "Terrestrial"))
# van Klink et al.
load("data/vanklink.RData")
`aax9931-vanKlink-SM-Data-S2` <- readRDS("data/aax9931-vanKlink-SM-Data-S2.rds")
meta_vk <- `aax9931-vanKlink-SM-Data-S2` %>% dplyr::select(Plot_ID, Datasource_ID, Realm)
meta_vk <- meta_vk %>% rename("DataSource_ID" = Datasource_ID)
vanklink <- left_join(vanklink, meta_vk, by = c("DataSource_ID", "Plot_ID"))
vanklink <- filter(vanklink, Duration > 4)
vanklink$StudyDur_ID <- paste0(vanklink$DataSource_ID, vanklink$Duration)

# Remove time series which are from BioTIME to avoid duplicates
bt_test <- `aax9931-vanKlink-SM-Data-S2` %>% filter(str_detect(Datasource_name, "BioTIME")) %>%
  dplyr::select(Datasource_ID) %>% distinct()

# Remove time series 63, 70, 249, 294, 301, 313 and 375

vanklink <- vanklink %>% filter(!DataSource_ID %in% c("63", "70", "249",
                                                      "294", "301", "313", "375"))
vanklink_simple <- vanklink %>% dplyr::select(Longitude, Latitude, Realm) %>% distinct()

world <- map_data("world")

(inv.map <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=eck4") +
    theme_map() +
    geom_point(data = inv_simple, 
               aes(x = rarefyID_x, y = rarefyID_y, colour = TAXA),
               alpha = 0.6, size = 2, shape = 19) +
    geom_point(data = vanklink_simple, 
               aes(x = Longitude, y = Latitude, colour = Realm),
               alpha = 0.6, size = 2, shape = 17) +
    scale_colour_manual(values = c("#59BAC0", "#d8b70a")) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.justification = "top"))

# The warning message is about time series which didn't have freshwater/terrestrial noted 
# those weren't used in the analyses

ggsave(inv.map, filename = "figures/inv_mapJune.pdf", height = 5, width = 8)
