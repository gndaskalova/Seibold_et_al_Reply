# Figures for the reply to Seibold et al.

# Note - the model scripts have to be run first to generate the objects needed here

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

# Theme ----
theme_inv <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14),
          axis.line.x = element_line(color = "black"), 
          axis.line.y = element_line(color = "black"),
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
load("data/biomass_terr_inv_m.RData")
load("data/biomass_fresh_inv_m.RData")
load("data/abundance_fresh_inv_m.RData")
load("data/abundance_terr_inv_m.RData") # Too big for GitHub
load("data/richness_terr_inv_m.RData")
load("data/richness_fresh_inv_m.RData")
load("data/duration_abundance2.RData")
load("data/duration_abundance.RData")
load("data/duration_biomass.RData")
load("data/duration_biomass2.RData")
load("rarefied_medians2018.Rdata") # Too big for GitHub

# ** Biomass ----
# Extract slopes and 95% CIs for each time series
biomass_fresh_eff <- as.data.frame(biomass_fresh_inv$Sol)[85:92]
biomass_fresh_eff <- biomass_fresh_eff %>% gather(rarefyID, slope, 1:8) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.03067))
biomass_fresh_eff$taxa <- "Freshwater"
biomass_fresh_intervals <- as.data.frame(HPDinterval(biomass_fresh_inv$Sol + biomass_fresh_inv$Sol[, 'year2']))
biomass_fresh_intervals <- biomass_fresh_intervals[85:92,]
biomass_fresh_intervals$rarefyID <- rownames(biomass_fresh_intervals)
biomass_fresh_eff <- left_join(biomass_fresh_eff, biomass_fresh_intervals)

biomass_terr_eff <- as.data.frame(biomass_terr_inv$Sol)[50:74]
biomass_terr_eff <- biomass_terr_eff %>% gather(rarefyID, slope, 1:25) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + 0.02255))
biomass_terr_eff$taxa <- "Terrestrial"
biomass_terr_intervals <- as.data.frame(HPDinterval(biomass_terr_inv$Sol + biomass_terr_inv$Sol[, 'year2']))
biomass_terr_intervals <- biomass_terr_intervals[50:74,]
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

(figure2a <- ggplot() +
    geom_point(data = biomass_effs, aes(x = duration, y = slope,
                                        colour = taxa, shape = category), size = 3, alpha = 0.5) +
    geom_point(data = biomass_effs[biomass_effs$category == "non-significant",], aes(x = duration, y = slope,
                                                                                     colour = taxa, shape = category), size = 3.2, alpha = 0.5) +
    geom_point(data = biomass_effs[biomass_effs$category == "non-significant",], aes(x = duration, y = slope,
                                                                                     colour = taxa, shape = category), size = 3.4, alpha = 0.5) +
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
    scale_shape_manual(values = c(21, 19)) +
    scale_x_continuous(limits = c(0, 101)) +
    scale_colour_manual(values = c("#59BAC0", "#d8b70a")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_y_continuous(limits = c(-0.2, 0.2)) +
    guides(colour = FALSE) +
    labs(x = "\nDuration (years)", y = "log(Biomass) / year\n", title = "a"))

ggsave(figure2a, filename = "figures/figure2a.png", height = 5, width = 5)

# ** Abundance ----

# Extract slopes and upper and lower 95% CI

# Freshwater
abundance_fresh_eff <- as.data.frame(abundance_fresh_inv$Sol)[207:278]
abundance_fresh_eff <- abundance_fresh_eff %>% gather(rarefyID, slope, 1:72) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.02934))
abundance_fresh_eff$taxa <- "Freshwater"
abundance_fresh_intervals <- as.data.frame(HPDinterval(abundance_fresh_inv$Sol + abundance_fresh_inv$Sol[, 'year2']))
abundance_fresh_intervals <- abundance_fresh_intervals[207:278,]
abundance_fresh_intervals$rarefyID <- rownames(abundance_fresh_intervals)
abundance_fresh_eff <- left_join(abundance_fresh_eff, abundance_fresh_intervals)
abundance_fresh_eff$rarefyID <- gsub("year2.rarefyID.", "", abundance_fresh_eff$rarefyID)

# Terrestrial
abundance_terr_eff <- as.data.frame(abundance_terr_inv$Sol)[308:397]
abundance_terr_eff <- abundance_terr_eff %>% gather(rarefyID, slope, 1:90) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.01051))
abundance_terr_eff$taxa <- "Terrestrial"
abundance_terr_intervals <- as.data.frame(HPDinterval(abundance_terr_inv$Sol + abundance_terr_inv$Sol[, 'year2']))
abundance_terr_intervals <- abundance_terr_intervals[308:397,]
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

(figure2b <- ggplot() +
    geom_point(data = abundance_effs, aes(x = duration, y = slope,
                                          colour = taxa, shape = category), size = 3, alpha = 0.5) +
    geom_point(data = abundance_effs[richness_effs$category == "non-significant",], aes(x = duration, y = slope,
                                                                                        colour = taxa, shape = category), size = 3.2, alpha = 0.5) +
    geom_point(data = abundance_effs[richness_effs$category == "non-significant",], aes(x = duration, y = slope,
                                                                                        colour = taxa, shape = category), size = 3.4, alpha = 0.5) +
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
    scale_shape_manual(values = c(21, 19)) +
    scale_x_continuous(limits = c(0, 101)) +
    scale_colour_manual(values = c("#59BAC0", "#d8b70a")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_y_continuous(limits = c(-0.2, 0.2)) +
    guides(colour = FALSE, shape = FALSE) +
    labs(x = "\nDuration (years)", y = "log(Abundance) / year\n", title = "b"))

# Note that the message about the one point excluded is the Lister and Garcia
# effect size which is too much of an outlier and doesn't fit on the graph

ggsave(figure2b, filename = "figures/figure2b.png", height = 5, width = 5)

# ** Species richness ----

# Extract slopes and 95% CIs for each time series

# Freshwater
richness_fresh_eff <- as.data.frame(richness_fresh_inv$Sol)[66:80]
richness_fresh_eff <- richness_fresh_eff %>% gather(rarefyID, slope, 1:15) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + -0.008648))
richness_fresh_eff$taxa <- "Freshwater"
richness_fresh_intervals <- as.data.frame(HPDinterval(richness_fresh_inv$Sol + richness_fresh_inv$Sol[, 'year2']))
richness_fresh_intervals <- richness_fresh_intervals[66:80,]
richness_fresh_intervals$rarefyID <- rownames(richness_fresh_intervals)
richness_fresh_eff <- left_join(richness_fresh_eff, richness_fresh_intervals)
richness_fresh_eff$rarefyID <- gsub("year2.rarefyID.", "", richness_fresh_eff$rarefyID)

# Terrestrial
richness_terr_eff <- as.data.frame(richness_terr_inv$Sol)[298:377]
richness_terr_eff <- richness_terr_eff %>% gather(rarefyID, slope, 1:80) %>%
  group_by(rarefyID) %>% summarise(slope = mean(slope + 0.007939))
richness_terr_eff$taxa <- "Terrestrial"
richness_terr_intervals <- as.data.frame(HPDinterval(richness_terr_inv$Sol + richness_terr_inv$Sol[, 'year2']))
richness_terr_intervals <- richness_terr_intervals[298:377,]
richness_terr_intervals$rarefyID <- rownames(richness_terr_intervals)
richness_terr_eff <- left_join(richness_terr_eff, richness_terr_intervals)
richness_terr_eff$rarefyID <- gsub("year2.rarefyID.", "", richness_terr_eff$rarefyID)

richness_effs <- bind_rows(richness_fresh_eff, richness_terr_eff)

# Add duration
duration <- rarefied_medians2018 %>% dplyr::select(rarefyID, duration) %>% distinct()
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
