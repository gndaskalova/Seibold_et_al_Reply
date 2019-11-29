# Invertebrate time series analysis

#  Load data ----
load("rarefied_mediansOct2017.Rdata")
load("inv_all_biomass26thNov_rarefy.RData")
load("inv_all_abundance26thNov_rarefy.RData")

# Libraries ----
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(ggalt)
library(lme4)
library(lmerTest)
library(MCMCglmm)

options(mc.cores = parallel::detectCores())

inv <- rarefied_medians %>% filter(TAXA %in% c("Terrestrial invertebrates",
                                               "Freshwater invertebrates")) %>% 
  filter(length(unique(YEAR)) > 3)

world <- map_data("world")

(inv.map <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=eck4") +
    theme_map() +
    geom_point(data = inv, 
               aes(x = rarefyID_x, y = rarefyID_y, colour = TAXA),
               alpha = 0.6, size = 2) +
    scale_colour_manual(values = c("#59BAC0", "#d8b70a")) +
    scale_y_continuous(limits = c(-80, 80)) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.justification = "top"))

ggsave(inv.map, filename = "figures/inv_map.png", height = 5, width = 8)

# Define priors ----
a <- 10000

pa_prior<-list(R=list(V=diag(1), nu=0.002), 
               G=list(	G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                       G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))

# Species richness models ----

# Freshwater
fresh.richness <- inv %>% filter(TAXA == "Freshwater invertebrates")
summary(fresh.richness$YEAR)  # starts in 1917, ends in 2014
fresh.richness$year2 <- fresh.richness$YEAR - 1994  # median of zero
summary(fresh.richness$year2)  # centered on zero

richness_fresh_inv <- MCMCglmm(S ~ year2, random = ~YEAR + us(1 + year2):rarefyID, 
                               prior = pa_prior, data = fresh.richness, 
                               family = "poisson", pr = TRUE, 
                               nitt = 100000, burnin = 10000)
summary(richness_fresh_inv)
plot(richness_fresh_inv)
save(richness_fresh_inv, file = "data/richness_fresh_inv_m.RData")

# Terrestrial
terr.richness <- inv %>% filter(TAXA == "Terrestrial invertebrates")
summary(terr.richness$YEAR)  # starts in 1898, ends in 2015
terr.richness$year2 <- terr.richness$YEAR - 1990  # median of zero
summary(terr.richness$year2)  # centered on zero

richness_terr_inv <- MCMCglmm(S ~ year2, random = ~YEAR + us(1 + year2):rarefyID, 
                              prior = pa_prior,  data = terr.richness, 
                              family = "poisson", pr = TRUE, 
                              nitt = 100000, burnin = 10000)

summary(richness_terr_inv)
plot(richness_terr_inv)
save(richness_terr_inv, file = "data/richness_terr_inv_m.RData")

# Biomass models ----

# Freshwater
fresh.biomass <- inv_all %>% filter(TAXA == "Freshwater invertebrates")
hist(fresh.biomass$log_biomass)
summary(fresh.biomass$YEAR)  # starts in 1969, ends in 2010
fresh.biomass$year2 <- fresh.biomass$YEAR - 1994  # median of zero
summary(fresh.biomass$year2)  # centered on zero


biomass_fresh_inv <- MCMCglmm(log_biomass ~ year2, random = ~YEAR + us(1 + year2):rarefyID, 
                              prior = pa_prior, data = fresh.biomass, 
                              pr = TRUE, nitt = 100000, burnin = 10000)
summary(biomass_fresh_inv)
plot(biomass_fresh_inv)
save(biomass_fresh_inv, file = "data/biomass_fresh_inv_m.RData")

# Calculate and save time series duration
duration_biomass2 <- fresh.biomass %>% group_by(rarefyID) %>%
  mutate(duration = max(YEAR) - min(YEAR)) %>%
  dplyr::select(rarefyID, duration) %>% distinct()
save(duration_biomass2, file = "data/duration_biomass2.RData")

# Terrestrial
terr.biomass <- inv_all %>% filter(TAXA == "Terrestrial invertebrates")
hist(terr.biomass$log_biomass)
summary(terr.biomass$YEAR)  # starts in 2004, ends in 2014
terr.biomass$year2 <- terr.biomass$YEAR - 2009  # median of zero
summary(terr.biomass$year2)  # centered on zero

biomass_terr_inv <- MCMCglmm(log_biomass ~ year2, random = ~YEAR + us(1 + year2):rarefyID, 
                             prior = pa_prior, data = terr.biomass, 
                             pr = TRUE, nitt = 100000, burnin = 10000)

summary(biomass_terr_inv)
plot(biomass_terr_inv)
save(biomass_terr_inv, file = "data/biomass_terr_inv_m.RData")

# Calculate and save time series duration
duration_biomass <- terr.biomass %>% group_by(rarefyID) %>%
  mutate(duration = max(YEAR) - min(YEAR)) %>%
  dplyr::select(rarefyID, duration) %>% distinct()
save(duration_biomass, file = "data/duration_biomass.RData")

# Abundance models ----

# Freshwater
fresh.abundance <- inv_all_abun %>% filter(TAXA == "Freshwater invertebrates")
hist(fresh.abundance$Abundance)
summary(fresh.abundance$YEAR)  # starts in 1916, ends in 2016
fresh.abundance$year2 <- fresh.abundance$YEAR - 1995  # median of zero
summary(fresh.abundance$year2)  # centered on zero

abundance_fresh_inv <- MCMCglmm(round(Abundance) ~ year2, random = ~YEAR + us(1 + year2):rarefyID, 
                                prior = pa_prior, data = fresh.abundance, 
                                family = "poisson", pr = TRUE, 
                                nitt = 100000, burnin = 10000)
summary(abundance_fresh_inv)
plot(abundance_fresh_inv)
save(abundance_fresh_inv, file = "data/abundance_fresh_inv_m.RData")

# Calculate and save time series duration
duration_abundance <- fresh.abundance %>% group_by(rarefyID) %>%
  mutate(duration = max(YEAR) - min(YEAR)) %>%
  dplyr::select(rarefyID, duration) %>% distinct()
save(duration_abundance, file = "data/duration_abundance.RData")

# Terrestrial
terr.abundance <- inv_all_abun %>% filter(TAXA == "Terrestrial invertebrates")
hist(terr.abundance$Abundance)
summary(terr.abundance$YEAR)  # starts in 1898, ends in 2015
terr.abundance$year2 <- terr.abundance$YEAR - 2003  # median of zero
summary(terr.abundance$year2)  # centered on zero

abundance_terr_inv <- MCMCglmm(round(Abundance) ~ year2, random = ~YEAR + us(1 + year2):rarefyID, 
                               prior = pa_prior, data = terr.abundance, 
                               family = "poisson", pr = TRUE, 
                               nitt = 100000, burnin = 10000)
summary(abundance_terr_inv)
plot(abundance_terr_inv)
save(abundance_terr_inv, file = "data/abundance_terr_inv_m.RData")

# Calculate and save time series duration
duration_abundance2 <- terr.abundance %>% group_by(rarefyID) %>%
  mutate(duration = max(YEAR) - min(YEAR)) %>%
  dplyr::select(rarefyID, duration) %>% distinct()
save(duration_abundance2, file = "data/duration_abundance2.RData")

# Model output table ----

# Function to extract MCMCglmm model summary outputs
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # adding the model name

# Load models
load("data/richness_terr_inv_m.RData")
load("data/richness_fresh_inv_m.RData")
load("data/abundance_terr_inv_m.RData")
load("data/abundance_fresh_inv_m.RData")
load("data/biomass_terr_inv_m.RData")
load("data/biomass_fresh_inv_m.RData")

# Creating a summary table of model outputs for models with random effects
dataList <- list(biomass_fresh_inv, biomass_terr_inv,
                 abundance_fresh_inv, abundance_terr_inv,
                 richness_fresh_inv, richness_terr_inv)

# Create a list of input model names
dataListNames <- list("Biomass (freshwater)", "Biomass (terrestrial)",
                      "Abundance (freshwater)", "Abundance (terrestrial)",
                      "Richness (freshwater)", "Richness (terrestrial)")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs <- as.data.frame(do.call(rbind, readyList))

# Write csv
write.csv(mcmc.outputs, file = "figures/biotime_mcmc_outputs_29thNov2019.csv")

# Tidy up table
colnames(mcmc.outputs)
mcmc.outputs <- mcmc.outputs %>%
  dplyr::select(modelName, variable,
                post.mean, l.95..CI,
                u.95..CI, eff.samp,
                pMCMC, effect) %>%
  mutate(variable = str_replace_all(variable, pattern = c('YEAR' = "year",
                                                  'year2' = "year",
                                                  'units' = 'sigma')))

mcmc.outputs$modelName[duplicated(mcmc.outputs$modelName)] <- " "

colnames(mcmc.outputs) <- c("Model", "Variable", "Post. mean", "Lowe 95% CI", "Upper 95% CI",
                            "Eff. sample", "pMCMC", "Effect")
library(stargazer)
stargazer(mcmc.outputs, type = "html", summary = FALSE, digits = 3)
