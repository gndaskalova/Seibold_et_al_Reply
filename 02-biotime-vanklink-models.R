# Invertebrate time series analysis

#  Load data ----
load("rarefied_mediansOct2017.Rdata")
load("D:/Gergana_Daskalova/LandUseHub/inv_all_biomass26thNov_rarefy.RData")
load("D:/Gergana_Daskalova/LandUseHub/inv_all_abundance26thNov_rarefy.RData")
bioTIMEmetadata <- read.csv("bioTIMEmetadataSept18.csv")


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
library(stringr)

options(mc.cores = parallel::detectCores())

inv <- rarefied_medians %>% filter(TAXA %in% c("Terrestrial invertebrates",
                                               "Freshwater invertebrates")) %>% 
  filter(length(unique(YEAR)) > 4)

# Define priors ----
a <- 10000

pa_prior <- list(R = list(V = diag(1), nu = 0.002),
                 G=list(G1 = list(V = diag(1), nu = 1, alpha.mu = c(0), alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = c(0), alpha.V = diag(1)*a),
                        G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*a)))

# Species richness models ----

# Freshwater
fresh.richness <- inv %>% filter(TAXA == "Freshwater invertebrates")
summary(fresh.richness$YEAR)  # starts in 1917, ends in 2014
fresh.richness$year2 <- fresh.richness$YEAR - 1994  # median of zero
summary(fresh.richness$year2)  # centered on zero

richness_fresh_inv <- MCMCglmm(S ~ year2, random = ~YEAR + BIOME_MAP:YEAR + us(1 + year2):rarefyID,
                               prior = pa_prior, data = fresh.richness, 
                               family = "poisson", pr = TRUE, 
                               nitt = 100000, burnin = 10000)
summary(richness_fresh_inv)
plot(richness_fresh_inv)
save(richness_fresh_inv, file = "data/richness_fresh_inv_m24June.RData")

# Terrestrial
terr.richness <- inv %>% filter(TAXA == "Terrestrial invertebrates")
summary(terr.richness$YEAR)  # starts in 1898, ends in 2015
terr.richness$year2 <- terr.richness$YEAR - 1990  # median of zero
summary(terr.richness$year2)  # centered on zero

richness_terr_inv <- MCMCglmm(S ~ year2, random = ~YEAR + BIOME_MAP:YEAR + us(1 + year2):rarefyID,
                              prior = pa_prior,  data = terr.richness, 
                              family = "poisson", pr = TRUE, 
                              nitt = 100000, burnin = 10000)

summary(richness_terr_inv)
plot(richness_terr_inv)
save(richness_terr_inv, file = "data/richness_terr_inv_m24June.RData")

# Biomass models ----

# Add biome information
meta <- bioTIMEmetadata %>% dplyr::select(STUDY_ID, BIOME_MAP)
meta$STUDY_ID <- as.factor(as.character(meta$STUDY_ID))

inv_all <- inv_all %>% mutate(STUDY_ID = rarefyID) %>%
  separate(STUDY_ID, c("STUDY_ID", "NA"))
inv_all$STUDY_ID <- as.factor(as.character(inv_all$STUDY_ID))
inv_all <- left_join(inv_all, meta, by = "STUDY_ID")

inv_all_abun <- inv_all_abun %>% mutate(STUDY_ID = rarefyID) %>%
  separate(STUDY_ID, c("STUDY_ID", "NA"))
inv_all_abun$STUDY_ID <- as.factor(as.character(inv_all_abun$STUDY_ID))
inv_all_abun <- left_join(inv_all_abun, meta, by = "STUDY_ID")

# Freshwater
dur_meta <- rarefied_medians %>% dplyr::select(rarefyID, duration) %>% distinct()
inv_all <-left_join(inv_all, dur_meta, by = "rarefyID")
inv_all <- inv_all %>% filter(duration > 4)
fresh.biomass <- inv_all %>% filter(TAXA == "Freshwater invertebrates")
hist(fresh.biomass$log_biomass)
summary(fresh.biomass$YEAR)  # starts in 1969, ends in 2010
fresh.biomass$year2 <- fresh.biomass$YEAR - 1993  # median of zero
summary(fresh.biomass$year2)  # centered on zero
summary(fresh.biomass)

biomass_fresh_inv <- MCMCglmm(log_biomass ~ year2, random = ~YEAR + BIOME_MAP:YEAR + us(1 + year2):rarefyID,
                              prior = pa_prior, data = fresh.biomass, 
                              pr = TRUE, nitt = 100000, burnin = 10000)
summary(biomass_fresh_inv)
plot(biomass_fresh_inv)
save(biomass_fresh_inv, file = "data/biomass_fresh_inv_m24June.RData")

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

biomass_terr_inv <- MCMCglmm(log_biomass ~ year2, random = ~YEAR + BIOME_MAP:YEAR + us(1 + year2):rarefyID,
                             prior = pa_prior, data = terr.biomass, 
                             pr = TRUE, nitt = 100000, burnin = 10000)

summary(biomass_terr_inv)
plot(biomass_terr_inv)
save(biomass_terr_inv, file = "data/biomass_terr_inv_m24June.RData")

# Calculate and save time series duration
duration_biomass <- terr.biomass %>% group_by(rarefyID) %>%
  mutate(duration = max(YEAR) - min(YEAR)) %>%
  dplyr::select(rarefyID, duration) %>% distinct()
save(duration_biomass, file = "data/duration_biomass.RData")

# Abundance models ----
inv_all_abun <-left_join(inv_all_abun, dur_meta, by = "rarefyID")
inv_all_abun <- inv_all_abun %>% filter(duration > 4)

# Freshwater
fresh.abundance <- inv_all_abun %>% filter(TAXA == "Freshwater invertebrates")
hist(fresh.abundance$Abundance)
summary(fresh.abundance$YEAR)  # starts in 1917, ends in 2014
fresh.abundance$year2 <- fresh.abundance$YEAR - 1995  # median of zero
summary(fresh.abundance$year2)  # centered on zero

abundance_fresh_inv <- MCMCglmm(round(Abundance) ~ year2, random = ~YEAR + BIOME_MAP:YEAR + us(1 + year2):rarefyID,
                                prior = pa_prior, data = fresh.abundance, 
                                family = "poisson", pr = TRUE, 
                                nitt = 100000, burnin = 10000)
summary(abundance_fresh_inv)
plot(abundance_fresh_inv)
save(abundance_fresh_inv, file = "data/abundance_fresh_inv_m24June.RData")

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

abundance_terr_inv <- MCMCglmm(round(Abundance) ~ year2, random = ~YEAR + BIOME_MAP:YEAR + us(1 + year2):rarefyID, 
                               prior = pa_prior, data = terr.abundance, 
                               family = "poisson", pr = TRUE, 
                               nitt = 100000, burnin = 10000)
summary(abundance_terr_inv)
plot(abundance_terr_inv)
save(abundance_terr_inv, file = "data/abundance_terr_inv_m24June.RData")

# Calculate and save time series duration
duration_abundance2 <- terr.abundance %>% group_by(rarefyID) %>%
  mutate(duration = max(YEAR) - min(YEAR)) %>%
  dplyr::select(rarefyID, duration) %>% distinct()
save(duration_abundance2, file = "data/duration_abundance2.RData")

# Models with van Klink data ----
load("data/vanklink.RData")
`aax9931-vanKlink-SM-Data-S2` <- readRDS("D:/Gergana_Daskalova/ArmageddonHub/data/aax9931-vanKlink-SM-Data-S2.rds")
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

vanklink_terr_biomass <- vanklink %>% filter(Realm == "Terrestrial" & MetricAB == "biomass")

# Terrestrial
hist(vanklink_terr_biomass$Number)
vanklink_terr_biomass$log_biomass <- log(vanklink_terr_biomass$Number)
hist(vanklink_terr_biomass$log_biomass)
summary(vanklink_terr_biomass$Year)  # starts in 1972, ends in 2017
vanklink_terr_biomass$year2 <- vanklink_terr_biomass$Year - 2009 # mediam of zero
summary(vanklink_terr_biomass$year2)
summary(vanklink_terr_biomass)

vanklink_terr_biomass <- vanklink_terr_biomass %>% drop_na(log_biomass) %>%
  filter(is.finite(log_biomass))

length(unique(vanklink_terr_biomass$DataSource_ID))
length(unique(vanklink_terr_biomass$StudyDur_ID))

vanklink_terr_biomass <- vanklink_terr_biomass %>% dplyr::select(log_biomass, year2, Year, 
                                                                 WWFecoRegion, StudyDur_ID, DataSource_ID) %>%
  distinct()

biomass_terr_inv2 <- MCMCglmm(log_biomass ~ year2, random = ~Year + WWFecoRegion:Year + us(1 + year2):StudyDur_ID,
                              prior = pa_prior, data = vanklink_terr_biomass, 
                              pr = TRUE, nitt = 100000, burnin = 10000)
summary(biomass_terr_inv2)
plot(biomass_terr_inv2)
save(biomass_terr_inv2, file = "data/biomass_terr_inv_m_vanklink24June.RData")

# Freshwater
vanklink_fresh_biomass <- vanklink %>% filter(Realm == "Freshwater" & MetricAB == "biomass")
hist(vanklink_fresh_biomass$Number)
vanklink_fresh_biomass$log_biomass <- log(vanklink_fresh_biomass$Number)
hist(vanklink_fresh_biomass$log_biomass)
summary(vanklink_fresh_biomass$Year)  # starts in 1925, ends in 2017
vanklink_fresh_biomass$year2 <- vanklink_fresh_biomass$Year - 2000 # median of zero
summary(vanklink_fresh_biomass$year2)
vanklink_fresh_biomass <- vanklink_fresh_biomass %>% drop_na(log_biomass) %>% 
  filter(is.finite(log_biomass))

summary(vanklink_fresh_biomass)

length(unique(vanklink_fresh_biomass$DataSource_ID))
length(unique(vanklink_fresh_biomass$StudyDur_ID))

vanklink_fresh_biomass <- vanklink_fresh_biomass %>% dplyr::select(log_biomass, year2, Year, 
                                                                 WWFecoRegion, StudyDur_ID, DataSource_ID) %>%
  distinct()

biomass_fresh_inv2 <- MCMCglmm(log_biomass ~ year2, random = ~Year + WWFecoRegion:Year + us(1 + year2):StudyDur_ID,
                               prior = pa_prior, data = vanklink_fresh_biomass, 
                               pr = TRUE, nitt = 100000, burnin = 10000)
summary(biomass_fresh_inv2)
plot(biomass_fresh_inv2)
save(biomass_fresh_inv2, file = "data/biomass_fresh_inv_m_vanklink24June.RData")

terr.abundance_vanklink <- vanklink %>% filter(Realm == "Terrestrial" & MetricAB == "abundance")
hist(terr.abundance_vanklink$Number)
summary(terr.abundance_vanklink$Year)  # starts in 1959, ends in 2018
terr.abundance_vanklink$year2 <- terr.abundance_vanklink$Year - 2002  # median of zero
summary(terr.abundance_vanklink$year2)  # centered on zero

length(unique(terr.abundance_vanklink$DataSource_ID))
length(unique(terr.abundance_vanklink$StudyDur_ID))

terr.abundance_vanklink <- terr.abundance_vanklink %>% dplyr::select(Number, year2, Year, 
                                                                 WWFecoRegion, StudyDur_ID, DataSource_ID) %>%
  distinct()

abundance_terr_inv2 <- MCMCglmm(round(Number) ~ year2, random = ~Year + WWFecoRegion:Year + us(1 + year2):StudyDur_ID, 
                                prior = pa_prior, data = terr.abundance_vanklink, 
                                family = "poisson", pr = TRUE, 
                                nitt = 100000, burnin = 10000)
summary(abundance_terr_inv2)
plot(abundance_terr_inv2)
save(abundance_terr_inv2, file = "data/abundance_terr_inv_m_vanklin24June.RData")

fresh.abundance_vanklink <- vanklink %>% filter(Stratum == "Water" & MetricAB == "abundance")
hist(fresh.abundance_vanklink$Number)
summary(fresh.abundance_vanklink$Year)  # starts in 1939, ends in 2017
fresh.abundance_vanklink$year2 <- fresh.abundance_vanklink$Year - 1999  # median of zero
summary(fresh.abundance_vanklink$year2)  # centered on zero

length(unique(fresh.abundance_vanklink$DataSource_ID))
length(unique(fresh.abundance_vanklink$StudyDur_ID))

fresh.abundance_vanklink <- fresh.abundance_vanklink %>% dplyr::select(Number, year2, Year, 
                                                                     WWFecoRegion, StudyDur_ID, DataSource_ID) %>%
  distinct()

abundance_fresh_inv2 <- MCMCglmm(round(Number) ~ year2, random = ~Year + WWFecoRegion:Year + us(1 + year2):StudyDur_ID, 
                                 prior = pa_prior, data = fresh.abundance_vanklink, 
                                 family = "poisson", pr = TRUE, 
                                 nitt = 100000, burnin = 10000)
summary(abundance_fresh_inv2)
plot(abundance_fresh_inv2)
save(abundance_fresh_inv2, file = "data/abundance_fresh_inv_m_vanklink24June.RData")

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

# Load models if previous code not run
load("richness_terr_inv_m24June.RData")
load("richness_fresh_inv_m24June.RData")
load("abundance_terr_inv_m24June.RData")
load("abundance_fresh_inv_m24June.RData")
load("biomass_terr_inv_m24June.RData")
load("biomass_fresh_inv_m24June.RData")
load("abundance_terr_inv_m_vanklink24June.RData")
load("abundance_fresh_inv_m_vanklink24June.RData")
load("biomass_terr_inv_m_vanklink24June.RData")
load("biomass_fresh_inv_m_vanklink24June.RData")

# Creating a summary table of model outputs for models with random effects
dataList <- list(biomass_fresh_inv, biomass_terr_inv,
                 abundance_fresh_inv, abundance_terr_inv,
                 richness_fresh_inv, richness_terr_inv,
                 biomass_fresh_inv2, biomass_terr_inv2,
                 abundance_fresh_inv2, abundance_terr_inv2)

# Create a list of input model names
dataListNames <- list("Biomass (freshwater, BioTIME)", "Biomass (terrestrial, BioTIME)",
                      "Abundance (freshwater, BioTIME)", "Abundance (terrestrial, BioTIME)",
                      "Richness (freshwater, BioTIME)", "Richness (terrestrial, BioTIME)",
                      "Biomass (freshwater, van Klink et al.)", "Biomass (terrestrial, van Klink et al.)",
                      "Abundance (freshwater, van Klink et al.", "Abundance (terrestrial, van Klink et al.)")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs <- as.data.frame(do.call(rbind, readyList))

# Write csv
write.csv(mcmc.outputs, file = "figures/biotime_vanklink_mcmc_outputs_26June2020.csv")

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

colnames(mcmc.outputs) <- c("Model", "Variable", "Post. mean", "Lower 95% CI", "Upper 95% CI",
                            "Eff. sample", "pMCMC", "Effect")
library(stargazer)
stargazer(mcmc.outputs, type = "html", summary = FALSE, digits = 3)
