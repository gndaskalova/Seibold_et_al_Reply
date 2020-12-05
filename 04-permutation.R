# Permutation analyses
# Code by Ally Phillimore

library(tidyverse)
library(lme4)
library(lmerTest)
library(scales)

# Data ----
seibold <- read.csv("data/25786.txt", sep = "\t", header = TRUE)

seibold$yearfactor <- as.factor(seibold$CollectionYear)
seibold$yearcenter <- seibold$CollectionYear-2012
seibold$region <- as.factor(substring(as.character(seibold$PlotID), 1, 2))
seibold$log_biomass <- log(seibold$biomass)

forest_all <- subset(seibold, Habitat_type == "forest")
forest<-subset(forest_all, Sampling_regime=="annual")
grass <- subset(seibold, Habitat_type == "grassland")


npermute<-1000


biomass_grass1 <- lmer(log_biomass ~ yearcenter + region + (1|PlotID), data = grass)
biomass_grass_coef<-coefficients(summary(biomass_grass1))[2,1]

nullcoef_biomass_grass <-c()

for (x in 1:npermute){
  grass$yearpermute<-sample(1:10)[grass$yearfactor]
  
  biomass_grass_random <- lmer(log_biomass ~ yearpermute + region + (1|PlotID), data = grass)
  nullcoef_biomass_grass[x]<-coefficients (summary(biomass_grass_random))[2,1]
}

2*length(which(nullcoef_biomass_grass<= biomass_grass_coef))/npermute
#two tailed p value 0.096


######
biomass_forest1 <- lmer(log_biomass ~ yearcenter + region + (1|PlotID), data = forest_all)
biomass_forest_coef<-coefficients(summary(biomass_forest1))[2,1]

nullcoef_biomass_forest <-c()

for (x in 1:npermute){
  forest_all$yearpermute<-sample(1:9)[forest_all$yearfactor]
  
  biomass_forest_random <- lmer(log_biomass ~ yearpermute + region + (1|PlotID), data = forest_all)
  nullcoef_biomass_forest[x]<-coefficients (summary(biomass_forest_random))[2,1]
}

2*length(which(nullcoef_biomass_forest<= biomass_forest_coef))/npermute
#two tailed p value 0.3684

#########
#0.46

abundance_grass1 <- glmer(abundance_identified ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
abundance_grass_coef<-coefficients(summary(abundance_grass1))[2,1]

nullcoef_abundance_grass <-c()

for (x in 1:npermute){
  grass$yearpermute<-sample(1:10)[grass$yearfactor]
  
  abundance_grass_random <- glmer(abundance_identified ~ yearpermute + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)

  nullcoef_abundance_grass[x]<-coefficients (summary(abundance_grass_random))[2,1]
}

2*length(which(nullcoef_abundance_grass<= abundance_grass_coef))/npermute


#########
#0.036

abundance_forest1 <- glmer(abundance_identified ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
abundance_forest_coef<-coefficients(summary(abundance_forest1))[2,1]

nullcoef_abundance_forest<-c()

for (x in 1:npermute){
  forest_all$yearpermute<-sample(1:9)[forest_all$yearfactor]
  
  abundance_forest_random <- glmer(abundance_identified ~ yearpermute + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)

  nullcoef_abundance_forest[x]<-coefficients (summary(abundance_forest_random))[2,1]
}

2*length(which(nullcoef_abundance_forest<= abundance_forest_coef))/npermute

#########
#1

species_grass1 <- glmer(species ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
species_grass_coef<-coefficients(summary(species_grass1))[2,1]

nullcoef_species_grass <-c()

for (x in 1:npermute){
  grass$yearpermute<-sample(1:10)[grass$yearfactor]
  
 species_grass_random <- glmer(species ~ yearpermute + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)

  nullcoef_species_grass[x]<-coefficients (summary(species_grass_random))[2,1]
}

2*length(which(nullcoef_species_grass<= species_grass_coef))/npermute
#0.122

#########

species_forest1 <- glmer(species ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
species_forest_coef<-coefficients(summary(species_forest1))[2,1]

nullcoef_species_forest<-c()

for (x in 1:npermute){
  forest_all$yearpermute<-sample(1:9)[forest_all$yearfactor]
  
  species_forest_random <- glmer(species ~ yearpermute + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)

  nullcoef_species_forest[x]<-coefficients (summary(species_forest_random))[2,1]
}

2*length(which(nullcoef_species_forest<= species_forest_coef))/npermute
#0.424




par(mfrow=c(3,2))
hist(nullcoef_biomass_grass,main="Grass biomass",xlab="Coefficient")
abline(v=biomass_grass_coef,col=2,lwd=2)

hist(nullcoef_biomass_forest,main="Forest biomass, all years",xlab="Coefficient")
abline(v=biomass_forest_coef,col=2,lwd=2)

hist(nullcoef_abundance_grass,main="Grass abundance",xlab="Coefficient")
abline(v=abundance_grass_coef,col=2,lwd=2)

hist(nullcoef_abundance_forest,main="Forest abundance, all years",xlab="Coefficient")
abline(v=abundance_forest_coef,col=2,lwd=2)

hist(nullcoef_species_grass,main="Grass species richness",xlab="Coefficient")
abline(v=species_grass_coef,col=2,lwd=2)

hist(nullcoef_species_forest,main="Forest species richness, all years",xlab="Coefficient")
abline(v=species_forest_coef,col=2,lwd=2)
