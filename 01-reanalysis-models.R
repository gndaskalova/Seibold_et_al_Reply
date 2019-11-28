# Reanalysis of Seibold et al. 2019
# Data available from Seibold et al. 2019


# Libraries ----
library(tidyverse)
library(lme4)
library(lmerTest)

# Data ----
seibold <- read.csv("data/25786.txt", sep = "\t", header = TRUE)

seibold$yearfactor <- as.factor(seibold$CollectionYear)
seibold$yearcenter <- seibold$CollectionYear-2012
seibold$region <- as.factor(substring(as.character(seibold$PlotID), 1, 2))
seibold$log_biomass <- log(seibold$biomass)

forest_all <- subset(seibold, Habitat_type == "forest")
forest<-subset(forest_all, Sampling_regime=="annual")
grass <- subset(seibold, Habitat_type == "grassland")

# Workflow ----

#In the models below 
#1 = an approximation of the analysis in Seibold, focusing solely on temporal trends.
#2 = addition of a year random term
#3 = addition of a random slope term across plots.

# Models in the next section below are run using the lme4 package as used by Seibold et al.

# Biomass model ----

# Forests
biomass_forest1 <- lmer(log_biomass ~ yearcenter + region + (1|PlotID), data = forest)
summary(biomass_forest1)

biomass_forest2 <- lmer(log_biomass ~ yearcenter + region + (1|yearfactor) + (1|PlotID), data = forest)
summary(biomass_forest2)

biomass_forest3 <- lmer(log_biomass ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID), data = forest)
summary(biomass_forest3)

biomass_forest_all1 <- lmer(log_biomass ~ yearcenter + region + (1|PlotID), data = forest_all)
summary(biomass_forest_all1)

biomass_forest_all2 <- lmer(log_biomass ~ yearcenter + region + (1|yearfactor) + (1|PlotID), data = forest_all)
summary(biomass_forest_all2)

biomass_forest_all3 <- lmer(log_biomass ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID), data = forest_all)
summary(biomass_forest_all3)

# Grasslands
biomass_grass1 <- lmer(log_biomass ~ yearcenter + region + (1|PlotID), data = grass)
summary(biomass_grass1)

biomass_grass2 <- lmer(log_biomass ~ yearcenter + region + (1|yearfactor) + (1|PlotID), data = grass)
summary(biomass_grass2)

biomass_grass3 <- lmer(log_biomass ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID), data = grass)
summary(biomass_grass3)


# Abundance model ----

# Forests
abundance_forest1 <- glmer(abundance_identified ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest)
summary(abundance_forest1)

abundance_forest2 <- glmer(abundance_identified ~ yearcenter + region + (1|yearfactor) + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest)
summary(abundance_forest2)

abundance_forest3 <- glmer(abundance_identified ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID) + (1| PlotIDYear), family = "poisson", data = forest)
summary(abundance_forest3)

abundance_forest_all1 <- glmer(abundance_identified ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
summary(abundance_forest_all1)

abundance_forest_all2 <- glmer(abundance_identified ~ yearcenter + region + (1|yearfactor) + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
summary(abundance_forest_all2)

abundance_forest_all3 <- glmer(abundance_identified ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
summary(abundance_forest_all3)

# Grasslands
abundance_grass1 <- glmer(abundance_identified ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
summary(abundance_grass1)

abundance_grass2 <- glmer(abundance_identified ~ yearcenter + region + (1|yearfactor) + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
summary(abundance_grass2)

abundance_grass3 <- glmer(abundance_identified ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
summary(abundance_grass3)


# Species richness model ----

# Forests
species_forest1 <- glmer(species ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest)
summary(species_forest1)

species_forest2 <- glmer(species ~ yearcenter + region + (1|yearfactor) + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest)
summary(species_forest2)

species_forest3 <- glmer(species ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID) + (1| PlotIDYear), family = "poisson", data = forest)
summary(species_forest3)

species_forest_all1 <- glmer(species ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
summary(species_forest_all1)

species_forest_all2 <- glmer(species ~ yearcenter + region + (1|yearfactor) + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
summary(species_forest_all2)

species_forest_all3 <- glmer(species ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID) + (1| PlotIDYear), family = "poisson", data = forest_all)
summary(species_forest_all3)

# Grasslands
species_grass1 <- glmer(species ~ yearcenter + region + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
summary(species_grass1)

species_grass2 <- glmer(species ~ yearcenter + region + (1|yearfactor) + (1|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
summary(species_grass2)

species_grass3 <- glmer(species ~ yearcenter + region + (1|yearfactor) + (1 + yearcenter|PlotID) + (1| PlotIDYear), family = "poisson", data = grass)
summary(species_grass3)

# Model summary table ----
coeflist<-list()

coeflist$coefsbf1<-summary(biomass_forest1)$coefficients
coeflist$coefsbf2<-summary(biomass_forest2)$coefficients
coeflist$coefsbf3<-summary(biomass_forest3)$coefficients
coeflist$coefsbfall1<-summary(biomass_forest_all1)$coefficients
coeflist$coefsbfall2<-summary(biomass_forest_all2)$coefficients
coeflist$coefsbfall3<-summary(biomass_forest_all3)$coefficients
coeflist$coefsbg1<-summary(biomass_grass1)$coefficients
coeflist$coefsbg2<-summary(biomass_grass2)$coefficients
coeflist$coefsbg3<-summary(biomass_grass3)$coefficients

coeflist$coefsaf1<-summary(abundance_forest1)$coefficients
coeflist$coefsaf2<-summary(abundance_forest2)$coefficients
coeflist$coefsaf3<-summary(abundance_forest3)$coefficients
coeflist$coefsafall1<-summary(abundance_forest_all1)$coefficients
coeflist$coefsafall2<-summary(abundance_forest_all2)$coefficients
coeflist$coefsafall3<-summary(abundance_forest_all3)$coefficients
coeflist$coefsag1<-summary(abundance_grass1)$coefficients
coeflist$coefsag2<-summary(abundance_grass2)$coefficients
coeflist$coefsag3<-summary(abundance_grass3)$coefficients

coeflist$coefssf1<-summary(species_forest1)$coefficients
coeflist$coefssf2<-summary(species_forest2)$coefficients
coeflist$coefssf3<-summary(species_forest3)$coefficients
coeflist$coefssfall1<-summary(species_forest_all1)$coefficients
coeflist$coefssfall2<-summary(species_forest_all2)$coefficients
coeflist$coefssfall3<-summary(species_forest_all3)$coefficients
coeflist$coefssg1<-summary(species_grass1)$coefficients
coeflist$coefssg2<-summary(species_grass2)$coefficients
coeflist$coefssg3<-summary(species_grass3)$coefficients

# Extract slope and se, proportional change, P value for each model

extractslopeandse<-function(x){return(paste(round(x["yearcenter", "Estimate"],3),"+/-",round(x["yearcenter", "Std. Error"],3)))}
extractpropdecline<-function(x){return(round(exp(x["yearcenter", "Estimate"]),3))}
extractpvslmer<-function(x){return(round(x["yearcenter", 5],3))}
extractpvsglmer<-function(x){return(round(x["yearcenter", 4],3))}

coefsandse<-unlist(lapply(coeflist,extractslopeandse))
propchange<-unlist(lapply(coeflist,extractpropdecline))
pvalues<-c(unlist(lapply(coeflist[1:9],extractpvslmer)),unlist(lapply(coeflist[10:27],extractpvsglmer)))

# Extract year variances and the significance of the year terms.

year_variance<-rep("-",27)
lik.ratio<-rep("-",27)
randomterm_pval<-rep("-",27)

#extract year variances for model types 2 and 3.
bfv2<-as.data.frame(VarCorr(biomass_forest2))
year_variance[2]<-round(bfv2[which(bfv2$grp=="yearfactor"),"vcov"],3)
bfv3<-as.data.frame(VarCorr(biomass_forest3))
year_variance[3]<-round(bfv3[which(bfv3$grp=="yearfactor"),"vcov"],3)

bfv2all<-as.data.frame(VarCorr(biomass_forest_all2))
year_variance[5]<-round(bfv2all[which(bfv2all$grp=="yearfactor"),"vcov"],3)
bfv3all<-as.data.frame(VarCorr(biomass_forest_all3))
year_variance[6]<-round(bfv3all[which(bfv3all$grp=="yearfactor"),"vcov"],3)

bgv2<-as.data.frame(VarCorr(biomass_grass2))
year_variance[8]<-round(bgv2[which(bgv2$grp=="yearfactor"),"vcov"],3)
bgv3<-as.data.frame(VarCorr(biomass_grass3))
year_variance[9]<-round(bfv3[which(bgv3$grp=="yearfactor"),"vcov"],3)

afv2<-as.data.frame(VarCorr(abundance_forest2))
year_variance[11]<-round(afv2[which(afv2$grp=="yearfactor"),"vcov"],3)
afv3<-as.data.frame(VarCorr(abundance_forest3))
year_variance[12]<-round(afv3[which(afv3$grp=="yearfactor"),"vcov"],3)

afv2all<-as.data.frame(VarCorr(abundance_forest_all2))
year_variance[14]<-round(afv2all[which(afv2all$grp=="yearfactor"),"vcov"],3)
afv3all<-as.data.frame(VarCorr(abundance_forest_all3))
year_variance[15]<-round(afv3all[which(afv3all$grp=="yearfactor"),"vcov"],3)

agv2<-as.data.frame(VarCorr(abundance_grass2))
year_variance[17]<-round(agv2[which(agv2$grp=="yearfactor"),"vcov"],3)
agv3<-as.data.frame(VarCorr(abundance_grass3))
year_variance[18]<-round(agv3[which(agv3$grp=="yearfactor"),"vcov"],3)

sfv2<-as.data.frame(VarCorr(species_forest2))
year_variance[20]<-round(sfv2[which(sfv2$grp=="yearfactor"),"vcov"],3)
sfv3<-as.data.frame(VarCorr(species_forest3))
year_variance[21]<-round(sfv3[which(sfv3$grp=="yearfactor"),"vcov"],3)

sfv2all<-as.data.frame(VarCorr(species_forest_all2))
year_variance[23]<-round(sfv2all[which(sfv2all$grp=="yearfactor"),"vcov"],3)
sfv3all<-as.data.frame(VarCorr(species_forest_all3))
year_variance[24]<-round(sfv3all[which(sfv3all$grp=="yearfactor"),"vcov"],3)

sgv2<-as.data.frame(VarCorr(species_grass2))
year_variance[26]<-round(sgv2[which(sgv2$grp=="yearfactor"),"vcov"],3)
sgv3<-as.data.frame(VarCorr(species_grass3))
year_variance[27]<-round(sgv3[which(sgv3$grp=="yearfactor"),"vcov"],3)

#extract likelihood ratio and LR test p value for comparisons of model types 1 and 2.

bfyeareffect<-anova(biomass_forest2,biomass_forest1)
lik.ratio[2]<-round(bfyeareffect$Chisq[2],3)
randomterm_pval[2]<-round(bfyeareffect$'Pr(>Chisq)'[2],3)

bfyeareffectall<-anova(biomass_forest_all2,biomass_forest_all1)
lik.ratio[5]<-round(bfyeareffectall$Chisq[2],3)
randomterm_pval[5]<-round(bfyeareffectall$'Pr(>Chisq)'[2],3)

bgyeareffect<-anova(biomass_grass2,biomass_grass1)
lik.ratio[8]<-round(bgyeareffect$Chisq[2],3)
randomterm_pval[8]<-round(bgyeareffect$'Pr(>Chisq)'[2],3)

afyeareffect<-anova(abundance_forest2,abundance_forest1)
lik.ratio[11]<-round(afyeareffect$Chisq[2],3)
randomterm_pval[11]<-round(afyeareffect$'Pr(>Chisq)'[2],3)

afyeareffectall<-anova(abundance_forest_all2,abundance_forest_all1)
lik.ratio[14]<-round(afyeareffectall$Chisq[2],3)
randomterm_pval[14]<-round(afyeareffectall$'Pr(>Chisq)'[2],3)

agyeareffect<-anova(abundance_grass2,abundance_grass1)
lik.ratio[17]<-round(agyeareffect$Chisq[2],3)
randomterm_pval[17]<-round(agyeareffect$'Pr(>Chisq)'[2],3)

sfyeareffect<-anova(species_forest2,species_forest1)
lik.ratio[20]<-round(sfyeareffect$Chisq[2],3)
randomterm_pval[20]<-round(sfyeareffect$'Pr(>Chisq)'[2],3)

sfyeareffectall<-anova(species_forest_all2,species_forest_all1)
lik.ratio[23]<-round(sfyeareffectall$Chisq[2],3)
randomterm_pval[23]<-round(sfyeareffectall$'Pr(>Chisq)'[2],3)

sgyeareffect<-anova(species_grass2,species_grass1)
lik.ratio[26]<-round(sgyeareffect$Chisq[2],3)
randomterm_pval[26]<-round(sgyeareffect$'Pr(>Chisq)'[2],3)

#extract likelihood ratio and LR test p value for comparisons of model types 2 and 3.
bfyearslopevariance<-anova(biomass_forest2,biomass_forest1)
lik.ratio[3]<-round(bfyearslopevariance$Chisq[2],3)
randomterm_pval[3]<-round(bfyearslopevariance$'Pr(>Chisq)'[2],3)

bfyearslopevarianceall<-anova(biomass_forest_all2,biomass_forest_all1)
lik.ratio[6]<-round(bfyearslopevarianceall$Chisq[2],3)
randomterm_pval[6]<-round(bfyearslopevarianceall$'Pr(>Chisq)'[2],3)

bgyearslopevariance<-anova(biomass_grass2,biomass_grass1)
lik.ratio[9]<-round(bgyearslopevariance$Chisq[2],3)
randomterm_pval[9]<-round(bgyearslopevariance$'Pr(>Chisq)'[2],3)

afyearslopevariance<-anova(abundance_forest2,abundance_forest1)
lik.ratio[12]<-round(afyearslopevariance$Chisq[2],3)
randomterm_pval[12]<-round(afyearslopevariance$'Pr(>Chisq)'[2],3)

afyearslopevarianceall<-anova(abundance_forest_all2,abundance_forest_all1)
lik.ratio[15]<-round(afyearslopevarianceall$Chisq[2],3)
randomterm_pval[15]<-round(afyearslopevarianceall$'Pr(>Chisq)'[2],3)

agyearslopevariance<-anova(abundance_grass2,abundance_grass1)
lik.ratio[18]<-round(agyearslopevariance$Chisq[2],3)
randomterm_pval[18]<-round(agyearslopevariance$'Pr(>Chisq)'[2],3)

sfyearslopevariance<-anova(species_forest2,species_forest1)
lik.ratio[21]<-round(sfyearslopevariance$Chisq[2],3)
randomterm_pval[21]<-round(sfyearslopevariance$'Pr(>Chisq)'[2],3)

sfyearslopevarianceall<-anova(species_forest_all2,species_forest_all1)
lik.ratio[24]<-round(sfyearslopevarianceall$Chisq[2],3)
randomterm_pval[24]<-round(sfyearslopevarianceall$'Pr(>Chisq)'[2],3)

sgyearslopevariance<-anova(species_grass2,species_grass1)
lik.ratio[27]<-round(sgyearslopevariance$Chisq[2],3)
randomterm_pval[27]<-round(sgyearslopevariance$'Pr(>Chisq)'[2],3)

habitat<-rep(c(rep("forest",6),rep("grassland",3)),3)
response<-c(rep("biomass",9),rep("abundance",9),rep("species richness",9))
modeltype<-rep(c(1,2,3),9)
plots<-rep(c(rep(30,3),rep(140,3),rep(150,3)),3)

# output the summary table.
tableoutput<-cbind(habitat,response,plots,modeltype,coefsandse,propchange,pvalues,year_variance,lik.ratio,randomterm_pval)
write.table(tableoutput,"figures/modeltable.txt",sep="\t",col.names=T,row.names=F)