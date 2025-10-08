
####Set working directory
#setwd()

###Load packages

##### Load JAGS libraries
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('progress')
library("tidyverse")


##### Load functions from Shocket et al. 202
source("R-scripts/working-versions-code/00_RSProjectFunctions.R")

##### Load TPC fittings data from constant temperatures calculated by Mordecai 2019, and transform t() data
#bite rate (a)
dengue.bite.rate.a <- read.csv("data-mordecai/AeaeDENV.a.csv", header = F)
dengue.bite.rate.a.trans <- t(dengue.bite.rate.a)

#fecundity (EFD)
denv.fecundity.efd <- read.csv("data-mordecai/AeaeDENV.EFD.csv", header = F)
denv.fecundity.efd.t <- t(denv.fecundity.efd)

#egg to adult development rate (MDR)
denv.dev.mdr <- read.csv("data-mordecai/AeaeDENV.MDR.csv", header = F)
denv.dev.mdr.t <- t(denv.dev.mdr)

#survival probability (pEA)
denv.surv.pea <- read.csv("data-mordecai/AeaeDENV.pEA.csv", header = F)
denv.surv.pea.t <- t(denv.surv.pea)

#average mosquito lifespan (lf)
zikv.lf <- read.csv("data-mordecai/AeaeZIKV.lf.csv", header = F)
zikv.lf.t <- t(zikv.lf)

#extrinsic incubation rate (PDR - Parasite development rate)
zikv.eir.pdr <- read.csv("data-mordecai/AeaeZIKV.PDR.csv", header = F)
zikv.eir.pdr.t <- t(zikv.eir.pdr)

#vector competence (b) (proportion of infected that successfully transmit)
zikv.vc.bc <- read.csv("data-mordecai/AeaeZIKV.b.csv", header = F)
zikv.vc.bc.t <- t(zikv.vc.bc)

#### Perform rate summation calculation on the TPC

# Temperature gradient that matches TPC predictions
Temp.gradient <- seq(5, 45, 0.1)

# Generate hourly temperature sequences across the temperature gradient
LPtemps_dtr9 <- LoganPartonCalc(dtr = 9, Temp.gradient)
LPtemps_dtr12 <- LoganPartonCalc(dtr = 12, Temp.gradient)

# Set negative temperature values to 5, since TPCs stop at 5 on the low end
LPtemps_dtr9[LPtemps_dtr9 < 5 ] <- 5
LPtemps_dtr12[LPtemps_dtr12 < 5 ] <- 5

# Set temperature values > 45 to 45, since TPCs stop at 45 on the high end
LPtemps_dtr9[LPtemps_dtr9 > 45 ] <- 45
LPtemps_dtr12[LPtemps_dtr12 > 45 ] <- 45


######################################################################################################################################
#### Apply rate summation - bite rate (dengue)

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_bite_rate_dengue <- as.data.frame(dengue.bite.rate.a.trans, col_names = Temp.gradient)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_bite_rate_rs9_dengue <- RSCalcTempGrad(predictions_bite_rate_dengue, LPtemps_dtr9, Temp.gradient)

predictions_bite_rate_rs12_dengue <- RSCalcTempGrad(predictions_bite_rate_dengue, LPtemps_dtr12, Temp.gradient)


#Process constant temp calculation output from TPC predictions
predictions_bite_rate_dengue_summary <- calcPostQuants(predictions_bite_rate_dengue, "bite_rate_dengue_const", Temp.gradient)
params_bite_rate_dengue_summary <- extractDerivedTPC(predictions_bite_rate_dengue, "bite_rate_dengue_const", Temp.gradient)

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_bite_rate_rs9_dengue_summary <- calcPostQuants(predictions_bite_rate_rs9_dengue, "bite_rate_rs09_dengue", Temp.gradient)
predictions_bite_rate_rs12_dengue_summary <- calcPostQuants(predictions_bite_rate_rs12_dengue, "bite_rate_rs12_dengue", Temp.gradient)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_bite_rate_rs9_summary <- extractDerivedTPC(predictions_bite_rate_rs9_dengue, "bite_rate_rs09_dengue", Temp.gradient)
params_bite_rate_rs12_summary <- extractDerivedTPC(predictions_bite_rate_rs12_dengue, "bite_rate_rs12_dengue", Temp.gradient)


### combine empirically measured constant and rate summation fluctuation treatments
bite_rate_dengue_predictions <- bind_rows(predictions_bite_rate_dengue_summary, predictions_bite_rate_rs9_dengue_summary, predictions_bite_rate_rs12_dengue_summary)

bite_rate_dengue_plot <- bite_rate_dengue_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Bite~rate~(day^-1)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2.5, 1, 1, 2.5), "lines"), axis.title.y = element_text(vjust = 1)) +
	coord_cartesian(clip = "off")

bite_rate_dengue_plot

######################################################################################################################################
#### Apply rate summation - fecundity (dengue)

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_fecundity_dengue <- as.data.frame(denv.fecundity.efd.t, col_names = Temp.gradient)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_fecundity_rs9_dengue <- RSCalcTempGrad(predictions_fecundity_dengue, LPtemps_dtr9, Temp.gradient)

predictions_fecundity_rs12_dengue <- RSCalcTempGrad(predictions_fecundity_dengue, LPtemps_dtr12, Temp.gradient)

#Process constant temp calculation output from TPC predictions
predictions_fecundity_dengue_summary <- calcPostQuants(predictions_fecundity_dengue, "fecundity_dengue_const", Temp.gradient)
params_fecundity_dengue_summary <- extractDerivedTPC(predictions_fecundity_dengue, "fecundity_dengue_const", Temp.gradient)

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_fecundity_rs9_dengue_summary <- calcPostQuants(predictions_fecundity_rs9_dengue, "fecundity_rs09_dengue", Temp.gradient)
predictions_fecundity_rs12_dengue_summary <- calcPostQuants(predictions_fecundity_rs12_dengue, "fecundity_rs12_dengue", Temp.gradient)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_fecundity_rs9_summary <- extractDerivedTPC(predictions_fecundity_rs9_dengue, "fecundity_rs09_dengue", Temp.gradient)
params_fecundity_rs12_summary <- extractDerivedTPC(predictions_fecundity_rs12_dengue, "fecundity_rs12_dengue", Temp.gradient)

### combine empirically measured constant and rate summation fluctuation treatments
fecundity_dengue_predictions <- bind_rows(predictions_fecundity_dengue_summary, predictions_fecundity_rs9_dengue_summary, predictions_fecundity_rs12_dengue_summary)

fecundity_dengue_plot <- fecundity_dengue_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Fecundity")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2, 1, 1, 2), "lines"), axis.title.y = element_text(vjust = 1)) +
	coord_cartesian(clip = "off")

fecundity_dengue_plot

######################################################################################################################
#### Apply rate summation - Egg to adult survival (dengue)

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_mdr_dengue <- as.data.frame(denv.dev.mdr.t, col_names = Temp.gradient)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_mdr_rs9_dengue <- RSCalcTempGrad(predictions_mdr_dengue, LPtemps_dtr9, Temp.gradient)

predictions_mdr_rs12_dengue <- RSCalcTempGrad(predictions_mdr_dengue, LPtemps_dtr12, Temp.gradient)

#Process constant temp calculation output from TPC predictions
predictions_mdr_dengue_summary <- calcPostQuants(predictions_mdr_dengue, "mdr_dengue_const", Temp.gradient)
params_mdr_dengue_summary <- extractDerivedTPC(predictions_mdr_dengue, "mdr_dengue_const", Temp.gradient)

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_mdr_rs9_dengue_summary <- calcPostQuants(predictions_mdr_rs9_dengue, "mdr_rs09_dengue", Temp.gradient)
predictions_mdr_rs12_dengue_summary <- calcPostQuants(predictions_mdr_rs12_dengue, "mdr_rs12_dengue", Temp.gradient)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_mdr_rs9_summary <- extractDerivedTPC(predictions_mdr_rs9_dengue, "mdr_rs09_dengue", Temp.gradient)
params_mdr_rs12_summary <- extractDerivedTPC(predictions_mdr_rs12_dengue, "mdr_rs12_dengue", Temp.gradient)

### combine empirically measured constant and rate summation fluctuation treatments
mdr_dengue_predictions <- bind_rows(predictions_mdr_dengue_summary, predictions_mdr_rs9_dengue_summary, predictions_mdr_rs12_dengue_summary)

mdr_dengue_plot <- mdr_dengue_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Egg~to~adult~development~rate")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2, 1, 1, 2), "lines"), axis.title.y = element_text(vjust = 1)) +
	coord_cartesian(clip = "off")

mdr_dengue_plot

######################################################################################################################################
#### Apply rate summation - Survival probability (dengue)

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_pea_dengue <- as.data.frame(denv.surv.pea.t, col_names = Temp.gradient)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_pea_rs9_dengue <- RSCalcTempGrad(predictions_pea_dengue, LPtemps_dtr9, Temp.gradient)

predictions_pea_rs12_dengue <- RSCalcTempGrad(predictions_pea_dengue, LPtemps_dtr12, Temp.gradient)

#Process constant temp calculation output from TPC predictions
predictions_pea_dengue_summary <- calcPostQuants(predictions_pea_dengue, "pea_dengue_const", Temp.gradient)
params_pea_summary <- extractDerivedTPC(predictions_pea_dengue, "pea_dengue_const", Temp.gradient)

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_pea_rs9_dengue_summary <- calcPostQuants(predictions_pea_rs9_dengue, "pea_rs09_dengue", Temp.gradient)
predictions_pea_rs12_dengue_summary <- calcPostQuants(predictions_pea_rs12_dengue, "pea_rs12_dengue", Temp.gradient)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_pea_rs9_summary <- extractDerivedTPC(predictions_pea_rs9_dengue, "pea_rs09_dengue", Temp.gradient)
params_pea_rs12_summary <- extractDerivedTPC(predictions_pea_rs12_dengue, "pea_rs12_dengue", Temp.gradient)

### combine empirically measured constant and rate summation fluctuation treatments
pea_dengue_predictions <- bind_rows(predictions_pea_dengue_summary, predictions_pea_rs9_dengue_summary, predictions_pea_rs12_dengue_summary)

pea_dengue_plot <- pea_dengue_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Survival~probability~(pEA)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2, 1, 1, 2), "lines"), axis.title.y = element_text(vjust = 1)) +
	coord_cartesian(clip = "off")

pea_dengue_plot


save.image()


######################################################################################################################################
#### Apply rate summation - Lifespan (Zika)

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_lf_zika <- as.data.frame(zikv.lf.t, col_names = Temp.gradient)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_lf_rs9_zika <- RSCalcTempGrad(predictions_lf_zika, LPtemps_dtr9, Temp.gradient)

predictions_lf_rs12_zika <- RSCalcTempGrad(predictions_lf_zika, LPtemps_dtr12, Temp.gradient)

#Process constant temp calculation output from TPC predictions
predictions_lf_zika_summary <- calcPostQuants(predictions_lf_zika, "zika_lf_const", Temp.gradient)
params_lf_summary <- extractDerivedTPC(predictions_lf_zika, "zika_lf_const", Temp.gradient)

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_lf_rs9_zika_summary <- calcPostQuants(predictions_lf_rs9_zika, "zika_lf_rs09", Temp.gradient)
predictions_lf_rs12_zika_summary <- calcPostQuants(predictions_lf_rs12_zika, "zika_lf_rs12", Temp.gradient)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_lf_rs9_summary <- extractDerivedTPC(predictions_lf_rs9_zika, "zika_lf_rs09", Temp.gradient)
params_lf_rs12_summary <- extractDerivedTPC(predictions_lf_rs12_zika, "zika_lf_rs12", Temp.gradient)

### combine empirically measured constant and rate summation fluctuation treatments
lf_zika_predictions <- bind_rows(predictions_lf_zika_summary, predictions_lf_rs9_zika_summary, predictions_lf_rs12_zika_summary)

lf_zika_plot <- lf_zika_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Lifespan~(days)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2, 1, 1, 2), "lines"), axis.title.y = element_text(vjust = 1)) +
	coord_cartesian(clip = "off")

lf_zika_plot


######################################################################################################################################
#### Apply rate summation - Extrinsic incubation rate (EIR)

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_eir_zika <- as.data.frame(zikv.eir.pdr.t, col_names = Temp.gradient)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_eir_rs9_zika <- RSCalcTempGrad(predictions_eir_zika, LPtemps_dtr9, Temp.gradient)

predictions_eir_rs12_zika <- RSCalcTempGrad(predictions_eir_zika, LPtemps_dtr12, Temp.gradient)

#Process constant temp calculation output from TPC predictions
predictions_eir_zika_summary <- calcPostQuants(predictions_eir_zika, "zika_eir_const", Temp.gradient)
params_eir_summary <- extractDerivedTPC(predictions_eir_zika, "zika_eir_const", Temp.gradient)

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_eir_rs9_zika_summary <- calcPostQuants(predictions_eir_rs9_zika, "zika_eir_rs09", Temp.gradient)
predictions_eir_rs12_zika_summary <- calcPostQuants(predictions_eir_rs12_zika, "zika_eir_rs12", Temp.gradient)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_eir_rs9_summary <- extractDerivedTPC(predictions_eir_rs9_zika, "zika_eir_rs09", Temp.gradient)
params_eir_rs12_summary <- extractDerivedTPC(predictions_eir_rs12_zika, "zika_eir_rs12", Temp.gradient)

### combine empirically measured constant and rate summation fluctuation treatments
eir_zika_predictions <- bind_rows(predictions_eir_zika_summary, predictions_eir_rs9_zika_summary, predictions_eir_rs12_zika_summary)

eir_zika_plot <- eir_zika_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Extrinsic~Incubation~rate")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2, 1, 1, 2), "lines"), axis.title.y = element_text(vjust = 1)) +
	coord_cartesian(clip = "off")

eir_zika_plot

######################################################################################################################################
#### Apply rate summation - Vector competence (bc) - Proportion of infected that successfully transmit

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames
predictions_vc_zika <- as.data.frame(zikv.vc.bc.t, col_names = Temp.gradient)

# Apply rate summation for DTR = 9 and 12 - these take about 1.5-2 minutes each to run
predictions_vc_rs9_zika <- RSCalcTempGrad(predictions_vc_zika, LPtemps_dtr9, Temp.gradient)

predictions_vc_rs12_zika <- RSCalcTempGrad(predictions_vc_zika, LPtemps_dtr12, Temp.gradient)

#Process constant temp calculation output from TPC predictions
predictions_vc_zika_summary <- calcPostQuants(predictions_vc_zika, "zika_vc_const", Temp.gradient)
params_vc_summary <- extractDerivedTPC(predictions_vc_zika, "zika_vc_const", Temp.gradient)

# Process the rate summation calculation output - generate a data frame of mean, median, lower CI, and upper CI across the temperature gradient 
predictions_vc_rs9_zika_summary <- calcPostQuants(predictions_vc_rs9_zika, "zika_vc_rs09", Temp.gradient)
predictions_vc_rs12_zika_summary <- calcPostQuants(predictions_vc_rs12_zika, "zika_vc_rs12", Temp.gradient)

# Process the rate summation calculation output - get summary statistics of Tmin, Tmax, Topt, and Tbreadth
params_vc_rs9_summary <- extractDerivedTPC(predictions_vc_rs9_zika, "zika_vc_rs09", Temp.gradient)
params_vc_rs12_summary <- extractDerivedTPC(predictions_vc_rs12_zika, "zika_vc_rs12", Temp.gradient)

### combine empirically measured constant and rate summation fluctuation treatments
vc_zika_predictions <- bind_rows(predictions_vc_zika_summary, predictions_vc_rs9_zika_summary, predictions_vc_rs12_zika_summary)

vc_zika_plot <- vc_zika_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), linewidth = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("Constant", "Rate summation - DTR9", "Rate summation - DTR12")) +
	ylab(parse(text = "Vector~competence")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2, 1, 1, 2), "lines"), axis.title.y = element_text(vjust = 1)) +
	coord_cartesian(clip = "off")

vc_zika_plot





save.image()





