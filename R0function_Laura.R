##### Load Packages
library(progress)
library(gridExtra)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot())




########## R0 function
####Tesla et al. 2018

# Creating a small constant to keep denominators from being zero.
#ec<-1/24

##################################
## Calculate the posterior distribution of each parameter and R0 vs. T
## Creating the function encoding the value of R0 as a function of the parameters
## Unlike in previous formulations, our b.samps fit includes bc (probability of infectiousness given exposure)

#myR0<-function(a, b, PDR, MDR, EFD, e2a, lf){
#mu = 1/(lf + ec)
#((a^2*b*(EFD*e2a*MDR/(mu)^2)*exp((-mu/(PDR+ec))))/(mu))^0.5
#}


##### Suitability function (Shocket et al. 2025)

# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
# ec.lf = 1/48

# R0 formulation where: 1) mosquito density (M) depends on lifetime fecundity (B); 2) gamma (y) is substituted for exp^(-1/(lf*PDR))expression
# R0eq = function(a, lf, B, y, bc, pEA, MDR) {
#	M = B * pEA * MDR * (lf+ec.lf)
#	R0 = (a^2 * bc * y * M * (lf+ec.lf) )^0.5
#	return(R0)
#}


########### Test R0 function (from Tesla paper)

#constant

ec<-1/24

myR0<-function(a, b, PDR, MDR, EFD, e2a, lf){
	mu = 1/(lf + ec)
	((a^2*b*(EFD*e2a*MDR/(mu)^2)*exp((-mu/(PDR+ec))))/(mu))^0.5
}


R0_constant <- myR0(predictions_bite_rate_dengue, predictions_vc_zika, predictions_eir_zika, predictions_mdr_dengue, 
					predictions_fecundity_dengue, predictions_pea_dengue, predictions_lf_zika)


R0_RS09 <- myR0(predictions_bite_rate_rs9_dengue, predictions_vc_rs9_zika, predictions_eir_rs9_zika, predictions_mdr_rs9_dengue, 
				predictions_fecundity_rs9_dengue, predictions_pea_rs9_dengue, predictions_lf_rs9_zika)


R0_RS12 <- myR0(predictions_bite_rate_rs12_dengue, predictions_vc_rs12_zika, predictions_eir_rs12_zika, predictions_mdr_rs12_dengue, 
				predictions_fecundity_rs12_dengue, predictions_pea_rs12_dengue, predictions_lf_rs12_zika)


#######Calculate quantiles for plotting

R0_constant_summary <- calcPostQuants(R0_constant, "R0_const", Temp.gradient)

R0_RS09_summary <- calcPostQuants(R0_RS09, "R0_dtr09", Temp.gradient)

R0_RS12_summary <- calcPostQuants(R0_RS12, "R0_dtr12", Temp.gradient)


# Calculate maximum predicted R0 value for each model

max(R0_constant_summary$median)
max(R0_RS09_summary$median)
max(R0_RS12_summary$median)



# Calculate % reduction in maximum R0 at Topt for each fluctuation model vs constant temperature model
1 - max(R0_RS09_summary$median)/max(R0_constant_summary$median) # 33%
1 - max(R0_RS12_summary$median)/max(R0_constant_summary$median) # 54%


# Calculate the max value of the median R0 from the constant temperature model
R0_scale <- max(R0_constant_summary$median)

# Get summary statistics of Tmin, Tmax, Topt, and Tbreadth 
params_R0_const <- extractDerivedTPC(R0_constant, "R0_const", Temp.gradient)

params_r0_rs09 <- extractDerivedTPC(R0_RS09, "r0_rs09", Temp.gradient)
params_r0_rs12 <- extractDerivedTPC(R0_RS12, "r0_rs12", Temp.gradient)

# #### Comparison: combine rate summation and constant
rseffect_r0_predictions <- bind_rows(R0_constant_summary, R0_RS09_summary, R0_RS12_summary)
rseffect_r0_params <- bind_rows(params_R0_const, params_r0_rs09, params_r0_rs12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))



#### Comparison: combine all treatments constant and rate summation

R0_plot <- rseffect_r0_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_rssuit09, c_rssuit12), labels = c("1: Constant", "2: Rate summation DTR 9", "2: Rate summation DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_rssuit09, ct_rssuit12), labels = c("1: Constant", "2: Rate summation DTR 9", "2: Rate summation DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1), labels = c("1: Constant", "2: Rate summation DTR 9", "2: Rate summation DTR 12")) +
	ylab("Suitability for transmission (R0)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank())

rseffect_r0_params_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rseffect_r0_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_constant, c_rssuit09, c_rssuit12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank())

plot_R0_combined <- R0_plot / rseffect_r0_params_plot + plot_layout(heights = c(3, 0.75))







#### Rate summation on R0_constant

### assuming dtr9

R0_const_dtr9 <- RSCalcTempGrad(R0_constant, LPtemps_dtr9, Temp.gradient)

plot(as.numeric(R0_const_dtr12[1,]) ~ Temp.gradient)

for(i in 1:20){
	lines(as.numeric(R0_const_dtr12[i,]) ~ Temp.gradient)
}

### assuming dtr12

R0_const_dtr12 <- RSCalcTempGrad(R0_constant, LPtemps_dtr12, Temp.gradient)


#######Calculate quantiles for plotting

R0_constant_summary <- calcPostQuants(R0_constant, "R0_const", Temp.gradient)

R0_const_dtr9_summary <- calcPostQuants(R0_const_dtr9, "R0_const_dtr09", Temp.gradient)

R0_const_dtr12_summary <- calcPostQuants(R0_const_dtr12, "R0_const_dtr12", Temp.gradient)

# Calculate maximum predicted R0 value for each model

max(R0_constant_summary$median)
max(R0_const_dtr9_summary$median)
max(R0_const_dtr12_summary$median)


# Calculate the max value of the median R0 from the constant temperature model
R0_scale <- max(R0_constant_summary$median)

# Get summary statistics of Tmin, Tmax, Topt, and Tbreadth 
params_R0_const <- extractDerivedTPC(R0_constant, "R0_const", Temp.gradient)

params_r0_const_dtr9 <- extractDerivedTPC(R0_const_dtr9, "r0_const_dtr9", Temp.gradient)
params_r0_const_dtr12 <- extractDerivedTPC(R0_const_dtr12, "r0_const_dtr12", Temp.gradient)

# #### Comparison: combine rate summation and constant
rseffect_r0_predictions_model2 <- bind_rows(R0_constant_summary, R0_const_dtr9_summary, R0_const_dtr12_summary)
rseffect_r0_params_model2 <- bind_rows(params_R0_const, params_r0_const_dtr9, params_r0_const_dtr12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))


#### Comparison: combine all treatments constant and rate summation

R0_plot_2 <- rseffect_r0_predictions_model2 %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_rssuit09, c_rssuit12), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_rssuit09, ct_rssuit12), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	ylab("Suitability for transmission (R0)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank())

rseffect_r0_params_plot_2 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rseffect_r0_params_model2, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_constant, c_rssuit09, c_rssuit12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank())

plot_R0_combined_2 <- R0_plot_2 / rseffect_r0_params_plot_2 + plot_layout(heights = c(3, 0.75))


R0_const_dtr12_summary

R0_plot_3 <- R0_const_dtr12_summary %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_rssuit09, c_rssuit12), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_rssuit09, ct_rssuit12), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	ylab("Suitability for transmission (R0)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank())


R0_plot_4 <- R0_const_dtr9_summary %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_rssuit09, c_rssuit12), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_rssuit09, ct_rssuit12), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1), labels = c("1: Constant", "2: Rate summation on R0 - DTR 9", "2: Rate summation on R0 - DTR 12")) +
	ylab("Suitability for transmission (R0)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank())

save.image()
