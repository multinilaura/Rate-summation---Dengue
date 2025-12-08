Zika Rate Summation - Preliminary analysis

This repository constains a preliminary analysis applying rate summation to Zika virus transmission by Aedes aegypti using previously published thermal performance cruves and modelling functions.

The analyses are based on functions published in "Mean daily temperatures predict the thermal limits of malaria transmission better than hourly rate summation" https://doi.org/10.1038/s41467-025-58612-w

The original code and data sources are available at:
	- Rate summation functions: https://github.com/mshocket/anopheles-rate-summation-release.
	- Thermal performance curves: https://github.com/mshocket/Thermal-Eco-MBD.


Software Requirements

Please use a version of R ≥ v4.4.1 - (available at https://www.r-project.org/) and JAGS ≥ v4.3.2 - (available at https://mcmc-jags.sourceforge.io/).

Instruction for Use

Follow the steps below to reproduce the analysis:

	- Access https://github.com/mshocket/anopheles-rate-summation-release and download the 'R-scipts folder'. This folder constains all R functions to run the scripts present in this repository.
	- Access https://github.com/mshocket/Thermal-Eco-MBD and download 'Trait Trajectories' folder, this folder contains the Meta Data for thermal performance curves fit for Mosquito-Borne Diseases Project. In this folder we will
	use trait data for AeaeDENV and AeaeZIKV.
	- Follow instructions for installing R and JAGS provided at the software links provided above.
	- Set working directory and place all downloaded folders and scripts into the same directory.
	- Run file Laura_TestAnalysis.R to generate TPC and apply rate summation for each life-history trait.
	- After all TCPs are generated, run file R0function_Laura.R to calculate and plot R0.
	
Notes

	- File paths may need to be adjusted depending on local system (Mac, PC, Linux, etc).
	- All scripts assume that required functions and trait data are correctly loaded into the working directory.
	
