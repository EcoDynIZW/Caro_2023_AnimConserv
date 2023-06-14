# Data and code for "Meta- and subpopulation estimation with disparate data: coconut crabs in the Western Indian Ocean.", Caro et al., Animal Conservation (in press)
---

All data and code necessary to (1) fit the integrated model combining capture-mark-recapture data and count data of coconut crabs collected across 29 sites in/near the Pemba archipelago, with and without covariates, (2) make posterior predictions of abundance across suitable habitat at all surveyed sites, (3) fit the GLMM to estimate the effect of covariates on catch per unit effort, and (4) fit the LMM to estimate the effect of covariates on coconut crab body mass. For details, see main text of publication.

## Description of the data and file structure

<ins> **This repository contains the following data files:** </ins> 

**(1) For the integrated model**  
*IntegratedModel_Inits.rds*: data to set initial values in integrated models with covariates  
*IntegratedModel_Nimble_Constants.rds*: data passed to Nimble as constants  
*IntegratedModel_Nimble_Data.rds*: data passed to Nimble as data

**(2) For the posterior prediction of abundance**  
*Posterior prediction site data.rds*: Site names and size (ha) of area with suitable habitat

**(3) For the catch per unit effort model**  
*CPUE data.rds*: constants, covariates and observations

**(4) For the weight model**  
*weight data.rds*: constants, covariates and observations

The content of each data file is described in detail at the beginning of the R script that makes use of it.


<ins>**In addition, the following R code files are available:**</ins>

*R code integrated model.R*: R code to load data, prep it for analysis and implement the integrated model, with and without covariates, and save model output

*R code posterior prediction abundance.R*: R code to load integrated model results, calculate abundance for each surveyed area and the Pemba archipelago overall, and compile all estimates in a table

*R code CPUE model.R*: R code to load data, prep it for analysis and implement the catch per unit effort model with covariates, save model output, and compile all covariate coefficients in a table

*R code weight model.R*: R code to load data, prep it for analysis and implement the weight model with covariates, save model output, and compile all covariate coefficients in a table


<ins>**Further, the following files with Nimble model code are provided**</ins>

*Nimble Integrated model.R*: integrated model, with and without covariate, as well as custom functions to vectorize dbern(), to ensure inclusion parameter psi is always <1, and to use correct detection parameter (from CMR model) in count model

*Nimble CPUE.R*: catch per unit effort model, with covariate

*Nimble weight.R*: weight model, with covariate


## Sharing/Access information

NA


## Code/Software

Data can be loaded into R and analyzed using the R package Nimble. Please see publication for version information and R scripts for additional packages used in processing data. 
