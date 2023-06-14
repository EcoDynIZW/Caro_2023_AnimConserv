################################################################################
#### Analyze coconut crab counts as function of covariates #####################

### From Caro et al. (in press). Meta- and subpopulation estimation with 
###      disparate data: coconut crabs in the Western Indian Ocean.
###      Animal Conservation

## Note that all sourced/read files need to be in the working directory, or else,
## their paths must be adjusted accordingly

###### The data object that is read in below contains the following information:
## $nsites: number of sites surveyed 
## $nsurveys: total number of surveys
## $site.index: index linking each survey to one of 29 sites
## $obs: number of crabs caught in each survey (response variable) 
## $effort: effort (in minutes spent searching/100) for each survey
## $covs: data frame with covariate information for each of the 29 sites 
##        (for covariate description, see main text; Shambas = agriculture)

rm(list=ls())
library(readxl)
library(writexl)
library(nimble)
library(MCMCvis)

source('Nimble CPUE.R')


##get prepped CPUE data
dat<-readRDS('CPUE data.rds')

#note: obs are aggregate counts over all visits within a survey

## pull out some pieces
nsurveys<-dat$nsurveys #number of surveys (here, data points)
nsites<-dat$nsites #number of areas
covs<-dat$covs #area level covariates
site.index<-dat$site.index #index linking survey to area

obs<-dat$obs #number of crabs caught per survey
effort<-dat$effort #total sampling effort per survey


models<-c("Fishers", "Shambas","Hotel","Inhabited","Islands","ProtectedIsle")
nmodels<-length(models)


##compile constants for Nimble
nimConsts<-list(nsites=nsites,
                nsurveys=nsurveys,
                #nvisits=nvisits,
                site.index=site.index) #,

#initial values
inits<-function(){list(r=runif(1, 1,2),
                       sd.eps=runif(1, 0.5,1),
                       b0=runif(1, 0.5, 1),
                       b=runif(1, -0.5, 0.5))}

##set parameters to monitor
params<-c('sd.eps', 'b0', 'b','r', 'epsilon')


##loop over all covariates and run models
for (jj in 1:nmodels){
  
  ##compile data for Nimble
  nimDat<-list(y=obs,
               effort=effort,
               COV=covs[,models[jj]])
  
  
  #(1) set up model
  model <- nimbleModel(Cpue.all, constants = nimConsts, 
                       data=nimDat, inits=inits(), check = FALSE)
  
  #(2) Compile model in c++
  cmodel <- compileNimble(model)       
  
  # (3) Configure MCMC - on an uncompiled model 
  conf.mcmc<-configureMCMC(model, monitors = params, thin=10)#

  # (4) Build the MCMC sampler based on configurations
  mcmc <- buildMCMC(conf.mcmc)
  
  # (5) Compile sampler in c++ together with compiled model
  cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)
  
  # (6) Run 
  samp <- runMCMC(cmcmc, niter = 75000, nburnin = 50000, nchains=3, 
                  inits = inits) 
  
  ## calculate summary and model WAIC
  summ<-MCMCsummary(samp, pg0=TRUE)
  
  ##write out posterior samples
  saveRDS(samp, paste('MCMC_CPUE_',models[[jj]] ,'.rds', sep=''))
  
} #end model loop


################################################################################
### get coefficients of all models and compile in results table ################

est.mat<-as.data.frame(matrix(NA, nrow=length(models), ncol=8))
colnames(est.mat)<-c('Model', 'Coefficient', 'SE', '5%', '95%','p>0', 
                     'Difference', 'SE(Difference)')

est.mat$Model<-models

for (jj in 1:length(models)){
  
  samp<-readRDS(paste('MCMC_CPUE_',models[jj] ,'.rds', sep=''))
  summ<-MCMCsummary(samp, probs = c(0.05, 0.95), pg0 = TRUE, round=2)
  
  est.mat[jj,c('Coefficient', 'SE', '5%', '95%','p>0')]<-unlist(summ['b',
                                                                         c("mean","sd","5%","95%","p>0")])
  ##add in avg difference in response
  b.post<-do.call(rbind,samp)
  b.post<-b.post[,c('b0','b')]
  #subtract predictor=0 from predictor=1, store posterior samples
  diff<-exp(b.post[,'b0']+b.post[,'b'])-exp(b.post[,'b0'])
  est.mat[jj,c('Difference', 'SE(Difference)')]<-round(c(mean(diff), sd(diff)), dig=2)
  
}


##order by p>0

pp0<-ifelse(est.mat$`p>0` >0.5, est.mat$`p>0`, 1-est.mat$`p>0`)
ord2<-order(pp0,decreasing = TRUE)
est.mat<-est.mat[ord2,]

###change to p!=0
est.mat$`p!=0`<-pp0[ord2]
est.mat.out<-est.mat[, c(1:5,9,7,8)]
write_xlsx(est.mat.out, 'CPUE Coefficient estimates.xlsx')

