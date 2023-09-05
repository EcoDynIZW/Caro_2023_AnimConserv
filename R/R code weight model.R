###############################################################################
### Analyze coconut crab body mass as function of covariates ##################

### From Caro et al. (in press). Meta- and subpopulation estimation with 
###      disparate data: coconut crabs in the Western Indian Ocean.
###      Animal Conservation

#Note: all sourced/read files must be in the working directory, or else, paths 
#      must be adjusted as appropriate

###### The data object that is read in below contains the following information:
## $nsites: number of sites with weight measurements
## $n: number of weight measurements
## $site.index: index linking each weight measurement to one of 23 sites
## $weight: individual weight weasurements (response variable); individuals with 
##          multiple weight measures were reduced to a single random measurement
## $sex: sex of individual for each weight measurement
## $covs: data frame with covariate information for each of the 23 sites 
##        (for covariate description, see main text; Shambas = agriculture)


rm(list=ls())
library(readxl)
library(writexl)
library(nimble)
library(MCMCvis)
source('Nimble Weight.R')

wgt<-readRDS('weight data.rds')
#Note: individuals with multiple weight measures were reduced to a single
#      random measurement


##pull out covariates
covs.wgt<-wgt$covs #area level covariates

models<-c("Fishers", "Shambas","Hotel","Inhabited","Islands","ProtectedIsle")
nmodels<-length(models)

##compile constants for Nimble
nimConsts<-list(nsites=wgt$nsites,
                n=wgt$n,
                site.index=wgt$site.index) 

##initial values
inits<-function(){list(l.sd.res=runif(wgt$nsites, 0.1,0.2),
                       mu.s=runif(1,0,0.2),
                       sig.s=runif(1, 0.5,1),
                       sd.eps=runif(1, 0.5,1),
                       b0=runif(1, 0.5, 1),
                       b=runif(1, -0.5, 0.5),
                       b.sex=runif(1, 0.5, 1))}

##set parameters to monitor
params<-c('sd.eps', 'b0', 'b','b.sex','sd.res', 'epsilon', 'mu.s', 'sig.s')

##loop over covariates and fit models
for (jj in 1:nmodels){
  
  ##compile data for Nimble (square root transform weights)
  ##note: sex=1 are females, sex=0 are males
  nimDat<-list(y=sqrt(wgt$weight),
               sex=wgt$sex,
               COV=covs.wgt[,models[jj]])
  
  
  #(1) set up model
  model <- nimbleModel(Weight.all, constants = nimConsts, 
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
  summ<-MCMCsummary(samp, pg0=TRUE,probs = c(0.05,0.95))
  
  ##write out sample and summary
  saveRDS(samp, paste('MCMC_Weight_',models[jj] ,'.rds', sep=''))
  saveRDS(summ, paste('Summary_Weight_',models[jj] ,'.rds', sep=''))
  
} #end model loop


  
################################################################################
### get coefficients 
  
  est.mat<-as.data.frame(matrix(NA, nrow=length(models), ncol=8))
  colnames(est.mat)<-c('Model', 'Coefficient', 'SE', '5%', '95%','p>0', 
                       'Difference', 'SE(Difference)')
  
  est.mat$Model<-models
  
  for (jj in 1:length(models)){
    
    samp<-readRDS(paste('MCMC_Weight_',models[jj] ,'.rds', sep=''))
    
    summ<-readRDS(paste('Summary_Weight_',models[jj] ,'.rds', sep=''))
    est.mat[jj,c('Coefficient', 'SE', '5%', '95%','p>0')]<-round(unlist(summ['b',
                                                                           c("mean","sd","5%","95%","p>0")]),
                                                                 dig=2)
    ##add in avg difference in response
    b.post<-do.call(rbind,samp)
    b.post<-b.post[,c('b0','b')]
    #subtract predictor=0 from predictor=1, store posterior samples
    diff<-(b.post[,'b0']+b.post[,'b'])^2-(b.post[,'b0'])^2
    est.mat[jj,c('Difference', 'SE(Difference)')]<-round(c(mean(diff), sd(diff)), dig=2)
    
  }
  
  
  ##order by p>0
  
  pp0<-ifelse(est.mat$`p>0` >0.5, est.mat$`p>0`, 1-est.mat$`p>0`)
  ord2<-order(pp0,decreasing = TRUE)
  est.mat<-est.mat[ord2,]
  
  ###change to p!=0
  est.mat$`p!=0`<-pp0[ord2]
  est.mat.out<-est.mat[, c(1:5,9,7,8)]
  write_xlsx(est.mat.out, 'Weight coefficient estimates.xlsx')
  
  
  