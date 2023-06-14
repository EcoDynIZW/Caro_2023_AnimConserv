################################################################################
### Fit integrated model with and without covariates ###########################

### From Caro et al. (in press). Meta- and subpopulation estimation with 
###      disparate data: coconut crabs in the Western Indian Ocean.
###      Animal Conservation

## All data and functions must be in the working directory, or else, paths in the
## code must be adjusted

###### The code reads in 3 data objects, which contain the following information:

#### 1. nimConsts
## $M: for CMR surveys, size of the augmented data set 
## $K.cmr: for CMR surveys, the number of sampling occasions/visits
## $nsurv.cmr: number of CMR surveys
## $site: index linking each CMR survey with the survey area; different observers
##        on Misali and Chumbe are treated as different 'sites' in this index; for 
##        detection submodel
## $nsites.cmr: number of site-observer combinations with CMR surveys
## $site.index.cmr: index linking CMR surveys to actual survey areas (regardless
##                  of observer identity) - for ecological sub-model
## $area.cmr: size of sampled area (ha) for CMR surveys
## $nsites: number of survey areas
## $nsurv.count: number of count surveys
## $K.count: for count surveys, the number of sampling occasions/visits
## $site.index.count: index linking count surveys to survey areas
## $area.count: size of sampled area (ha) for count surveys
## $cmr.index: index showing which CMR site a count survey was conducted at (-999
##             if count survey happened at a site without any CMR surveys); for 
##             assigning appropriate CMR based detection parameters to count surveys

#### 2. nimDatAll
## $y.cmr: array with CMR data; max. data augmentation x number of CMR surveys x
##         max. number of sampling occasions
## $effort.cmr: sampling effort (in minutes searching/100) for all CMR surveys; 
##              number of CMR surveys x max. number of sampling occasions
## $y.count: matrix with count data; number of count surveys x max. number of 
##           sampling occasions
## $effort.count: sampling effort (in minutes searching/100) for all count surveys; 
##              number of surveys x max. number of sampling occasions
## $cov.cmr: area-level covariates for all CMR surveys (surveys at the same site
##          have the same covariate value; for description of covariates, see main
##          text; shambas = agriculture; non-covariate columns can be ignored)
## $cov.count: area-level covariates for all count surveys
## $z: binary indicator for CMR data marking if an individual was observed (z=1)
##     or augmented (z=0)
## $Nest: estimated abundance for CMR surveys from ML model; used to set initial
##        values in the null model

#### 3. ins
## List of 6, with parameter estimates from a model treating all surveys as count
## surveys; 1 element for each covariate; used to set initial values for integrated
## models with covariates. Each element is a data frame with the following data:
##  columns: posterior summary statistics of model parameters
##  rows: parameters
##        - b and b0: covariate coefficient and intercept of density model
##        - epsilon: area level random effect values of density model
##        - lp.b, p.int: coefficients of effort, intercepts (for 3 observers) 
##                       for detection model
##        - r: negative binomial overdispersion parameter
##        - sd.eps: standard deviation of area random effect

rm(list=ls())

library(nimble)
library(MCMCvis)
library(writexl)
library(readxl)

#load model code
source('Nimble Integrated model.R')

##function to get posterior mode
get.mode<-function(x){
  xx<-density(x)
  out<-xx$x[xx$y == max(xx$y)]
  return(out)
}

##read in some objects that hold data, constants and priors for model

##constants, for covariate and null model(s)
nimConsts<-readRDS('IntegratedModel_Nimble_Constants.rds')

##data (observations, effort, covariates), needs to be further subset to run models
nimDatAll<-readRDS('IntegratedModel_Nimble_Data.rds')

##info to set initial values for covariate models
ins<-readRDS('IntegratedModel_Inits.rds')

##pull out some pieces
covs.cmr<-nimDatAll$cov.cmr
M.site<-nimConsts$M
nsites.cmr<-nimConsts$nsites.cmr
covs.count<-nimDatAll$cov.count

################################################################################
### Covariate models ###########################################################

###covariates to model
models<-c("Fishers", "Shambas","Hotel","Inhabited","Islands","ProtectedIsle")
nmodels<-length(models)


##create initial values for z, for CMR part of model
M.max<-nrow(nimDatAll$z)
npop<-ncol(nimDatAll$z)
z<-nimDatAll$z
z.in<-matrix(NA, M.max, npop)
z.in[is.na(z)]<-0


##initial values for N per survey as max obs +1 (for N-mixture part)
N.in<-apply(nimDatAll$y.count,1,function(x)max(x, na.rm=TRUE)+1)


##loop over the 6 covariates
##note: this will run for several hours

for (jj in 1:nmodels){

  ##set initial values based on simplified model estimates of abundance
  ##and covariate coefficient
  ests<-ins[[jj]]$mean

  nest<-exp(ests[pmatch('b0',rownames(ins[[jj]]))]+ests[pmatch('b',rownames(ins[[jj]]))]*covs.cmr[,models[jj]]+
    log(covs.cmr$`Sampled area`)+ests[grep('epsilon',rownames(ins[[jj]]))][covs.cmr$Area])
  
  inits<-function(){list(z=z.in,
                         p.int=runif(nsites.cmr, 0, 0.1),
                         eps.site=ests[grep('epsilon',rownames(ins[[jj]]))],
                         lp.b=runif(nsites.cmr,0.5, 0.75),
                         N.count=N.in, 
                         r=runif(1, 1,2),
                         sd.eps=ests[pmatch('sd.eps',rownames(ins[[jj]]))],
                         sd.psi=runif(1, 0.01,0.02),
                         psi=runif(npop, nest/M.site-0.05, nest/M.site+0.05),
                         b0=ests[pmatch('b0',rownames(ins[[jj]]))],
                         b=ests[pmatch('b',rownames(ins[[jj]]))])} #
  
  
nimDat<-list(y.cmr=nimDatAll$y.cmr, 
              effort.cmr=nimDatAll$effort.cmr,
              y.count = nimDatAll$y.count,
              effort.count=nimDatAll$effort.count,
             cov.cmr=covs.cmr[,models[jj]],
             cov.count=covs.count[,models[jj]],
             z=z)


params<-c('sd.eps','sd.psi',  'b0', 'b','r',
           'N.cmr', 'lp.int', 'lp.b',
           'psi','eps.site')

#(1) set up model
model <- nimbleModel(Integrated.model, constants = nimConsts, 
                     data=nimDat, inits=inits(), check = FALSE)

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model 
conf.mcmc<-configureMCMC(model, monitors = params, thin=10)#
##use block sampling for b and b0
conf.mcmc$removeSamplers(c('b', 'b0'))
conf.mcmc$addSampler(target=c('b', 'b0'), type = 'RW_block')

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run 
# To speed up computation, run in 3 separate instances of R, save each chain
# separately; adjust name under saveRDS (_1, _2, _3) so as to not overwrite 
# chains
# Note that at these settings, this may take several hours to run
system.time(
  (samp <- runMCMC(cmcmc, niter = 100000, nburnin = 50000, nchains=1, 
                   inits = inits) )
)
saveRDS(samp, paste('MCMC_',models[jj] ,'_1','.rds', sep=''))

}


#################################################################################
### write out table with coefficient estimates and avg differences in response ##

est.mat<-as.data.frame(matrix(NA, nrow=length(models), ncol=8))
colnames(est.mat)<-c('Model', 'Coefficient', 'SE', '5%', '95%','p>0', 
                     'Difference', 'SE(Difference)')

est.mat$Model<-models

for (jj in 1:length(models)){
  
  ###compile 3 chains
  samp<-list()
  for (i in 1:3){
    samp[[i]]<-readRDS(paste('MCMC_',models[jj] ,'_',i,'.rds', sep=''))
  }
  summ<-MCMCsummary(samp, pg0=TRUE, probs = c(0.05, 0.95),round=2)
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
write_xlsx(est.mat.out, 'Coefficient estimates.xlsx')


################################################################################
### Null model #################################################################
################################################################################

nimDat0<-list(y.cmr=nimDatAll$y.cmr, 
             effort.cmr=nimDatAll$effort.cmr,
             y.count = nimDatAll$y.count,
             effort.count=nimDatAll$effort.count,
             z=z)

inits0<-function(){list(z=z.in,
                         p.int=runif(nsites.cmr, 0, 0.1),
                         #mu.lpb=runif(1, 0.2, 0.7),
                         #sig.lpb=runif(1, 0.5,1),
                        eps.site=runif(nimConsts$nsites, -0.01, 0.01),
                        lp.b=runif(nsites.cmr,0.5, 0.75),
                       ##count model
                       N.count=N.in, 
                       r=runif(1, 1,2),
                       ##shared
                       sd.eps=runif(1, 0.01,0.02),
                       sd.psi=runif(1, 0.01,0.02),
                       psi=runif(npop, nimDatAll$Nest/M.site-0.1, nimDatAll$Nest/M.site+0.1),
                       #epsilon=eps,
                       b0=runif(1, 0.8, 1.2))} #b=runif(1, -0.5, 0.5)


params0<-c('sd.eps','sd.psi',  'b0', 'r',
           'N.cmr', 'lp.int', 'lp.b',
           'psi','eps.site')


#(1) set up model
model <- nimbleModel(Integrated.model0, constants = nimConsts, 
                     data=nimDat0, inits=inits0(), check = FALSE)

#(2) Compile model in c++
cmodel <- compileNimble(model)       

# (3) Configure MCMC - on an uncompiled model 
conf.mcmc<-configureMCMC(model, monitors = params0, thin=10)#

# (4) Build the MCMC sampler based on configurations
mcmc <- buildMCMC(conf.mcmc)

# (5) Compile sampler in c++ together with compiled model
cmcmc <- compileNimble(mcmc, project = cmodel, resetFunctions = TRUE)

# (6) Run 
# Faster than covariate models, so run can be run in parallel
# Note that at these settings, this may take several hours
system.time(
(samp <- runMCMC(cmcmc,  niter = 100000, nburnin = 50000, nchains=3, 
                inits = inits0) )
)
summ<-MCMCsummary(samp, pg0 = TRUE, func=get.mode)
MCMCtrace(samp, filename='Integrated MCMC noran.pdf')
##save full chains for posterior prediction of abundances
saveRDS(samp, 'MCMC_Null.rds')
saveRDS(summ, 'Summary_Null.rds')
