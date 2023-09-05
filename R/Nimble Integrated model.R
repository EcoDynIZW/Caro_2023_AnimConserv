
Integrated.model0<-nimbleCode({
  
  ##priors shared density model
  ##intercept density
  b0~dnorm(0, sd=2)
  
  ## random site effect for avg density
  sd.eps~dgamma(0.1, 0.1)
  sd.psi~dgamma(0.1, 0.1)
  
  for (j in 1:nsites){
    eps.site[j]~dnorm(0, sd=sd.eps)
  }
  
  #priors, site-specific detection parameters for CMR models
  for (j in 1:nsites.cmr){
    p.int[j]~dunif(0, 1)
    lp.int[j]<-log(p.int[j]/(1-p.int[j]))
    lp.b[j]~dnorm(0, sd=2)
  }
  
  # prior, overdispersion in count model
  r ~ dunif(0,50)
  
  ##############################################################################
  #### CMR model
  
  ## loop over CMR surveys
  for (j in 1:nsurv.cmr){
    #priors inclusion probability per survey
    #calculated from expected abundance and data augmentation M
    #random site effect and b0 shared with count model
    #site.index.cmr links CMR site to overall site ID
    #and is the same for all Misali/Chumbe surveys regardless of observer
    log(lam.surv.cmr[j])<-b0 + log(area.cmr[j]) + eps.site[site.index.cmr[j]]
    
    mu.psi[j]<-get.psi(lam=lam.surv.cmr[j],
                    M=M[j])
    psi[j]~dbeta(mean=mu.psi[j], sd=sd.psi)
    
    ##site is really site-observer and references the correct detection parameter
    logit(p.cmr[j,1:K.cmr[j]])<-lp.int[site[j]]+lp.b[site[j]]*effort.cmr[j,1:K.cmr[j]]
    
    for ( i in 1:M[j]){ #augment only as needed
      z[i,j]~dbern(psi[j])
      y.cmr[i,j,1:K.cmr[j]]~dbern_V(p=p.cmr[j,1:K.cmr[j]],
                            indicator=z[i,j])
      
    }
    N.cmr[j]<-sum(z[1:M[j],j])
  }
  
  ##############################################################
  ### count model

  for (j in 1:nsurv.count){
    ##get survey-specific detection parameters:
    ##For sites/surveys with CMR survey, use correposnding CMR estimate
    ##for other sites, use average over Tim's CMR surveys only
    ##calculated in custom function below
    det.par[j,1:2]<-get.det.parms(cmr.index = cmr.index[j],
                               lp.int=lp.int[1:nsites.cmr],
                               lp.b=lp.b[1:nsites.cmr])

    # survey specific mean density depends on area-level covariate and RE value,
    # and survey specific area sampled
    # site index goes from 1 to 29 (number of sites) 
    log(lam.surv.count[j])<-b0 + log(area.count[j]) + eps.site[site.index.count[j]]
                                                      
    # abundance for a survey is NegBin RV
    pp[j]<-r/(r+lam.surv.count[j]) 
    N.count[j]~dnegbin(pp[j],r)
    
    for (k in 1:K.count[j]){
      # detection probability at a visit depends on effort
      logit(p.count[j,k])<-det.par[j,1]+det.par[j,2]*effort.count[j,k]
      # crabs observed is Binomial RV
      y.count[j,k]~dbinom(p.count[j,k], N.count[j])
    }
  }
})

################################################################################
Integrated.model<-nimbleCode({
  
  ##priors shared density model
  ##intercept, cov effect density
  b0~dnorm(0, sd=2)
  b~dnorm(0, sd=2)
  
  ## random site effect for avg density
  sd.eps~dgamma(0.1, 0.1)
  sd.psi~dgamma(0.1, 0.1)
  
  for (j in 1:nsites){
    eps.site[j]~dnorm(0, sd=sd.eps) 
  }
  
  #priors, site-specific detection parameters for CMR models
  for (j in 1:nsites.cmr){
    p.int[j]~dunif(0, 1)
    lp.int[j]<-log(p.int[j]/(1-p.int[j]))
    lp.b[j]~dnorm(0, sd=2)
  }
  
  # prior, overdispersion in count model
  r ~ dunif(0,50)
  
  ##############################################################################
  #### CMR model
  
  ## loop over CMR surveys
  for (j in 1:nsurv.cmr){
    #priors inclusion probability per survey
    #calculated from expected abundance and data augmentation M
    #random site effect and b0 shared with count model
    #site.index.cmr links CMR site to overall site ID
    #and is the same for all Misali/Chumbe surveys regardless of observer
    log(lam.surv.cmr[j])<-b0 + b*cov.cmr[j] + log(area.cmr[j]) + eps.site[site.index.cmr[j]]
    
    mu.psi[j]<-get.psi(lam=lam.surv.cmr[j],
                       M=M[j])
    psi[j]~dbeta(mean=mu.psi[j], sd=sd.psi)
    
    ##site is really site-observer and references the correct detection parameter
    logit(p.cmr[j,1:K.cmr[j]])<-lp.int[site[j]]+lp.b[site[j]]*effort.cmr[j,1:K.cmr[j]]
    
    for ( i in 1:M[j]){ #augment only as needed
      z[i,j]~dbern(psi[j])
      y.cmr[i,j,1:K.cmr[j]]~dbern_V(p=p.cmr[j,1:K.cmr[j]],
                                    indicator=z[i,j])
      
    }
    N.cmr[j]<-sum(z[1:M[j],j])
  }
  
  ##############################################################
  ### count model
  
  for (j in 1:nsurv.count){
    ##get survey-specific detection parameters:
    ##For sites/surveys with CMR survey, use corresponding CMR estimate
    ##for other sites, use average over Tim's CMR surveys only
    ##calculated in custom function below
    det.par[j,1:2]<-get.det.parms(cmr.index = cmr.index[j],
                                  lp.int=lp.int[1:nsites.cmr],
                                  lp.b=lp.b[1:nsites.cmr])

    # survey specific mean depends on area-level covariate and RE value,
    # and survey specific area sampled
    # site index goes from 1 to 29 (number of sites) 
    log(lam.surv.count[j])<-b0 + b*cov.count[j] + log(area.count[j]) + eps.site[site.index.count[j]]
    
    # abundance for a survey is NegBin RV
    pp[j]<-r/(r+lam.surv.count[j]) 
    N.count[j]~dnegbin(pp[j],r)
    
    for (k in 1:K.count[j]){
      # detection probability at a visit depends on effort
      logit(p.count[j,k])<-det.par[j,1]+det.par[j,2]*effort.count[j,k]
      # crabs observed is Binomial RV
      y.count[j,k]~dbinom(p.count[j,k], N.count[j])
    }
  }
})


################################################################################
#### custom function ###########################################################

## vectorize dbern to increase efficiency

dbern_V <- nimbleFunction( run = function(x = double(1), ##observations
                                          p=double(1),   ##p
                                          indicator = double(0), ##psi
                                          log = double(0)) {
  
  returnType(double())
  
  #shortcut for individuals not alive
  if(indicator == 0){
    if (sum(x)==0){
      if (log == 0)return(1.0)
      else return(0.0)
    } else {
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  nocc<-length(x)
  
  prob1<-x*p  
  prob0<-(1-x)*(1-p) 
  prob<-prod(prob1+prob0)
  
  if(log) return(log(prob))
  return(prob)
  
})

rbern_V <- nimbleFunction( run = function(n = integer(0), ##observations
                                          p=double(1),   ##p
                                          indicator = double(0)) {
  
  returnType(double(1))
  k<-length(p)
  #shortcut for individuals not alive
  if(indicator == 0){
    return(rep(0,k))
  }
  
  
  x<-rbinom(k,1,p)
  return(x)
  
})

registerDistributions(list(dbern_V = list(
  BUGSdist = "dbern_V(p, indicator)",
  Rdist = "dbern_V(p, indicator)",
  types = c('value = double(1)',
            'p = double(1)',
            'indicator = double(0)'
  )
))
)

#### function to get psi and restrict it to <1
# ensures that algorithm does not produce an error, but still requires checking that M 
# is sufficiently large not to constrain estimates of psi/abundance; look at estimates
# of psi in model output and make sure even 97.5th percentile of posterior is well below 1

get.psi<-nimbleFunction(run = function(lam=double(0),   ##expected abundance
                                      M = double(0)) {
  
  returnType(double(0))
  x<-lam/M
  if(x>=1)x<-0.9999
  
  return(x)
  
})


### function to get survey-specific detection parameters
# site-specific CMR estimate (Tim) if count survey took place at a site with CMR survey
# average over Tim's CMR estimate if count survey took place at a site without CMR survey

get.det.parms <- nimbleFunction( run = function(cmr.index=double(0),  ##index which CMR site survey corresponds to, -999 if non-CMR site (here, site means site-observer combination)
                                                lp.int=double(1),     ##detection intercept, logit scale
                                                lp.b=double(1)       ##effect of effort on p
) {
  
  returnType(double(1))
  
  #for non-CMR sites, use average over Tim's surveys
  if (cmr.index == -999){
    lpind<-mean(lp.int[1:5])
    lpb<-mean(lp.b[1:5])
  } else {
    
    lpind<-lp.int[cmr.index]
    lpb<-lp.b[cmr.index]
  }
  
  return(c(lpind, lpb))
  
})
