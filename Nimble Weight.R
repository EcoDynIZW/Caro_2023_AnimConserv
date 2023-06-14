Weight.all<-nimbleCode({
  
  ##priors, intercept and coefficient
  b0~dnorm(0, sd=2)
  b~dnorm(0, sd=2)
  b.sex~dnorm(0, sd=2) #effect of being female
  
  ## random site effect for avg weight
  sd.eps~dgamma(0.1, 0.1)
  ## all surveys at a site get same RE value, epsilon
  for (k in 1:nsites){
    epsilon[k]~dnorm(0, sd=sd.eps)
    ##residual SD by site, to account for heteroscedasticity
    ##in sqrt-transformed data
    l.sd.res[k]~dnorm(mu.s, sd=sig.s)
    sd.res[k]<-exp(l.sd.res[k])
  }
  mu.s~dnorm(0, sd=1)
  sig.s~dgamma(0.1, 0.1)

  
  ##loop over all surveys 
  ##covariates are the same for all surveys within a site
  
  for (j in 1:n){
    # mean depends on area-level covariate and site RE value,
    # plus sex of individual
    # site index goes from 1 to 29 (number of sites) 
    
    mu[j]<-b0 + b*COV[site.index[j]] + b.sex*sex[j] + epsilon[site.index[j]] 
    y[j]~dnorm(mu[j],sd=sd.res[site.index[j]])
    
  }
  
})
