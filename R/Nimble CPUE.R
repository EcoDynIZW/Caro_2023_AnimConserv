Cpue.all<-nimbleCode({
  
  ##priors, intercept and coefficient
  b0~dnorm(0, sd=2)
  b~dnorm(0, sd=2)
  
  ## random site effect for avg CPUE
  sd.eps~dgamma(0.1, 0.1)
  
  ## all surveys at a site get same RE value, epsilon
  for (k in 1:nsites){
    epsilon[k]~dnorm(0, sd=sd.eps)
    }

  ##NB overdispersion
  r~dunif(0,50)
  
  for (j in 1:nsurveys){
    # survey specific mean depends on area-level covariate and RE value,
    # site index goes from 1 to 29 (number of sites) and needs to be used 
    # instead of the actual site number, as the max site number exceeds 
    # number of sites in the model
    # log(effort)=offset transforming count into count/unit effort

    log(lam.survey[j])<-b0 + b*COV[j] + epsilon[site.index[j]] +
                        +log(effort[j])
    
    #Ccount is NegBin RV
    pp[j]<-r/(r+lam.survey[j])
    y[j]~dnegbin(pp[j], r)

    # #per visit data - leads to estimation problems
    # for (k in 1:nvisits[j]){
    #   log(lam.visit[j,k])<-log(lam.survey[j])+log(effort[j,k])
    #   pp[j,k]<-r/(r+lam.visit[j,k])
    #   y[j,k]~dnegbin(pp[j,k],r)
    # }
    
  }
  
  })
