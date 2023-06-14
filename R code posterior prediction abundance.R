########################################################################################
##### Posterior predictions of coconut crab population size in all suitable habitat ####

### From Caro et al. (in press). Meta- and subpopulation estimation with 
###      disparate data: coconut crabs in the Western Indian Ocean.
###      Animal Conservation

## Use expected density from integrated model for all sites, CMR and non-CMR; 
## extrapolate based on Null model (but estimates are very consistent across models)
## Script requires running integrated model first (script 'R code integrated model.R')
## Also, all necessary data objects must be in working directory, or else, paths need
## to be adjusted accordingly

###### The data object that is read in below is a data frame with the following columns:
## $Number: numerical index (1-29) for all surveyed sites
## $Name: name of all surveyed sites
## $Total area: total area of suitable crab habitat (ha) for each site (see main
##              text for details)

rm(list=ls())
library(MCMCvis)
library(writexl)
library(readxl)

get.mode<-function(x){
  xx<-density(x)
  out<-xx$x[xx$y == max(xx$y)]
  return(out)
}

##load MCMC samples of integrated model without covariates 
##(produced with script 'R code integrated model.R')
samp0<-readRDS('MCMC_Null.rds')

##get site level area data and names (for results table)
site.data<-readRDS('Posterior prediction site data.rds')

nsites<-nrow(site.data) # number of sampled sites

################################################################################

## Step 1: Extrapolate abundance to area of suitable habitat for all sites

samp.mat<-do.call(rbind, samp0)
nkeep<-5000 #perform posterior prediction for 5000 random posterior samples
set.seed(123)
keep<-sample(1:nrow(samp.mat),nkeep, replace = FALSE)

##calculate density based on integrated model parameters, from that, extrapolate total N
##works for all sites except Chumbe (see below)
D.post0<-matrix(NA, nkeep, nsites)
A.post0<-matrix(NA, nkeep, nsites)

for (j in 1:nsites){

  #get total area (ha) of suitable habitat for sampled site
  arr<-site.data[j,"Total area"]
  #calculate expected density from integrated model parameters
  D.post0[,j]<-exp(samp.mat[keep,'b0'] + 
                     samp.mat[keep,paste('eps.site[', j, ']', sep='')])
  #multiply with area to obtain abundance
  A.post0[,j]<-D.post0[,j]*arr
  
}

##get psoterior summaries for abundance and density, including posterior mode
D.summ0<-MCMCsummary(list(D.post0), Rhat=FALSE,func=get.mode,
                     probs = c(0.025,0.05,0.25, 0.5, 0.75,0.95, 0.975))

A.summ0<-MCMCsummary(list(A.post0), Rhat=FALSE,func=get.mode,
                     probs = c(0.025,0.05,0.25, 0.5, 0.75,0.95, 0.975))



################################################################################

## Step 2: Calculate N for Chumbe (area=9) based on reported proportions of 
##         individuals caught in lodge area vs rest of the island by Killströmer 
##         & Bergwall 2013 (for details see Methods)

D.chumbe<-D.post0[,9]
pKitchen<-0.738 #according to paper, proportion of abundance in Kitchen area
pOther<-1-pKitchen
area.other<-site.data[9,'Total area'] - 3.7772 #3.78 is sampled area = lodge area

A.chumbe<-D.chumbe*3.7772 + #abundance lodge area
  D.chumbe*3.7772*(pOther/pKitchen) #abundance rest of island
A.chumbe.sum<-MCMCsummary(list(matrix(A.chumbe, nrow=nkeep, ncol=1)), Rhat=FALSE,
                          probs = c(0.025,0.05,0.25, 0.5, 0.75,0.95, 0.975), func=get.mode)

#replace Chumbe results in overall summary, and in posterior samples
A.summ0[9,]<-A.chumbe.sum
A.post0[,9]<-A.chumbe


################################################################################

##Step 3: Make output table and get total population for Pemba

##prep output table
full.results<-data.frame(Name=site.data$Name, 
                         Area=round(site.data$`Total area`, dig=2),
                         Pemba=1)

#mark sites not on Pemba
not.pemba<-which(full.results$Name %in% c('Chumbe', "Bongoyo (Dar es Salaam)",
                                          "Mbudya (Dar es Salaaam)"))
full.results[not.pemba, 'Pemba']<-0
new.cols<-c('D','5%.D', '95%.D','N','5%.N', '95%.N')
full.results[,new.cols]<-NA

##fill with estimates from A.summ0 and D.summ0
full.results[, 'D']<-round(D.summ0$func, dig=2)
full.results[, '5%.D']<-round(D.summ0$`5%`, dig=2)
full.results[, '95%.D']<-round(D.summ0$`95%`, dig=2)
full.results[, 'N']<-round(A.summ0$func, dig=2)
full.results[, '5%.N']<-round(A.summ0$`5%`, dig=2)
full.results[, '95%.N']<-round(A.summ0$`95%`, dig=2)

##get N for Pemba
N.total<-apply(A.post0[,full.results$Pemba ==1],1,sum) 

N.t.summ<-MCMCsummary(list(matrix(N.total, nkeep, 1)), Rhat=FALSE,func=get.mode,
                      probs = c(0.025,0.05,0.25, 0.5, 0.75,0.95, 0.975))

res.full<-data.frame(matrix(nrow=1, ncol=ncol(full.results)))
colnames(res.full)<-colnames(full.results)
res.full[,'Name']<-'Pemba archipelago'
res.full[,c("N","5%.N","95%.N")]<-N.t.summ[,c("func","5%","95%")]
full.results<-rbind(full.results, res.full)
write_xlsx(full.results, 'Abundance and density.xlsx')

