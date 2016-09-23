library(R2jags)
library(MCMCglmm)

## Function to create initial values
##############################################%%
F_initPlsparse<-function(MaxLitter=4,mums=mums,mother=mother){
  mothernb<-length(mums)
  pld<-matrix(data=0,nrow=nind,ncol=MaxLitter)
  pli<-matrix(data=NA,nrow=nind,ncol=MaxLitter)
  
  for(m in 1:mothernb)#this assumes that the first individual must come from the first litter!! (as it should be indeed now that data are re-ordered)
  {
    ind<-which(mother==m)
    for (i in 1:length(ind))
    {
      if (i==1)
      {
        pld[ind[i],1:MaxLitter]<-c(1,rep(x=0,times=MaxLitter-1))
        pli[ind[i],1:MaxLitter]<-NA
      }else{
        pld[ind[i],1:MaxLitter]<-NA
        pli[ind[i],1:MaxLitter]<-runif(MaxLitter,0,1)
      }
    }
  }
  
  initmeanLB<-matrix(data = 0,nrow = mothernb,ncol = MaxLitter)
  
  for (i in 1:nrow(initmeanLB))
  {
    initmeanLB[i,1]<-runif(1,0,40)
    for (j in 2: ncol(initmeanLB))
    {
      initmeanLB[i,j]<-runif(1,initmeanLB[i,j-1]+20,initmeanLB[i,j-1]+60)
    }
  }
  toreturn<-list(pld=pld,pli=pli,initmeanLB=initmeanLB)
  return(toreturn)
}

### JAGS model ###
##############################
sink("models/GrowthViability")
cat("
    model {
    ###################################
    ######priors and constraints#######
    ###################################
    #####survival######
    meanmu~dnorm(0,0.001)
    mean.phi<-1/(1+exp(-meanmu))#logit-1 transform
    BetaA~dnorm(0,0.001)#selection coefficient for asymptotic size
    BetaD~dnorm(0,0.001)#selection coefficient for time to first snow fall
    BetaAD~dnorm(0,0.001)#selection coefficient for interaction asymptotic size and time to first snow fall
    BetaS~dnorm(0,0.001)#sex difference in survival
    
    #####growth#####
    for (i in 1:nind)#loop over individuals
    {
    A[i]~dunif(mA-5,mA+16)#individual asymptotic size
    k[i]~dnorm(mk,tauk)T(0.01,0.025) #individual growth rate
    }
    tauk<-pow(sdk,-2)#precision of individual growth rates
    sdk~dunif(0,0.001)#standard deviation of individual growth rates
    
    tau<-pow(sdMass,-2)#precision in mass measurment
    sdMass<-2.05#standard deviation in mass measurment, estimated from animals measured multiple times
    
    ## Birth dates
    for (mum in 1:mothernb)#loop over mothers*year
    {
    #first litter
    meanLB[mum,1]~dunif(-30,maxbd)  #date of first litter of the mother on this year (between -30, i.e. April 26th and October 5th)
    for (i in 2:MaxLitter)#loop over litters 2 to 5
    {
    meanLB[mum,i]~dunif(gbegin[mum,i],gend[mum,i])#date of successive litters
    gbegin[mum,i]<-meanLB[mum,i-1]+20 #with a minimum of 20 days between successive litters...
    gend[mum,i]<-meanLB[mum,i-1]+120#and a maximum of 120 days
    }
    }
    for (i in 1:nind)
    {
    for(j in 1:MaxLitter)
    {
    pl[i,j]~dgamma(1,1)
    plG[i,j]<-pl[i,j]/sum(pl[i,])
    }
    }
    ###################################
    ############likelihood#############
    ###################################
    for (obs in 1:nobs)
    {
    M[obs]~dnorm(mf[obs],tau)
    mf[obs]<-A[whichind[obs]]*(1-exp(-k[whichind[obs]]*(date[obs])))
    date[obs]<-t[obs]-meanLB[mother[whichind[obs]],litter[whichind[obs]]]
    }
    for (i in 1:nind)
    {
    DeltaWinter[i]<-endSeason[IndCohort[i]]-meanLB[mother[i],litter[i]]
    litter[i] ~ dcat(plG[i,])
    phi[i]~dbin(p[i],1)
    logit(p[i])<-meanmu+BetaA*A[i]+BetaD*DeltaWinter[i]+BetaAD*A[i]*DeltaWinter[i]+BetaS*sex[i]
    pa[i]<-phi[i]*A[i]
    }    
    }
    ",fill = TRUE)
sink()


sparceinitPl<-F_initPlsparse(MaxLitter = MaxLitter,mums = mums,mother = mother)
pld<-sparceinitPl$pld
pli<-sparceinitPl$pli
initmeanLB<-sparceinitPl$initmeanLB

dataGrowthAKlphi2<-list(mothernb=mothernb,M=mass$Weight,endSeason=endSeason,
                        t=mass$RelativeJulian,IndCohort=IndCohort,sex=sex,
                        nind=length(unique(mass$id)),nobs=length(mass$Weight),
                        pl=pld,MaxLitter=MaxLitter,mother=mother,
                        whichind=whichind,maxbd=maxbd,mA=mA,mk=mk,phi=phis,A=Aobs)
initsGAKlphi2 <- function() list(A=initA,k=ki,pl=pli,meanLB=initmeanLB,BetaA=0,BetaD=0,dBetaA=0,BetaS=0,meanmu=0,mu=0.5,tauphi2=0.1)

paramsGAKlphi2 <- c("A","k","litter","meanLB","BetaA","BetaD","BetaAD","meanmu","BetaS","mf","p")
# MCMC settings

GrowthAKlphi2b<-jags(dataGrowthAKlphi2,initsGAKlphi2,paramsGAKlphi2,"models/GrowthViability",
                     n.chains = 3, n.thin = 6000, n.iter = 6300000, n.burnin = 300000, working.directory = getwd())
sumGAK2b<-print(GrowthAKlphi2b)
mean(sumGAK2b$summary[,"Rhat"])
sumGAK2b$summary[which(sumGAK2b$summary[,"Rhat"]>1.1),"Rhat"]

traceplot(GrowthAKlphi2b ,varname="deviance")

GrowthAKlphi2c<-jags(dataGrowthAKlphi2,initsGAKlphi2,paramsGAKlphi2,"models/GrowthViability",
                     n.chains = 3, n.thin = 3000, n.iter = 3300000, n.burnin = 300000, working.directory = getwd())


######################## Posterior Predictive Checks
ppcmean<-vector(length = 3000)
count<-1
for (i in 1:1000)
  for(j in 1:3)
  {
    ppcmean[count]<-mean(mass$Weight) - mean(GrowthAKlphi2c$BUGSoutput$sims.array[i,j,3000+(1:1225)])
    count<-count+1
  }
hist(ppcmean)
mean(ppcmean>0)# posterior predictive check for mean mass p-value= 0.794


ppc<-vector(length = 3000*length(mass$Weight))
count<-1
for (ind in 1: length(mass$Weight))
  for (i in 1: 1000)
    for (j in 1:3)
    {
      ppc[count]<-mass$Weight[ind]-GrowthAKlphi2c$BUGSoutput$sims.array[i,j,3000+ind] 
      count<-count+1
    }

plot(ppc[1:10000])
mean(ppc>0)#posterior predictive check p-value for mass 0.459

names(GrowthAKlphi2c$BUGSoutput$sims.array[1,1,3000:4838])

ppcmeanPhi<-vector(length = 3000)
count<-1
for (i in 1:1000)
  for(j in 1:3)
  {
    ppcmeanPhi[count]<-mean(phis) - mean(GrowthAKlphi2c$BUGSoutput$sims.array[i,j,3000+(1227:1838)])
    count<-count+1
  }
hist(ppcmeanPhi)
mean(ppcmeanPhi>0)# posterior predictive check for mean survival p-value= 0.499