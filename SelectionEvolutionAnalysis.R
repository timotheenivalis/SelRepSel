
####Packages and data #####
library(lme4)
library(reshape)#for merge
library(car)
library(MCMCglmm)
library(boot)
library(pedantics)

AllM<-read.table(file = "AllM.txt",header=T)

####Phenotypic selection estimation####

####Phenotypic selection test####


#####Quantitative genetics####
#### Blups analysis ####

#### Main model for Robertson-Price equation ####
mgBivLRS2<-MCMCglmm(cbind(Weight,rLRS2) ~ trait-1 + trait:(f)+at.level(trait,1):(Sex*Age*StDate+StDateSq),
                    random =~ us(trait):animal + us(trait):id + us(trait):Mother + idh(trait):Year, 
                    rcov=~us(trait):units, prior = priorgBiv1, data=AllM,pedigree=ped, family=c("gaussian","gaussian"),
                    verbose=TRUE,nitt=3000000,burnin=1000000,thin=1000)
summary(mgBivLRS2)

#mode and credibility interval of the additive genetic covariance mass-fitness
posterior.mode(mgBivLRS2$VCV[,"Weight:rLRS2.animal"]) ; HPDinterval(mgBivLRS2$VCV[,"Weight:rLRS2.animal"])

#### Robertson-Price decomposition for fitness components ####
#Adult males ARS
mARSadMNorm<-MCMCglmm(cbind(Weight,rARS)~trait-1 + at.level(trait,1):(StDate+StDateSq),
                      random=~us(trait):animal+us(trait):id+us(trait):Mother+idh(trait):Year, rcov=~us(trait):units,
                      pedigree=ped,data=AllM[which( AllM$Age=="A" & AllM$Sex=="Male" ),], 
                      prior = priorARSad,family=c("gaussian","gaussian"),verbose=TRUE,nitt=110000,burnin=10000,thin=100)
summary(mARSadMNorm)

#Adult females ARS
mARSadFNorm<-MCMCglmm(cbind(Weight,rARS)~trait-1 + at.level(trait,1):(StDate+StDateSq),
                      random=~us(trait):animal+us(trait):id+us(trait):Mother+idh(trait):Year, rcov=~us(trait):units,
                      pedigree=ped,data=AllM[which( AllM$Age=="A" & AllM$Sex=="Female" ),], 
                      prior = priorARSad,family=c("gaussian","gaussian"),verbose=TRUE,nitt=110000,burnin=10000,thin=100)
summary(mARSadFNorm)

#Juvenile survival
mRelphimjuvNorm<-MCMCglmm(cbind(Weight,rPhi)~trait-1 + at.level(trait,1):(Sex*StDate+StDateSq+Age3+f)+ at.level(trait,2):(Sex+f),
                          random=~us(trait):animal+us(at.level(trait,1)):id+us(trait):Mother+idh(trait):Year, rcov=~us(trait):units,
                          pedigree=ped,data=AllM[which( AllM$Age=="J"),],
                          prior = priorPhimadNorm,family=c("gaussian","gaussian"),verbose=TRUE,nitt=110000,burnin=10000,thin=100)
summary(mRelphimjuvNorm)

#Adult survival
mRelphimadNorm<-MCMCglmm(cbind(Weight,rPhi)~trait-1 + at.level(trait,1):(Sex*StDate+StDateSq+Age3+f)+ at.level(trait,2):(Sex+f),
                         random=~us(trait):animal+us(at.level(trait,1)):id+us(trait):Mother+idh(trait):Year, rcov=~us(trait):units,
                         pedigree=ped,data=AllM[which(!duplicated(paste(AllM$id,AllM$Year)) & AllM$Age=="A"),],
                         prior = priorPhimadNorm,family=c("gaussian","gaussian"),verbose=TRUE,nitt=110000,burnin=10000,thin=100)
summary(mRelphimadNorm)