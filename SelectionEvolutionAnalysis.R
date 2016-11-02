
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
priorgBiv1<-list(G=list(G1=list(V=diag(2), nu=2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000),
                        G2=list(V=diag(2), nu=2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000),
                        G3=list(V=diag(2), nu=2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000),
                        G4=list(V=diag(2), nu=2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)),
                 R=list(V=diag(2), nu=2))

mgBivLRS2<-MCMCglmm(cbind(Weight,rLRS2) ~ trait-1 + trait:(f)+at.level(trait,1):(Sex*Age*StDate+StDateSq),
                    random =~ us(trait):animal + us(trait):id + us(trait):Mother + idh(trait):Year, 
                    rcov=~us(trait):units, prior = priorgBiv1, data=AllM,pedigree=ped, family=c("gaussian","gaussian"),
                    verbose=TRUE,nitt=3000000,burnin=1000000,thin=1000)
summary(mgBivLRS2)

#mode and credibility interval of the additive genetic covariance mass-fitness
posterior.mode(mgBivLRS2$VCV[,"Weight:rLRS2.animal"]) ; HPDinterval(mgBivLRS2$VCV[,"Weight:rLRS2.animal"])

#### Full model with all traits (long to run, good to confirm by pair-wise bivariate models like the one above) ####
b <- 1e+10
B1_V <- diag(32)*b; B1_mu <- rep(0,32)
B1_V[4,4]<-0.01
B1_mu[4]<-0.5
B1_V[8,8]<-1
B1_mu[8]<-0
priorg0<-list(G=list(G1=list(V=diag(4), nu=4, alpha.mu=rep(0,4), alpha.V=VA*1000),
                     G2=list(V=diag(4), nu=4, alpha.mu=rep(0,4), alpha.V=Vid*1000),
                     G3=list(V=diag(4), nu=4, alpha.mu=rep(0,4), alpha.V=VM*1000),
                     G4=list(V=diag(3), nu=3, alpha.mu=rep(0,3), alpha.V=VY*1000)),
              R=list(V=VR, nu=4),
              B=list(V=B1_V, mu=B1_mu))

mgLRS0<-MCMCglmm(cbind(Weight,Body_Length,Tail_Length,rLRS) ~ trait-1 + trait:(f)+at.level(trait,c(1:3)):(Sex*Age*StDate+StDateSq),
                 random =~ us(trait):animal + us(trait):id + us(trait):Mother + idh(at.level(trait,c(1:3))):Year, 
                 rcov=~us(trait):units, prior = priorg0,start=startg0, data=AllM,pedigree=ped, family=c("gaussian","gaussian","gaussian","gaussian"),
                 verbose=TRUE,nitt=1300000,burnin=300000,thin=1000,pr=FALSE)

#### Robertson-Price decomposition for fitness components ####
#Adult males ARS
priorARSad <- list(G=list(G1=list(V=diag(2), nu=2),
                        G2=list(V=diag(2), nu=2),
                        G3=list(V=diag(2), nu=2),
                        G4=list(V=diag(2), nu=2)),
                 R=list(V=diag(2),nu=2))


mARSadMNorm <- MCMCglmm(cbind(Weight,rARS)~trait-1 + at.level(trait,1):(StDate+StDateSq),
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
priorPhimadNorm<-list(G=list(G1=list(V=diag(2), nu=2),
                             G2=list(V=diag(1), nu=1),
                             G3=list(V=diag(2), nu=2),
                             G4=list(V=diag(2), nu=2)),
                      R=list(V=diag(2),nu=2,fix=2))

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