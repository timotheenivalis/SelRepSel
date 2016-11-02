
#### Graphe 1A : phenotypic trend ####
par(cex=1.5,mar=c(5, 4.5, 4, 2) + 0.1)
hbar<-0.1
plot(yearmean[,1],x=yearmean$year,ylim=c(min(yearmean[,1:3]),max(yearmean[,1:3])),xlab="Year",ylab="",pch=16,xaxt="n",yaxt="n",main="\\textbf{\\Large A}")
mtext(text = c("Age/sex corrected mass (g)"),side = 2,line = 3.5,cex=1.9)
axis(side = 1,at = yearmean$year)
axis(side = 2,at = 37:42,labels = 37:42,las=1,padj = 0.2)
segments(x0 = yearmean$year,y0 = yearmean[,2],y1 = yearmean[,3])
segments(x0 = yearmean$year+hbar,y0 = yearmean[,2],x1 = yearmean$year-hbar)
segments(x0 = yearmean$year+hbar,y0 = yearmean[,3],x1 = yearmean$year-hbar)

lines(y=newdatRef$Weight,x=newdatRef$Year,lwd=6,lty=2,col=mycolors3[1])
lines(y=newdatRef$plo,x=newdatRef$Year,lwd=3,lty=2,col=mycolors3[1])
lines(y=newdatRef$phi,x=newdatRef$Year,lwd=3,lty=2,col=mycolors3[1])

#### Graphe 1B : trend in breeding values ####
par(cex=1.5,mar=c(5, 4.5, 4, 2) + 0.1)
hbar<-0.1
plot(x=-1:8,y=posterior.mode(t(yearsBV)),ylim=c(min(hpdi),max(hpdi)),xaxt="n",yaxt="n",xlab="Cohort",ylab="",pch=16,main="\\textbf{ \\Large B}")
axis(side = 1,at = -1:8,label=FALSE)

axis(side=1,at=0:8,label=c(2006:2014),tick = FALSE)
axis(side = 2,at = seq(-1,1,0.5),las=1,padj=0)
mtext(text = c("Mean breeding value for mass (g)"),side = 2,line = 3.5,cex=1.9)
segments(x0=-1:8,y0 = hpdi[,1],y1=hpdi[,2])
segments(x0 = -1:8+hbar,y0 = hpdi[,1],x1 = -1:8-hbar)
segments(x0 = -1:8+hbar,y0 = hpdi[,2],x1 = -1:8-hbar)
lines(x=mmy[,2],y=bvnewdat$newdat,lwd=6,col=rgb(255,200,5,maxColorValue = 255))
lines(x=mmy[,2],y=bvnewdat$plo,lwd=3,col=rgb(255,200,5,maxColorValue = 255))
lines(x=mmy[,2],y=bvnewdat$phi,lwd=3,col=rgb(255,200,5,maxColorValue = 255))

#### Graphe 2A : All predictions / estimations of evolution ####
par(mar=c(4,6,2,2)+0.1,cex=1.5)
hbar<-0.1

pred<-c(R[1],Ri[1],posterior.mode(h2S),fixef(mwytrend)["Year"],posterior.mode(lmBV[,2]),posterior.mode(mgBivLRS2$VCV[,"Weight:rLRS2.animal"]/gentime))
CI<-matrix(c(HPDinterval(Rm)-mean(Rm)+R[1],HPDinterval(Rm)-mean(Rm)+Ri[1],HPDinterval(h2S),fixef(mwytrend)["Year"]-1.96*sqrt(vcov(mwytrend)["Year","Year"]),fixef(mwytrend)["Year"]+1.96*sqrt(vcov(mwytrend)["Year","Year"]),HPDinterval(lmBV[,2]),HPDinterval(mgBivLRS2$VCV[,"Weight:rLRS2.animal"]/gentime)),ncol=2,byrow=TRUE)
x=barplot(pred,ylim=c(-0.5,0.4),names.arg = c("MBE","MBE$_{\\rho=0}$","UBE","PT","TPBV","GCPE"),axes = FALSE,col = mycolors7,border=mycolors7,main="\\textbf{\\Large A}")
segments(x0=x,y0=CI[,1],y1=CI[,2],lwd=3)
segments(x0=x+hbar,y0=CI[,1],x1=x-hbar,lwd=3)
segments(x0=x+hbar,y0=CI[,2],x1=x-hbar,lwd=3)
axis(side = 2,at = round(seq(-0.5,0.4,0.2),1),las=2,padj=0)
mtext(side = 2, "Predicted change in mean mass (g/year)", line = 4,cex = 2)
abline(h=0)

#### Graphe 2B : Robertson-Price decomposition of fitness components #####
SrLRS<-mgBivLRS2$VCV[,"Weight:rLRS2.animal"]+mgBivLRS2$VCV[,"Weight:rLRS2.id"]+mgBivLRS2$VCV[,"Weight:rLRS2.Mother"]+mgBivLRS2$VCV[,"Weight:rLRS2.units"]
ErLRS<-mgBivLRS2$VCV[,"Weight:rLRS2.id"]+mgBivLRS2$VCV[,"Weight:rLRS2.Mother"]+mgBivLRS2$VCV[,"Weight:rLRS2.units"]
GrLRS<-mgBivLRS2$VCV[,"Weight:rLRS2.animal"]

SARSm<-mARSadMNorm$VCV[,"rARS:Weight.animal"]+mARSadMNorm$VCV[,"rARS:Weight.id"]+mARSadMNorm$VCV[,"rARS:Weight.Mother"]+mARSadMNorm$VCV[,"rARS:Weight.units"]
EARSm<-mARSadMNorm$VCV[,"rARS:Weight.id"]+mARSadMNorm$VCV[,"rARS:Weight.Mother"]+mARSadMNorm$VCV[,"rARS:Weight.units"]
GARSm<-mARSadMNorm$VCV[,"rARS:Weight.animal"]

SARSf<-mARSadFNorm$VCV[,"rARS:Weight.animal"]+mARSadFNorm$VCV[,"rARS:Weight.id"]+mARSadFNorm$VCV[,"rARS:Weight.Mother"]+mARSadFNorm$VCV[,"rARS:Weight.units"]
EARSf<-mARSadFNorm$VCV[,"rARS:Weight.id"]+mARSadFNorm$VCV[,"rARS:Weight.Mother"]+mARSadFNorm$VCV[,"rARS:Weight.units"]
GARSf<-mARSadFNorm$VCV[,"rARS:Weight.animal"]

Sphij<-mRelphimjuvNorm$VCV[,"rPhi:Weight.animal"]+mRelphimjuvNorm$VCV[,"rPhi:Weight.Mother"]+mRelphimjuvNorm$VCV[,"rPhi:Weight.units"]
Ephij<-mRelphimjuvNorm$VCV[,"rPhi:Weight.Mother"]+mRelphimjuvNorm$VCV[,"rPhi:Weight.units"]
Gphij<-mRelphimjuvNorm$VCV[,"rPhi:Weight.animal"]

Sphia<-mRelphimadNorm$VCV[,"rPhi:Weight.animal"]+mRelphimadNorm$VCV[,"rPhi:Weight.Mother"]+mRelphimadNorm$VCV[,"rPhi:Weight.units"]
Ephia<-mRelphimadNorm$VCV[,"rPhi:Weight.Mother"]+mRelphimadNorm$VCV[,"rPhi:Weight.units"]
Gphia<-mRelphimadNorm$VCV[,"rPhi:Weight.animal"]

selectcomp<-matrix(data = sapply(list(SrLRS,ErLRS,GrLRS,SARSm,EARSm,GARSm,SARSf,EARSf,GARSf,Sphij,Ephij,Gphij,Sphia,Ephia,Gphia),mean),nrow= 3)

CI<-matrix(data = sapply(list(SrLRS,ErLRS,GrLRS,SARSm,EARSm,GARSm,SARSf,EARSf,GARSf,Sphij,Ephij,Gphij,Sphia,Ephia,Gphia),HPDinterval),ncol= 2,byrow=TRUE)


par(mar=c(4,6,4,2)+0.1,cex=1.5)
hbar<-0.1
x=barplot(selectcomp,cex.names = 1.4,cex.lab=1.4,cex.axis = 1.4,beside=TRUE,names.arg = c("LRS","ARS$_\\textrm{\\male}$","ARS$_\\textrm{\\female}$","$\\phi_{\\mathrm{Juv}}$","$\\phi_{\\mathrm{Ad}}$"),las=1,ylim=c(-3.25,5.4),ylab="Selection differential (g)",col = mycolors3, border=mycolors3, main="\\textbf{ \\Large B}")
legend(x = 14,y=3.8,legend=c("Phenotypic","Environmental","Genetic"),fill = mycolors3,bty="n",xjust = 0,yjust= 0.5, x.intersp = 0.5)
segments(x0=x,y0=CI[,1],y1=CI[,2],lwd=3)
segments(x0=x+hbar,y0=CI[,1],x1=x-hbar,lwd=3)
segments(x0=x+hbar,y0=CI[,2],x1=x-hbar,lwd=3)
abline(h=0,lwd=1)


#### Graphe 3B #####
par(cex=2,mar=c(5, 5, 4, 1) + 0.1, cex.axis=1.3, cex.lab=1.3)
plot(x = 2006:2013,y=rep(0,8),ylim=c(90,350),pch=16,col="white",xlab="\\textbf{Year}",ylab="\\textbf{Julian day}",las=1,xaxt="n",yaxt="n",main="\\textbf{\\Large B}")
axis(side = 1,at = seq(2006,2013,2),cex=1.4)
axis(side=2,at=seq(100,350,50),padj = 0.5,cex=1.4,las=1)
firstfree0<-2006:2013
for (y in 2006:2013)
{
  thisyear<-rawsnow[which(rawsnow$year==y ),]#& rawsnow$julian>firstfree[y-2005] & rawsnow$julian<lastfree[y-2005]),]
  firstfree0[y-2005]<-min(thisyear$julian[which(thisyear$hto000j0==0)])
  points(thisyear$julian,x=rep(y,nrow(thisyear)),
         col=rgb(red = 0,green = 0,blue = 1,alpha =(1-(max(rawsnow$hto000j0,na.rm=T)-thisyear$hto000j0)/max(rawsnow$hto000j0,na.rm=T))^(1/3)),pch=16,cex=0.5+5*(1-(max(rawsnow$hto000j0,na.rm=T)-thisyear$hto000j0)/max(rawsnow$hto000j0,na.rm=T)))
  }

points(x=IndBD$cohort,y=IndBD$bd+min(AllMj$Julian))

#### Graphe 3C #####
af<-c(114,128,126,128,114,109,130,113,114,99,130,123,115,92,121,114,110)
al<-c(318,276,308,325,308,277,324,324,325,312,276,286,292,261,288,282,294)

par(cex=2,mar=c(5, 5, 4, 1) + 0.1, cex.axis=1.3, cex.lab=1.3)
freeSnowSpan<-al-af
years<-1998:2014
plot(freeSnowSpan[-c(1,2,3)],x=years[-c(1,2,3)],xlab="\\textbf{Year}",ylab="\\textbf{Snow-free season length (days)}",pch=16,xaxt="n",yaxt="n",cex.lab=1.4,ylim=c(128,215),main="\\textbf{\\Large C}")
axis(side = 1,at = seq(2001,2014,4),cex=1.4)
axis(side=2,at=seq(110,210,20),padj = 0.5,cex=1.4,las=1)
abline(v = 2005.5,lty=2,lwd=5,col=gray(0.5))
arrows(x0 = 2007.7,y0 = 140,x1 = 2014.3,code = 3,lwd=5)
arrows(x0 = 2000.7,y0 = 140,x1 = 2007.3,code = 3,lwd=5)
text(x=c(2004,2011,2007.7),y=c(140,140,134)-5,labels=c("\\textbf{Period a}","\\textbf{Period b}","\\textbf{a-b=$\\boldsymbol{-23}$ days, p=0.02}"))

arrows(x0 = 2005.7,y0 = 195,x1 = 2014.3,code = 3,lwd=5,col=gray(0.2))
text(x=c(2010),y=c(195)-5,labels=c("\\textbf{Snow vole monitoring}"),col=gray(0.2))


#### Graphe 3D #####
par(cex=2,mar=c(4,6,2,2)+0.1, cex.axis=1.3, cex.lab=1.3)
fullArray<-rbind(GrowthAKlphi2b$BUGSoutput$sims.array[,1,],GrowthAKlphi2b$BUGSoutput$sims.array[,2,],GrowthAKlphi2b$BUGSoutput$sims.array[,3,])
vcovbeta<-var(fullArray[,c("meanmu","BetaA","BetaD","BetaS","BetaAD")])
fixmeans<-matrix(data=as.numeric(GrowthAKlphi2b$BUGSoutput$mean[c("meanmu","BetaA","BetaD","BetaS","BetaAD")]))
newdat <- expand.grid(A=30:52,D=0:183,S=0:1)

matrix(data=as.numeric(GrowthAKlphi2b$BUGSoutput$sd[c("meanmu","BetaA","BetaD","BetaS","BetaAD")]))

mm <- as.matrix(data.frame(intercept=1,newdat,AD=newdat$A*newdat$D))#design matrix
newdat$logitphi <-  mm %*% fixmeans
pvar1 <- diag(mm %*% tcrossprod(vcovbeta,mm))
newdat <- data.frame(
  newdat
  , plo = newdat$logitphi-1.96*sqrt(pvar1)
  , phi = newdat$logitphi+1.96*sqrt(pvar1)
)
newdatRef<-newdat[which(newdat$S==0 & newdat$D==81),]

newdatRefnoSel<-newdat[which(newdat$S==0 & newdat$D==120),]


upperMassJW=max(SurvWinter$MassCorected[which(SurvWinter$Age=="J")])
lowerMassJW=min(SurvWinter$MassCorected[which(SurvWinter$Age=="J")])

xvJW<-seq(lowerMassJW,upperMassJW,0.1)
xvJWS<-(xvJW-mean(SurvWinter$MassCorected))/sd(SurvWinter$MassCorected)

pj<-(fixef(mlrwinter)[1]+(fixef(mlrwinter)[2]+fixef(mlrwinter)[4])*xvJW+fixef(mlrwinter)[3])
yvj<-exp(pj)/(1+exp(pj))

newdatM <- expand.grid(survival=0,MassCorected=seq(from=8,to=62),Age=c("A","J"))
mm <- model.matrix(terms(mlrwinter),newdatM)
newdatM$survival <- mm %*% fixef(mlrwinter)
pvar1 <- diag(mm %*% tcrossprod(vcov(mlrwinter),mm))

newdatM <- data.frame(
  newdatM
  , plo = newdatM$survival-1.96*sqrt(pvar1)
  , phi = newdatM$survival+1.96*sqrt(pvar1)
) 
newdatM$survival<-inv.logit(newdatM$survival)
newdatM$plo<-inv.logit(newdatM$plo)
newdatM$phi<-inv.logit(newdatM$phi)

rangeJ<-round(lowerMassJW):upperMassJW
faclwd <- 1.2

plot(x=xvJW,y=yvj,col="blue",ylim=c(-0.01,1.01),xlim=c(5,52),type="n",lwd=8,xlab="\\textbf{Real or Asymptotic mass (g)}",ylab="\\textbf{Survival probability}",las=1, main="\\textbf{ \\large D}")
lines(x=xvJW,y=yvj,col=mycolors3[1],lwd=8*faclwd,type="l")
lines(rangeJ, newdatM$plo[which(newdatM$Age=="J" & newdatM$MassCorected>min(xvJW) & newdatM$MassCorected<max(xvJW))],col=mycolors3[1],type="l",lty=2,lwd=4*faclwd)
lines(rangeJ,newdatM$phi[which(newdatM$Age=="J" & newdatM$MassCorected>min(xvJW) & newdatM$MassCorected<max(xvJW) )],col=mycolors3[1],type="l",lty=2,lwd=4*faclwd)
points(x = SurvWinter$MassCorected,y=SurvWinter$survival-0.01,col=mycolors3[1])

legend(x = 2,y=0.8,legend=c("Real mass","Asymptote $_{2006-2007}$","Asymptote $_{2008-2014}$"),fill = c(mycolors3[1],"dark red","red"),xjust = 0,yjust= 0.5,bty="n", x.intersp = 0.2)

lines(1/(1+exp(-newdatRef$logitphi)),x=newdatRef$A,col="red",lty=1,lwd=8*faclwd)
lines(1/(1+exp(-newdatRef$plo)),x=newdatRef$A,col="red",lty=2,lwd=4*faclwd)
lines(1/(1+exp(-newdatRef$phi)),x=newdatRef$A,col="red",lty=2,lwd=4*faclwd)
points(x = GrowthAKlphi2b$BUGSoutput$mean$A,y=phis+0.01,col="red")

lines(1/(1+exp(-newdatRefnoSel$logitphi)),x=newdatRefnoSel$A,col="dark red",lty=1,lwd=8*faclwd)
lines(1/(1+exp(-newdatRefnoSel$plo)),x=newdatRefnoSel$A,col="dark red",lty=2,lwd=4*faclwd)
lines(1/(1+exp(-newdatRefnoSel$phi)),x=newdatRefnoSel$A,col="dark red",lty=2,lwd=4*faclwd)
