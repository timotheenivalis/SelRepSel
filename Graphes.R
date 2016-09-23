
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