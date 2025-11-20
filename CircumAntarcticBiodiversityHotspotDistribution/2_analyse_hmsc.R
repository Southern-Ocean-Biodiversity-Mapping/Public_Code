##############################################################################################################
library(Hmsc)
library(bayesplot)
library(corrplot)
library(coda)
##############################################################################################################

################################
##### load and setup data ######
################################

#### load data
# ## model runs
# load("...Rdata")
# ## model fit and predicted values
# load("_Mf.Rdata")
# ## crossvalidation
# load("..._5foldcv.Rdata")
# ## predicted vs observed
# load("..._PredVsObs.Rdata")

## reorder species to alphabetical
sp.v <- order(models[[2]]$spNames)

## remove first bit of the species labels
colnames(models[[1]]$Y) <- gsub('Biota - ', '', colnames(models[[1]]$Y))
colnames(models[[2]]$Y) <- gsub('Biota - ', '', colnames(models[[2]]$Y))
spNames.ab <- "totalabundance"
spNames.pa <- gsub('Biota - ', '', models[[2]]$spNames)

## extracting the posterior distribution from the model object
mpost <- list()
for(i in 1:2){
  mpost[[i]] <- convertToCodaObject(models[[i]])
}

################################
##### MCMC convergence #########
################################
#### graphical overview
## abundance
mcmc.dat <- mpost[[1]]$Beta[,1:15]
for(k in 1:4){colnames(mcmc.dat[[k]]) <- colnames(models[[1]]$X)}
plot(mcmc_trace(mcmc.dat))
## presence-absence
for(i in sp.v){
  v <- (1:15)-15+15*i
  cat(models[[2]]$spNames[i])
  mcmc.dat <- mpost[[2]]$Beta[,v]
  for(k in 1:4){colnames(mcmc.dat[[k]]) <- colnames(models[[2]]$X)}
  plot(mcmc_trace(mcmc.dat))
}
### If good, the following is true:
### - different chain yielding the same results
### - chains rise and fall rapidly without apparent autocorrelation
### - the first half looks essentially identical to the second half

#### quantitative overview
### effective size - abundance
es.dat <- effectiveSize(mpost[[1]]$Beta)
es.dat.ab <- matrix(es.dat, ncol=15, byrow=TRUE)
colnames(es.dat.ab) <- colnames(models[[1]]$X)
rownames(es.dat.ab) <- spNames.ab
kable(es.dat.ab)
### effective size - presence-absence
es.dat <- effectiveSize(mpost[[2]]$Beta)
es.dat.pa <- matrix(es.dat, ncol=15, byrow=TRUE)
colnames(es.dat.pa) <- colnames(models[[2]]$X)
rownames(es.dat.pa) <- spNames.pa
kable(es.dat.pa[sp.v,]) ## ordered by names
## If good, the following is true:  
## - effective sample size not too far away from sample size  
## - potential scale reduction factors close to 1, indicating the multiple chains give consistent results  

##### potential scale reduction factor
## the closer to 1 the better
gd.dat.ab <- gelman.diag(mpost[[1]]$Beta)$psrf
gd.dat.ab.pe <- matrix(gd.dat.ab[,1], ncol=15, byrow=TRUE)
gd.dat.ab.uci <- matrix(gd.dat.ab[,2], ncol=15, byrow=TRUE)
colnames(gd.dat.ab.pe) <- colnames(models[[1]]$X)
rownames(gd.dat.ab.pe) <- spNames.ab
colnames(gd.dat.ab.uci) <- colnames(models[[1]]$X)
rownames(gd.dat.ab.uci) <- spNames.ab

##### potential scale reduction factor {.tabset}
## the closer to 1 the better
gd.dat.pa <- gelman.diag(mpost[[2]]$Beta)$psrf
gd.dat.pa.pe <- matrix(gd.dat.pa[,1], ncol=15, byrow=TRUE)
gd.dat.pa.uci <- matrix(gd.dat.pa[,2], ncol=15, byrow=TRUE)
colnames(gd.dat.pa.pe) <- colnames(models[[2]]$X)
rownames(gd.dat.pa.pe) <- spNames.pa
colnames(gd.dat.pa.uci) <- colnames(models[[2]]$X)
rownames(gd.dat.pa.uci) <- spNames.pa

################################
##### parameter estimates ######
################################
#### parameter significance
postBeta.ab = getPostEstimate(models[[1]], parName="Beta")
postBeta.pa = getPostEstimate(models[[2]], parName="Beta")

source("https://raw.githubusercontent.com/Southern-Ocean-Biodiversity-Mapping/Public_Code/refs/heads/main/Hmsc_plotBetaSimple.R")
## response directions
plotBetaSimple(models[[1]], post=postBeta.ab, mar=c(8, 30, 2, 1), param="Support", spSort="reverse", main="Response direction")
plotBetaSimple(models[[2]], post=postBeta.pa, mar=c(8, 30, 2, 1), param="Support", spSort="reverse", main="Response direction")

## mean response
plotBetaSimple(models[[1]], post=postBeta.ab, mar=c(8, 30, 2, 1), param="Mean", spSort="reverse", main="Mean response estimate")
plotBetaSimple(models[[2]], post=postBeta.pa, mar=c(8, 30, 2, 1), param="Mean", spSort="reverse", main="Mean response estimate")

#### species associations
OmegaCor = computeAssociations(models[[2]])
supportLevel = 0.95
for (r in 1:models[[2]]$nr){
  plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
  toPlot = ((OmegaCor[[r]]$support>supportLevel) +
              (OmegaCor[[r]]$support<(1-supportLevel))>0)*OmegaCor[[r]]$mean
  par(xpd=T)
  colnames(toPlot)=rownames(toPlot)=gsub("_"," ",x=colnames(toPlot))
  corrplot(toPlot[plotOrder,plotOrder], method = "color",
           col=colorRampPalette(c("blue","white","red"))(200),
           title="",type="lower",tl.col="black",tl.cex=.4, mar=c(0,0,6,0))
}


################################
##### explanatory power ########
################################
#### Model Fit
# From the Hmsc vignette:  
#   All measures of model fit are based on comparing the posterior predictive distribution (predY)) to the observed values (hM$Y). The predicted distribution is first summarized to a single matrix of predicted values by taking the posterior mean (for normal and probit models) or posterior median (for Poisson models). All measures of model fit are given as vectors with one value for each species.  
# 
# The kinds of measures of model fit depend on the type of response variable.  
# - For all types of response variables, root-mean-square error (RMSE) between predicted and observed values is computed.  
# - For normal models, R2 is computed as squared pearson correlation between observed and predicted values, times the sign of the correlation.  
# - For probit models, Tjur R2 and AUC are computed.  
# - For Poisson models, a pseudo-R2 is computed as squared spearman correlation between observed and predicted values, times the sign of the correlation (SR2).  
# - For Poisson models, the observed and predicted data are also truncated to occurrences (presence-absences), for which the same measures are given as for the probit models (O.RMSE, O.AUC and O.TjurR2).  
# - For Poisson models, the observed and predicted data are also subsetted to conditional on presence, for which the root-mean-square error and pseudo-R2 based on squared spearman correlation are computed (C.RMSE, C.SR2).  
# 
# The measures O.RMSE, O.AUC, O.TjurR2, C.RMSE and C.SR2 can be computed only if the option expected=FALSE has been used when making the predictions If the model includes a mixture of response variable types, the resulting measures of model fit contain NA’s for those response variables for which they cannot be computed.  

##### RMSE - ab
MF[[1]]$RMSE
MF.cv[[1]]$RMSE

##### SR2 - ab
MF[[1]]$SR2
MF.cv[[1]]$SR2

par(mfrow=c(2,2), mar=c(5,4,4,1))
##### RMSE - pa
sel <- order(MF[[2]]$RMSE)
## model vs crossvalidation
plot(MF[[2]]$RMSE[sel], pch=16, main="presence-absence model vs Crossvalidation", ylim=c(0.1,0.5), ylab="RMSE")
points(MF.cv[[2]]$RMSE[sel], pch=16, col="red")
legend("topleft", pch=16, col=c("black","red"), legend=c(paste0("p-a model (mean RMSE = ",round(mean(MF[[2]]$RMSE),2),")"), paste0("CV             (mean RMSE = ",round(mean(MF.cv[[2]]$RMSE),2),")")))

##### AUC - pa
sel <- order(MF[[2]]$AUC)
## model vs crossvalidation
plot(MF[[2]]$AUC[sel], pch=16, main="presence-absence model vs Crossvalidation", ylim=c(0.1,1), ylab="AUC")
points(MF.cv[[2]]$AUC[sel], pch=16, col="red")
legend("bottomright", pch=16, col=c("black","red"), legend=c(paste0("p-a model (mean AUC = ",round(mean(MF[[2]]$AUC),2),")"), paste0("CV             (mean AUC = ",round(mean(MF.cv[[2]]$AUC),2),")")))

##### TjurR2 - pa
# Tjur R2 compares the average fitted probability of the two response outcomes. The difference between the average fitted probability for success and the average fitted probability for the failure.
sel <- order(MF[[2]]$TjurR2)
## model vs crossvalidation
plot(MF[[2]]$TjurR2[sel], pch=16, main="presence-absence model vs Crossvalidation", ylim=c(0,1), ylab="Tjur R^2")
points(MF.cv[[2]]$TjurR2[sel], pch=16, col="red")
legend("topleft", pch=16, col=c("black","red"), legend=c(paste0("p-a model (mean TjurR2 = ",round(mean(MF[[2]]$TjurR2),2),")"), paste0("CV             (mean TjurR2 = ",round(mean(MF.cv[[2]]$TjurR2),2),")")))

#### predicted vs observed
preds.cv.ab.median = as.data.frame(apply(preds.cv[[1]], 1:2, median))
preds.cv.pa.median = as.data.frame(apply(preds.cv[[2]], 1:2, median))

pred.ab    <- as.vector(predY.ab.median)[[1]]
pred.ab.cv <- as.vector(preds.cv.ab.median)[[1]]
obs.ab     <- as.vector(models[[1]]$Y)#/5430*540

pred.pa    <- rowSums(predY.pa.median)
pred.pa.cv <- rowSums(preds.cv.pa.median)
obs.pa     <- rowSums(models[[2]]$Y)

par(mfrow=c(2,2), mar=c(5,4,4,1))
## abundance
plot(pred.ab, obs.ab, pch=16, cex=0.2, main="Abundance model",
     xlab="predicted", ylab="observed", xlim=c(0,1600), ylim=c(0,1600))
abline(0,1, lty=2, col="grey")
legend("right", legend=paste0("RMSE = ", round(mean(MF[[1]]$RMSE),2),
                              "\nSR2 = ",round(mean(MF[[1]]$SR2),2)), bty="n", cex=0.5)
plot(pred.ab.cv, obs.ab, pch=16, cex=0.2, main="Abundance model - CV",
     xlab="predicted", ylab="observed", xlim=c(0,1600), ylim=c(0,1600))
abline(0,1, lty=2, col="grey")
legend("right", legend=paste0("RMSE = ", round(mean(MF.cv[[1]]$RMSE),2),
                              "\nSR2 = ",round(mean(MF.cv[[1]]$SR2),2)), bty="n", cex=0.5)
## presence-absence
plot(pred.pa, obs.pa, pch=16, cex=0.3, main="p-a model",
     xlab="predicted", ylab="observed", xlim=c(0,45), ylim=c(0,45))
abline(0,1, lty=2, col="grey")
legend("right", legend=paste0("RMSE = ", round(mean(MF[[2]]$RMSE),2),
                              "\nAUC = ",round(mean(MF[[2]]$AUC),2),
                              "\nTjurR2 = ",round(mean(MF[[2]]$TjurR2),2)), bty="n", cex=0.5)
plot(pred.pa.cv, obs.pa, pch=16, cex=0.3, main="p-a model - CV",
     xlab="predicted", ylab="observed", xlim=c(0,45), ylim=c(0,45))
abline(0,1, lty=2, col="grey")
legend("right", legend=paste0("RMSE = ",round(mean(MF.cv[[2]]$RMSE),2),
                              "\nAUC = ",round(mean(MF.cv[[2]]$AUC),2),
                              "\nTjurR2 = ",round(mean(MF.cv[[2]]$TjurR2),2)), bty="n", cex=0.5)


################################
##### variance explained #######
################################
## what proportion of variance is explained by the models and their parameters
groupnames = c("bathymetry", "waterproperties","food","effort")
group = c(1,1,1,1,1,1,2,2,2,2,3,3,3,4)
VP <- list()
for(i in 1:2){
  VP[[i]] = computeVariancePartitioning(models[[i]],group=group, groupnames=groupnames)
}

#### presence-absence
sel <- order(colnames(VP[[2]]$vals))
dat <- VP[[2]]
hM <- models[[2]]
ng = dim(dat$vals)[1]
colnames(dat$vals) <- spNames.pa
##
leg = dat$groupnames
for (r in 1:hM$nr) {
  leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
}
means = round(100 * rowMeans(dat$vals), 1)
for (i in 1:ng) {
  leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                 ")", sep = "")
}
##
par(mfrow=c(1,4), mar=c(10,2,1,1), oma=c(0,15,0,0))
plot.new()
barplot(dat$vals[,rev(sel)], horiz=TRUE, col=heat.colors(ng), las=1, legend=leg, main="variance partitioning", args.legend=list(x=1, y=-6, bty = "n"))
barplot(MF[[2]]$AUC[rev(sel)], horiz=TRUE, las=1, xlim=c(0,1), main="AUC Model vs CV", legend=c("Model","CV"), args.legend=list(x=0.6, y=-6, bty = "n", fill=c("black","grey50")))
barplot(MF.cv[[2]]$AUC[rev(sel)], horiz=TRUE, add=TRUE, col="grey50")
abline(v=0.2, lty=2)
abline(v=0.4, lty=2)
abline(v=0.6, lty=2)
abline(v=0.8, lty=2)
abline(v=1, lty=2)
barplot(MF[[2]]$TjurR2[rev(sel)], horiz=TRUE, las=1, xlim=c(0,1), main="TjurR2 Model vs CV")
barplot(MF.cv[[2]]$TjurR2[rev(sel)], horiz=TRUE, add=TRUE, col="grey50")
abline(v=0.2, lty=2)
abline(v=0.4, lty=2)
abline(v=0.6, lty=2)
abline(v=0.8, lty=2)
abline(v=1, lty=2)

#### abundance
dat <- VP[[1]]
hM <- models[[1]]
ng = dim(dat$vals)[1]
colnames(dat$vals) <- spNames.ab
##
leg = dat$groupnames
for (r in 1:hM$nr) {
  leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
}
means = round(100 * rowMeans(dat$vals), 1)
for (i in 1:ng) {
  leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                 ")", sep = "")
}
##
par(mfrow=c(1,4), mar=c(10,2,1,1), oma=c(0,15,0,0))
plot.new()
barplot(dat$vals, horiz=TRUE, col=heat.colors(ng+1), las=1, legend=leg, main="variance partitioning", args.legend=list(x=1.1, y=-0.4, bty = "n"))
barplot(MF[[1]]$SR2, horiz=TRUE, las=1, xlim=c(0,1), main="SR2 vs CV", args.legend=list(x=1.1, y=-0.4, bty = "n", fill=c("black","grey50")))
barplot(MF.cv[[1]]$SR2, horiz=TRUE, add=TRUE, col="grey50")
abline(v=0.2, lty=2)
abline(v=0.4, lty=2)
abline(v=0.6, lty=2)
abline(v=0.8, lty=2)
abline(v=1, lty=2)
barplot(MF.cv[[1]]$RMSE, horiz=TRUE, col="grey50", las=1, xlim=c(0,200), main="RMSE vs CV", args.legend=list(x=1, y=0, bty = "n", fill=c("black","grey50")))
barplot(MF[[1]]$RMSE, horiz=TRUE, add=TRUE)
abline(v=0.2, lty=2)
abline(v=0.4, lty=2)
abline(v=0.6, lty=2)
abline(v=0.8, lty=2)
abline(v=1, lty=2)


