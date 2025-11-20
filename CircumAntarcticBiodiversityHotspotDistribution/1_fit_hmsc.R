##############################################################################################################
library(Hmsc)
##############################################################################################################

################################
##### load and setup data ######
################################

#### load data
dat.bio <-    read.csv("https://raw.githubusercontent.com/Southern-Ocean-Biodiversity-Mapping/Public_Code/refs/heads/main/CircumAntarcticBiodiversityHotspotDistribution/data_biodiversity_cells_2km.csv")
dat.scaled <- read.csv("https://raw.githubusercontent.com/Southern-Ocean-Biodiversity-Mapping/Public_Code/refs/heads/main/CircumAntarcticBiodiversityHotspotDistribution/data_environmental_scaled_cells_2km.csv")

#### check and remove data-points for which environmental data is incomplete:
waom.na.sel <- which(is.na(dat.scaled$seafloorcurrents_mean))
npp.na.sel <- which(is.na(dat.scaled$npp_mean))
## combine to one vector
na.sel <- unique(c(waom.na.sel,npp.na.sel))
## remove nas from data
dat.bio <- dat.bio[-na.sel,]
dat.scaled <- dat.scaled[-na.sel,]
## change columns to factors where relevant
dat.scaled[c(1,8:12,14:16)] <- lapply(dat.scaled[c(1,8:12,14:16)], factor) 

#### presence-absence data
## first calculate presence-absence data from cover data
cov_pa.raw <- dat.bio[,-1]
cov_pa.raw[cov_pa.raw>0] <- 1
## remove rare species,we want to only analyse species that occur at least at 10 sites
cov_pa.raw2 <- cov_pa.raw[,-which(colSums(cov_pa.raw)<=9)]
## remove substrate categories, unidentified morphotypes and unscorable points
cov_pa <- cov_pa.raw2[,-c(grep("Physical",names(cov_pa.raw2)),
                          grep("Unsco",names(cov_pa.raw2)),
                          grep("General.Unkn",names(cov_pa.raw2))[1:2])]

#### total abundance data
## aggregate abundance of all morphospecies together
cov_ab.raw <- dat.bio[,-1]
## remove substrate categories, unidentified morphotypes and unscorable points
cov_ab.raw <- cov_ab.raw[,-c(grep("Physical",names(cov_ab.raw)),
                             grep("Unsco",names(cov_ab.raw)),
                             grep("General.Unkn",names(cov_ab.raw))[1:2])]
## total cover
cov_ab <- rowSums(cov_ab.raw)

rm(cov_ab.raw, cov_pa.raw, cov_pa.raw2, na.sel, npp.na.sel, waom.na.sel)

##############################################################################################################
#### setup model parameters
model_vars <- c("depth","depth2","logslope","tpi","distance2canyons","distance2canyons2",
                "seafloortemperature","seafloorcurrents_mean","seafloorcurrents_residual","seafloorsalinity","npp_mean","log.flux.mean","sed.mean","cover_points_scorable")
XData <- dat.scaled[,which(names(dat.scaled)%in%model_vars)]
## study design
studyDesign <- data.frame(cellID=dat.scaled$cellID, surveyID=dat.scaled$cover_cells_survey)#, transectID=dat.scaled$cover_cells_transect1)
## hmsc model formula
XFormula = ~depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+seafloortemperature+seafloorcurrents_mean+seafloorcurrents_residual+seafloorsalinity+npp_mean+log.flux.mean+sed.mean+log(cover_points_scorable)
## random survey effect
rL.s = HmscRandomLevel(units=levels(dat.scaled$cover_cells_survey))
rL.s$nfMax=10
## hmsc models
m.pa = Hmsc(Y = cov_pa, XData = XData, XFormula = XFormula, distr = "probit",
            studyDesign = studyDesign, ranLevels = list(surveyID=rL.s))
m.ab = Hmsc(Y = cov_ab, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
            studyDesign = studyDesign, ranLevels = list(surveyID=rL.s))

################################
##### HMSC run modelling #######
################################

#### run hmsc
models <- list(m.ab, m.pa)
thin = 10  ## keeping every 10th iteration (the higher the less correlated the samples are but the longer it takes)
samples = 800 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
nChains = 4
set.seed(2)
ptm = proc.time()
for(i in 1:2){
  print(i)
  print(proc.time())
  models[[i]] <- sampleMcmc(models[[i]], samples = samples, thin = thin, transient = transient,
                            nChains = nChains, nParallel = nChains,
                            #initPar = "fixed effects",
                            updater = list(GammaEta = FALSE))
}
(computational.time = proc.time() - ptm)
#save(models, file = "...Rdata")

################################
##### evaluating model fit #####
################################
MF <- list()
preds = list()
for(i in 1:2){
  preds[[i]] = computePredictedValues(models[[i]])
  MF[[i]] = evaluateModelFit(hM=models[[i]], predY = preds[[i]])
}
MF
## calculate median, mean, sd and CIs for each prediction:
p5  <- function(x, z=1.96) {mean(x) - z * sd(x) / sqrt(length(x))}
p95 <- function(x, z=1.96) {mean(x) + z * sd(x) / sqrt(length(x))}
MF.preds <- list()
## pa model
MF.preds$pa.mean   = as.data.frame(apply(preds[[1]], 1:2, mean))
MF.preds$pa.median = as.data.frame(apply(preds[[1]], 1:2, median))
MF.preds$pa.sd     = as.data.frame(apply(preds[[1]], 1:2, sd))
MF.preds$pa.5      = as.data.frame(apply(preds[[1]], 1:2, p5))
MF.preds$pa.95     = as.data.frame(apply(preds[[1]], 1:2, p95))
## ab model
MF.preds$ab.mean   = as.data.frame(apply(preds[[2]], 1:2, mean))
MF.preds$ab.median = as.data.frame(apply(preds[[2]], 1:2, median))
MF.preds$ab.sd     = as.data.frame(apply(preds[[2]], 1:2, sd))
MF.preds$ab.5      = as.data.frame(apply(preds[[2]], 1:2, p5))
MF.preds$ab.95     = as.data.frame(apply(preds[[2]], 1:2, p95))
## save files
#save(MF, MF.preds, file = "..._MF.Rdata")

################################
##### cross-validation #########
################################
models[[1]]$studyDesign$transectID=dat.scaled$cover_cells_transect1
models[[2]]$studyDesign$transectID=dat.scaled$cover_cells_transect1

set.seed(2)
partition.ab = createPartition(models[[1]], nfolds=5, column="transectID")
partition.pa = createPartition(models[[2]], nfolds=5, column="transectID")

## This takes 5-times longer than the standard model fit
library(doParallel)
library(foreach)
parallel::detectCores()
#UseCores = parallel::detectCores() - 1
UseCores = 32
c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
registerDoParallel(c1)
getDoParWorkers()

ptm = proc.time()
MF.cv = list()
preds.cv = list()
for(i in 1:2){
  print(i)
  if(i == 1) {partition=partition.ab
  }else partition=partition.pa 
  preds.cv[[i]] = pcomputePredictedValues(models[[i]], partition=partition, nParallel=32)
  MF.cv[[i]] = evaluateModelFit(hM=models[[i]], predY = preds.cv[[i]])
}
(computational.time = proc.time() - ptm)
parallel::stopCluster(cl = c1)

#save(MF.cv, preds.cv, file = "..._5foldcv.Rdata")

################################
##### predicted vs observed ####
################################
p5 <- function(x, z=1.96) {
  mean(x) - z * sd(x) / sqrt(length(x))
}
p95 <- function(x, z=1.96) {
  mean(x) + z * sd(x) / sqrt(length(x))
}
# ## load model runs
# load("...")
# ab <- models[[1]]
# pa <- models[[2]]
# rm(models)

## environmental data to predict on: 
XData.grid <- pa$XData
## locations for the random effects:
xy.grid <- pa$ranLevels$cellID[]$s

## setup prediction - pa
predY.pa <- predict(pa, XData=XData.grid, expected=TRUE)
#Gradient.pa = prepareGradient(pa, XDataNew = XData.grid, sDataNew = list(cellID = xy.grid))
#predY.pa <- predict(pa, Gradient=Gradient.pa, expected=TRUE) ## this gives probabilities instead of integer outcomes
## setup prediction - abund 
predY.ab <- predict(ab, XData=XData.grid, expected=TRUE)

## Multiply matrices
predY.pa.array = simplify2array(predY.pa)
predY.ab.array = simplify2array(predY.ab)

## Get posterior median
predY.pa.mean =   as.data.frame(apply(predY.pa.array, 1:2, mean))
predY.pa.median = as.data.frame(apply(predY.pa.array, 1:2, median))
predY.ab.mean =   as.data.frame(apply(predY.ab.array, 1:2, mean))
predY.ab.median = as.data.frame(apply(predY.ab.array, 1:2, median))

## Posterior width
predY.pa.5 =  as.data.frame(apply(predY.pa.array, 1:2, p5))
predY.pa.95 = as.data.frame(apply(predY.pa.array, 1:2, p95))
predY.ab.5 =  as.data.frame(apply(predY.ab.array, 1:2, p5))
predY.ab.95 = as.data.frame(apply(predY.ab.array, 1:2, p95))

## save
# save(predY.pa.mean, predY.pa.median, predY.pa.5, predY.pa.95,
#      predY.ab.mean, predY.ab.median, predY.ab.5, predY.ab.95, 
#      XData.grid, xy.grid,
#      file="..._PredVsObs.Rdata")

#############################################################



