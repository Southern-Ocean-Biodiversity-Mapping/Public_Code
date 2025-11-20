## fitting an hmsc using Otsos book, course scripts and : https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13345&file=mee313345-sup-0002-AppendixS2.pdf


## this is to fit a single total abundance model
## or a functional groups hmsc with only a few groups
##############################################################################################################


library(Hmsc)
library(terra)
library(dplyr)
'%!in%' <- function(x,y)!('%in%'(x,y))

##############################################################################################################

#res <- "500m"
res <- "2km"

##############################################################################################################
## load data
load(file=paste0("~/Cell_level_env_",res,"_202412.Rdata"))
load(file=paste0("~/Cell_level_bio_2pc_",res,"_202501.Rdata"))

## check for NAs in the data:
## find NAs in waom and npp data
waom.na.sel <- which(is.na(cell_metadata_env_scaled$seafloorcurrents_absolute))
npp.na.sel <- which(is.na(cell_metadata_env_scaled$npp_mean))
## remove seamount transects (tan1802 & tan1901)
seamount.transects <- c("TAN1802_160","TAN1802_170","TAN1802_179","TAN1802_180","TAN1802_184","TAN1802_185","TAN1802_191","TAN1802_193",
                        "TAN1802_195","TAN1802_196","TAN1802_197","TAN1802_207","TAN1802_208","TAN1802_209","TAN1802_213","tan1901_209")
seamount.transects.sel <- which(cell_metadata_env_scaled$cover_cells_transect1%in%seamount.transects)
## combine to one vector
na.sel <- unique(c(waom.na.sel,npp.na.sel,seamount.transects.sel))

## remove nas from data
cover_cells <- cover_mod.2km[-na.sel,]
metadat <- cell_metadata_env_scaled[-na.sel,]
metadat$cellID <- factor(metadat$cellID)

## remove substrates, noid and unscorable
cov_ab.raw <- cover_cells[,-1]
if(res=="2km") cov_ab.raw <- cov_ab.raw[,-c(grep("Physical",names(cov_ab.raw)),
                                         grep("Unsco",names(cov_ab.raw)),
                                         grep("General Unkn",names(cov_ab.raw))[1:2])]
# cov_ab.raw <- cov_ab.raw[,-c(grep("Sub",names(cov_ab.raw)),
#                 grep("Unsco",names(cov_ab.raw)),
#                 grep("NoID",names(cov_ab.raw))[1:2])]

## total cover
cov_ab <- rowSums(cov_ab.raw)

# ## FG cover (if I run this, NEED TO FIX WHICH NUMBERS BELONG TO WHICH CATEGORY)
# sel.SF <- c(1,3:9,12,14:35,37:38,40:41,45:52,54,57:59,65,67,69:74,76,78,80:84,86:87,89:90,92,96:97)
# sel.DF <- c(13,36,39,56)
# sel.Oth <- c(2,10,11,42:44,53,55,60:64,66,68,75,77,79,85,88,91,93:95)
# cov_ab.1 <- rowSums(cov_ab.raw[,sel.SF]) ## Suspension feeders
# cov_ab.2 <- rowSums(cov_ab.raw[,sel.DF]) ## Deposit feeders
# cov_ab.3 <- rowSums(cov_ab.raw[,sel.Oth]) ## Other
# cov_ab <- cbind(cov_ab.1, cov_ab.2, cov_ab.3)
# names(cov_ab) <- c("SF", "DF", "Oth")

## divide by number of scorable points per cell
cov_ab_perc <- (cov_ab/metadat$cover_points_scorable)*100

cov_ab_perc_scaledlog <- scale(log(cov_ab_perc))
hist(cov_ab_perc_scaledlog)
# par(mfrow=c(4,4))
# for(i in 1:ncol(cov_ab_perc_scaledlog)){
#   hist(cov_ab_perc_scaledlog[,i], main=colnames(cov_ab_perc_scaledlog)[i])
# }
hist(as.matrix(cov_ab))


###########################
##### set up the data #####
## Spong erect simple (#80) is a problem when choosing initpar="fixed effects" in sampleMcmc later
Y <- cov_ab#_perc_scaledlog#[,-80]

## XData only the variables we choose in XFormula below
model_vars <- c("depth","depth2","logslope","tpi","distance2canyons","distance2canyons2",
                "seafloortemperature","seafloorcurrents_mean","seafloorcurrents_residual","seafloorsalinity","npp_mean","log.flux.mean","sed.mean","cover_points_scorable")
#XData <- dplyr::select(metadat, model_vars)
XData <- metadat[,which(names(metadat)%in%model_vars)]
# XData$offset <- log(metadat$cover_points_scorable)
  
############################
##### set up the model #####
## study design
studyDesign <- data.frame(cellID=metadat$cellID, surveyID=metadat$cover_cells_survey)
# studyDesign <- data.frame(cellID=metadat$cellID, surveyID=metadat$cover_cells_survey, transectID=metadat$cover_cells_transect1)

## spatial random effect
xy <- metadat[,4:5]
colnames(xy) = c("x","y")
sRL = xy
rownames(sRL) = metadat$cellID

########## knots at 200km distance, 250km min distance:
## add points between AP and Ross Sea
xy.knots <- rbind(xy,
                  c(-2143647,498436),
                  c(-1916289,355285),
                  c(-2000000,100000),
                  c(-1983655,-368890),
                  c(-1739456,-638350),
                  c(-1621567,-975176),
                  c(-1360527,-1160431),
                  c(-900000,-1300000),
                  c(-627930,-1345685))
## add points to East Antarctica
xy.knots <- rbind(xy.knots,
                  c(710952,-2154067),
                  c(1115144,-2288798),
                  c(2100359,-1724614),
                  c(2529812,-1017280),
                  c(2757170,-629930),
                  c(2824535,-318366),
                  c(2672963,-57326),
                  c(2656122,254237),
                  c(2487709,498436),
                  c(2386661,742635),
                  c(2268772,1256295),
                  c(2125621,1685748),
                  c(1645644,1812057),
                  c(820421,2056256),
                  c(1200000,2000000))
## add a point to Weddell Sea
xy.knots <- rbind(xy.knots,
                  c(-1503678,1054199),
                  c(-1436312,1300000),
                  c(-1958393,1298398))
Knots = constructKnots(xy.knots, knotDist = 200000, minKnotDist = 250000)
rL = HmscRandomLevel(sData=sRL, sMethod='GPP', sKnot=Knots)
rL$nfMax=10

## survey effect:
rL.s = HmscRandomLevel(units=levels(metadat$cover_cells_survey))
rL.s$nfMax=10

# ## transect effect: 322 transects which might be a bit much. Also, hmsc gives an error
# rL.t = HmscRandomLevel(units=levels(metadat$cover_cells_transect1))
# rL.t$nfMax=10

##
# XFormula = ~depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+seafloortemperature+seafloorcurrents_mean+seafloorsalinity+npp_mean+cover_points_scorable
XFormula = ~depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+seafloortemperature+seafloorcurrents_mean+seafloorcurrents_residual+seafloorsalinity+npp_mean+log.flux.mean+sed.mean+log(cover_points_scorable)

mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
             studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "lognormal poisson",
            studyDesign = studyDesign, ranLevels = list(surveyID=rL.s))
mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "lognormal poisson",
              studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))

# mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "normal",
#              studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
# mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "normal",
#             studyDesign = studyDesign, ranLevels = list(surveyID=rL.s))
# mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "normal",
#               studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))


#######################################
##### run MCMC and save the model #####
#######################################
## 5 days on the VM to fit all three 2km models with 8000 iterations and 75 species

models <- list(mFULL, mENV, mSPACE)
model = 1
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 800 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
nChains = 4
filename.string <- paste(res,"_model_cells_", as.character(model), "_",
                          "totalabundance", 
                          "_chains_",as.character(nChains),
                          "_thin_", ... = as.character(thin),
                          "_samples_", as.character(samples), sep = "")

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
computational.time = proc.time() - ptm
filename = file.path(paste(filename.string, ".Rdata", sep = ""))
save(models, file=filename, computational.time)

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
## full model
MF.preds$fm.mean   = as.data.frame(apply(preds[[1]], 1:2, mean))
MF.preds$fm.median = as.data.frame(apply(preds[[1]], 1:2, median))
MF.preds$fm.sd     = as.data.frame(apply(preds[[1]], 1:2, sd))
MF.preds$fm.5      = as.data.frame(apply(preds[[1]], 1:2, p5))
MF.preds$fm.95     = as.data.frame(apply(preds[[1]], 1:2, p95))
## env-only
MF.preds$eo.mean   = as.data.frame(apply(preds[[2]], 1:2, mean))
MF.preds$eo.median = as.data.frame(apply(preds[[2]], 1:2, median))
MF.preds$eo.sd     = as.data.frame(apply(preds[[2]], 1:2, sd))
MF.preds$eo.5      = as.data.frame(apply(preds[[2]], 1:2, p5))
MF.preds$eo.95     = as.data.frame(apply(preds[[2]], 1:2, p95))
## save files
filename2 = file.path(paste(filename.string, "_MF.Rdata", sep = ""))
save(MF, MF.preds, file=filename2)


####################################################### TAKES X-TIMES (FOLDS) LONGER THAN THE MODEL FITTING!!!
load(paste(filename.string, ".Rdata", sep = ""))
## add transect ID to model
models[[2]]$studyDesign$transectID <- metadat$cover_cells_transect1
##### evaluating model fit using cross validation #####
set.seed(2)
#partition = createPartition(models[[1]], nfolds=5, column="surveyID") ## use column to partition according to different hierarchies (e.g. leave an entire region out)
partition = createPartition(models[[2]], nfolds=5, column="transectID")
# preds.cv = computePredictedValues(m, partition = partition)
# MF.cv = evaluateModelFit(hM=m, predY = preds.cv)
# MF.cv

## This takes a long time, 10 days for the full model & 1h for the environment only model
## parallel processing: PER CELL that contains values
library(doParallel)
library(foreach)
parallel::detectCores()
#UseCores = parallel::detectCores() - 1
UseCores = 20
c1<-makeCluster(UseCores, outfile="", type="FORK") ## "FORK" is faster than "PSOCK", but only works on linux/mac
registerDoParallel(c1)
getDoParWorkers()

ptm = proc.time()
preds.cv = list()
MF.cv = list()
for(i in 1:2){
  preds.cv[[i]] = pcomputePredictedValues(models[[i]], partition=partition, nParallel=20)
  MF.cv[[i]] = evaluateModelFit(hM=models[[i]], predY = preds.cv[[i]])
}
computational.time = proc.time() - ptm
parallel::stopCluster(cl = c1)
filename2 = file.path(paste(res,"_model_cells_", as.character(model), "_",
                            "totalabundance", 
                            "_chains_",as.character(nChains),
                            "_thin_", ... = as.character(thin),
                            "_samples_", as.character(samples), "_5foldcv.Rdata", sep = ""))
save(MF.cv, preds.cv, partition, file=filename2, computational.time)
# for(i in 2){
#   preds = computePredictedValues(models[[i]], partition=partition)
#   MF.cv[[i]] = evaluateModelFit(hM=models[[i]], predY = preds)
# }
# save(MF.cv, file=filename2, computational.time)














