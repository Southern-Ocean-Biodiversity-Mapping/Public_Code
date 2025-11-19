## fitting an hmsc using Otsos book, course scripts and : https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13345&file=mee313345-sup-0002-AppendixS2.pdf

##############################################################################################################


library(Hmsc)
library(terra)
'%!in%' <- function(x,y)!('%in%'(x,y))

##############################################################################################################

#res <- "500m"
res <- "2km"

##############################################################################################################
env.derived <- "/pvol3TB/data_environmental/"

## load scaled environmental rasters:
env_stack_scaled <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_scaled.tif"))

## load data
load(file=paste0("Cell_level_env_",res,"_202412.Rdata"))
load(file=paste0("Cell_level_bio_2pc_",res,"_202501.Rdata"))

# r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
# r2 <- r.stack$depth

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

## presence absence data:
cov_pa.raw <- cover_cells[,-1]
cov_pa.raw[cov_pa.raw>0] <- 1

## remove rare species
cov_pa.raw2 <- cov_pa.raw[,-which(colSums(cov_pa.raw)<=9)]

## combine UBS_B with Bryozoan_Hard_Branching_Antler
# cov_pa.raw2$Bryozoan_Hard_Branching_Antler <- cov_pa.raw2$Bryozoan_Hard_Branching_Antler+cov_pa.raw2$UBS_B
# cov_pa.raw2$Bryozoan_Hard_Branching_Antler[cov_pa.raw2$Bryozoan_Hard_Branching_Antler>1] <- 1
# cov_pa.raw3 <- cov_pa.raw2[,-grep("UBS_B",names(cov_pa.raw2))]
cov_pa.raw2$'Biota - Bryozoa - Hard - Branching - Antler' <- cov_pa.raw2$'Biota - Bryozoa - Hard - Branching - Antler'+cov_pa.raw2$'Biota - Matrix - UBS_B TemporaryLabelDoNotUse'
cov_pa.raw2$'Biota - Bryozoa - Hard - Branching - Antler'[cov_pa.raw2$'Biota - Bryozoa - Hard - Branching - Antler'>1] <- 1
cov_pa.raw3 <- cov_pa.raw2[,-grep("Biota - Matrix - UBS_B TemporaryLabelDoNotUse",names(cov_pa.raw2))]

## remove substrates, noid and unscorable
if(res=="2km") cov_pa <- cov_pa.raw3[,-c(grep("Physical",names(cov_pa.raw3)),
                                         grep("Unsco",names(cov_pa.raw3)),
                                         grep("General Unkn",names(cov_pa.raw3))[1:2])]

#save(cover_cells, cov_pa, metadat, file=paste0("2km_model_cells_data.Rdata"))
## combine UBS with Hydroid_Matrix?


###########################
##### set up the data #####

# ## for simplicity, start by analysing only 9 species at 200 sites
# s <- sample(1:nrow(cell_metadata_env_clean_scaled),200)
# s2 <- s[order(s)]
# metadat <- metadat[s2,]
# metadat$cellID <- factor(metadat$cellID)
# Y <- dat_cov_pa[s2,c(4,6,7:11,13,15)]  ## species data

## only Bryozoans
# Y <- dat_cov_pa[,c(4,5,13,16,17,20,21,30,31,74,61,97,104,105,107)] ##
#Y <- dat_cov_pa[,4]

## or go with the full dataset here:
# Y <- dat_cov_pa  ## species data

## COMPARE THIS LIST TO THE Species list EXCEL FILE WITH THE 2% CUTOFF!
# colSums(dat_cov_pa)[which((colSums(dat_cov_pa)/961)<0.018)]
#Y <- dat_cov_pa[,-which((colSums(dat_cov_pa)/nrow(dat_cov_pa))<0.018)]
Y <- cov_pa

## XData only the variables we choose in XFormula below
model_vars <- c("depth","depth2","logslope","tpi","distance2canyons","distance2canyons2",
                "seafloortemperature","seafloorcurrents_mean","seafloorcurrents_residual","seafloorsalinity","npp_mean","log.flux.mean","sed.mean","cover_points_scorable")
#XData <- dplyr::select(metadat, model_vars)
XData <- metadat[,which(names(metadat)%in%model_vars)]

############################
##### set up the model #####

## study design - a random spatial effect at the sample level (raster-cell), and a random effect at the survey level to address year and gear
#studyDesign <- data.frame(cellID=metadat$cellID)
studyDesign <- data.frame(cellID=metadat$cellID, surveyID=metadat$cover_cells_survey, transectID=metadat$cover_cells_transect1,
                          gear=metadat$gear, year=as.factor(metadat$year))

## survey effect:
rL.s = HmscRandomLevel(units=levels(metadat$cover_cells_survey))
## transect effect:
rL.t = HmscRandomLevel(units=levels(metadat$cover_cells_transect1))
## gear effect:
rL.g = HmscRandomLevel(units=levels(metadat$gear))
## year effect:
rL.y = HmscRandomLevel(units=levels(as.factor(metadat$year)))

## spatial random effect
xy <- metadat[,4:5]
colnames(xy) = c("x","y")
sRL = xy
rownames(sRL) = metadat$cellID

### using the standard algorithm:
## 5min per iteration to fit the spatial model, which is way too long (~1 month for 10k iterations)
#rL = HmscRandomLevel(sData=sRL)

### trying NNGP:
## doesn't work, the error message is: "Failed updaters and their counts in chain 1  ( 15  attempts)"
# rL = HmscRandomLevel(sData=sRL, sMethod="NNPG")
# rL = setPriors(rL,nfMin=1,nfMax=1)

### trying GPP:
## ~ 40s for 10 iterations; 60s for 50 iterations, 84s for 100 iterations; using knots at 500km distance -> ~1h20min for 10k iterations
## ~ 3.5min for 10 iterations; 5min for 50 iterations; 6.6min for 100 iterations; using knots at 250km distance -> ~2h40min for 10k iterations
## ~ 7min for 2 iterations; 15min for 100 iterations using knots at 200km distance -> ~13h20min for 10k iterations
## BUT, on a 250km grid, with the full dataset, the predictions will take 41 days!!!
## first specifying knots on a grid
# xy.knots <- rbind(xy,c(2900000,0)) ## add a point to the right to allow mapping of East Antarctica
# Knots = constructKnots(xy.knots, knotDist = 250000, minKnotDist = 2500000)
# #Knots = constructKnots(xy, knotDist = 50000, minKnotDist = 2000000)
# plot(xy.knots[,1],xy.knots[,2],pch=18, asp=1)
# points(Knots[,1],Knots[,2],col='red',pch=18)

## knots at 200km distance, 250km min distance:
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
rL.s$nfMax=10
rL.t$nfMax=10
rL.g$nfMax=10
rL.y$nfMax=10

###
XFormula = ~depth+depth2+logslope+tpi+distance2canyons+distance2canyons2+seafloortemperature+seafloorcurrents_mean+seafloorcurrents_residual+seafloorsalinity+npp_mean+log.flux.mean+sed.mean+log(cover_points_scorable)

mFULL = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
             studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))
mENV = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
            studyDesign = studyDesign, ranLevels = list(surveyID=rL.s))
mSPACE = Hmsc(Y = Y, XData = XData, XFormula = ~1, distr = "probit",
              studyDesign = studyDesign, ranLevels = list(cellID=rL, surveyID=rL.s))

#######################################
##### run MCMC and save the model #####
#######################################
## 5 days on the VM to fit all three 2km models with 8000 iterations and 83 species

models <- list(mFULL, mENV, mSPACE)
modeltype = 1
model = 1
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 800 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
nChains = 4
filename.string <- paste(res,"_model_cells_", as.character(model), "_",
                         c("pa","abundance")[modeltype], 
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
filename = file.path(paste("/pvol3TB/2_fitting_and_running_models/",filename.string, ".Rdata", sep = ""))
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
# ## full model
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
filename2 = file.path(paste("/pvol3TB/3_model_analysis/",filename.string, "_MF.Rdata", sep = ""))
save(MF, MF.preds, file=filename2)


####################################################### TAKES X-TIMES (FOLDS) LONGER THAN THE MODEL FITTING!!!
##### evaluating model fit using cross validation #####
load("/pvol3TB/2_fitting_and_running_models/2km_model_cells_1_pa_chains_4_thin_10_samples_800.Rdata")
set.seed(2)
partition = createPartition(models[[2]], nfolds=5, column="transectID") ## use column to partition according to different hierarchies (e.g. leave an entire region out)
# partition = createPartition(models[[1]], nfolds=5, column="surveyID") ## use column to partition according to different hierarchies (e.g. leave an entire region out)
# partition = createPartition(models[[1]], nfolds=10, column="surveyID") ## use column to partition according to different hierarchies (e.g. leave an entire region out)
# preds.cv = computePredictedValues(m, partition = partition)
# MF.cv = evaluateModelFit(hM=m, predY = preds.cv)
# MF.cv

## This takes a long time, 5 days for the full model & 1h for the environment only model
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
for(i in 2:1){
  print(i)
  preds.cv[[i]] = pcomputePredictedValues(models[[i]], partition=partition, nParallel=32)
  MF.cv[[i]] = evaluateModelFit(hM=models[[i]], predY = preds.cv[[i]])
}
computational.time = proc.time() - ptm
parallel::stopCluster(cl = c1)
filename3 = file.path(paste("/pvol3TB/3_model_analysis/",filename.string, "_5foldcv.Rdata", sep = ""))
#filename3 = file.path(paste(filename.string, "_10foldcv.Rdata", sep = ""))
## NOTE THAT PREVIOUSLY I DIDN'T SAVE THE PREDS.CV FILE!!!
save(MF.cv, preds.cv, partition, computational.time, file=filename3)
# for(i in 2){
#   preds = computePredictedValues(models[[i]], partition=partition)
#   MF.cv[[i]] = evaluateModelFit(hM=models[[i]], predY = preds)
# }
# save(MF.cv, file=filename2, computational.time)














