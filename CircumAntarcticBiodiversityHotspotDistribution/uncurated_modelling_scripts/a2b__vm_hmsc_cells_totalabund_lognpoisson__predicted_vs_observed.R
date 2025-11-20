## Calculating predicted abundance vs observed abundance at the sampling sites
## NOTE THAT SURVEY EFFECT IS DISREGARDED (CANNOT BE PREDICTED FOR...), so predictions are for a "typical" survey
## This takes about 1h to complete for the spatial model
## Using Gradient in the predictions changes the outcome in unexpected ways, not sure why. I have removed that step.

library(Hmsc)
library(terra)
'%!in%' <- function(x,y)!('%in%'(x,y))

p5 <- function(x, z=1.96) {
  mean(x) - z * sd(x) / sqrt(length(x))
}
p95 <- function(x, z=1.96) {
  mean(x) + z * sd(x) / sqrt(length(x))
}

##############################################################################################################
#res <- "500m"
res <- "2km"
##############################################################################################################
env.derived <- "/pvol3TB/data_environmental/"

dir <- "~/"
dir2 <- "/pvol3TB/2_fitting_and_running_models/"
dir3 <- "/pvol3TB/3_model_analysis/"

r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth
## create an empty raster to fill for mapping
empty.ra <- rast(r2)
empty.ra[] <- NA

## specify model to load
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 800 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
nChains = 4

## abundance model
model = 1
ab_string <- "_totalabundance"
#ab_string <- "_FGabundance"
#load(paste0("~/",res,"_model_cells_",model,ab_string,"_chains_4_thin_10_samples_800.Rdata"))
load(paste0(dir2,res,"_model_cells_",model,ab_string,"_chains_4_thin_10_samples_800.Rdata"))
ab <- models[[2]]

rm(models)

#############################################################
## we only need to predict on the sites we have observations from, and predict on the number of images/points we observed

## environmental data to predict on: 
XData.grid <- ab$XData

## locations for the random effects:
xy.grid <- ab$ranLevels$cellID[]$s

#############################################################
ptm = proc.time()
## setup environmental data
XData.grid.loop <- XData.grid
xy.grid.loop <- xy.grid

## setup prediction - abund
#GRADIENT IS DOING WEIRD STUFF TO THE PREDICTIONS: #Gradient.ab = prepareGradient(ab, XDataNew = XData.grid.loop, sDataNew = list(cellID = xy.grid.loop))
predY.loop.ab <- predict(ab, XData = XData.grid.loop, expected=TRUE)   #predict(ab, Gradient=Gradient.ab, expected=TRUE) ## this gives probabilities instead of integer outcomes
#rm(Gradient.ab)  

## Multiply matrices
predY.loop = simplify2array(predY.loop.ab)

## Get posterior median
predY.mean = as.data.frame(apply(predY.loop, 1:2, mean))
predY.median = as.data.frame(apply(predY.loop, 1:2, median))

## Posterior width
predY.5 = as.data.frame(apply(predY.loop, 1:2, p5))
predY.95 = as.data.frame(apply(predY.loop, 1:2, p95))

## save-string
dat.name <- paste0(dir3,res,"_model_cells_",model,ab_string,"_chains_4_thin_10_samples_800_PredVsObs.Rdata")
## save and then remove objects
save(predY.mean, predY.median, predY.5, predY.95,
     XData.grid.loop, xy.grid.loop,
     file=dat.name)

computational.time = proc.time() - ptm

#############################################################
## CV predictions



