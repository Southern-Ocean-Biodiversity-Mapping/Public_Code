## assemble hmsc model predictions


## fitting an hmsc using Otsos book, course scripts and : https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13345&file=mee313345-sup-0002-AppendixS2.pdf
library(Hmsc)
library(terra)
library(bayesplot)
'%!in%' <- function(x,y)!('%in%'(x,y))

dir <- "~/"
dir2 <- "/pvol3TB/2_fitting_and_running_models/"
dir3 <- "/pvol3TB/3_model_analysis/"
dir4 <- "/pvol3TB/4_model_prediction/"

###############################
#res <- "500m"
res <- "2km"

## select which model to run
#model.sel <- 1 ## full model
model.sel <- 2 ## environment only model

if(model.sel==2){
  modelspec <- "_envonly"
} else modelspec <- ""

###############################

## specify model to load
thin = 10  ## a value of 10 means every 10th iteration is kept (the higher the less correlated the samples are but the longer it takes)
samples = 800 ## how many total samples we want
transient = ceiling(0.5*samples*thin)
nChains = 4

## presence absence model
modeltype = 1
model = 1
load(paste0(dir2,res,"_model_cells_",model,"_pa_chains_4_thin_10_samples_800.Rdata"))
pa <- models[[model.sel]]
load(paste0(dir3,res,"_model_cells_",model,"_pa_chains_4_thin_10_samples_800_MF.Rdata"))
pa.MF <- MF[[model.sel]]

## pa character string to save the files:
pa.base.str <- paste0("/pvol3TB/biodiversity_prediction/",res,
                   "_model_cells_", as.character(model), "_",
                   c("pa","abundance")[modeltype], 
                   "_chains_",as.character(nChains),
                   "_thin_", ... = as.character(thin),
                   "_samples_", as.character(samples),modelspec,"_")

rm(models, MF)

## abundance model
model = 1
load(paste0(dir2,res,"_model_cells_",model,"_abundcondpres_chains_4_thin_10_samples_800.Rdata"))
ab <- models[[model.sel]]
load(paste0(dir3,res,"_model_cells_",model,"_abundcondpres_chains_4_thin_10_samples_800_MF.Rdata"))
ab.MF <- MF[[model.sel]]

## ab character string to save the files:
ab.base.str <- paste0("/pvol3TB/biodiversity_prediction/",res,
                   "_model_cells_", as.character(model), "_",
                   "abundcondpres", 
                   "_chains_",as.character(nChains),
                   "_thin_", ... = as.character(thin),
                   "_samples_", as.character(samples),modelspec,"_")

rm(models, MF)

## reorder species to alphabetical
sp.v <- order(pa$spNames)

#### environmental stuff
env.derived <- "/pvol3TB/data_environmental/"
r.stack <- rast(paste0(env.derived,"Circumpolar_EnvData_",res,"_shelf_mask_unscaled_variables.tif"))
r2 <- r.stack$depth
empty.ra <- rast(r2)
empty.ra[] <- NA

#### load model stuff
# m <- models[[1]]
load(paste0("/pvol3TB/biodiversity_prediction/hmsc_",res,"_model_cell_sel.Rdata"))
load(paste0("/pvol3TB/biodiversity_prediction/hmsc_",res,"_model_cell_grid.Rdata"))

# load(file=paste0("/pvol/biodiversity_prediction/",res,"_model_50km_cells_with_data.Rdata"))
load(file=paste0("/pvol3TB/biodiversity_prediction/",res,"_model_100km_cells_with_data.Rdata"))

################

## identify prediction files that we need to load
pred.dir <- paste0("/pvol3TB/biodiversity_prediction/pred_files",modelspec,"/",res,"/abundance/")
pred.files.raw <- list.files(pred.dir)
pred.files <- pred.files.raw[grep("abundcondpres", pred.files.raw)]
# pred.files <- pred.files[-grep("denmanpreds.Rdata", pred.files)]
# pred.files <- pred.files[grep("800_pred",pred.files)]
# pred.files <- pred.files[grep("2_abundance",pred.files)]
# pred.files.pa <- pred.files[grep("pa_only",pred.files)]
pred.files.pa <- pred.files
pred.files.ab <- pred.files
# pred.files.ab <- pred.files[-grep("pa_only",pred.files)]

#########################################

###### PRESENCE/ABSENCE DATA
#### load all prediction files into single files, starting with the first file and then loop
#### THIS CAN BE DONE MUCH FASTER BY SPLITTING THE JOBS RATHER THAN ONE LONG RUN (3h)
load(paste0(pred.dir,pred.files.pa[1]))
## which cells do we need to fill with data
row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
all.sel.ra <- sel.not.na[sel[row.numbers]]
## fill objects for each value
all.predY.mean <- predY.pa.mean
all.predY.median <- predY.pa.median
all.predY.se <- predY.pa.se
all.predY.5 <- predY.pa.5
all.predY.95 <- predY.pa.95
## loop through all other species
for(i in 2:length(pred.files.pa)){
  print(i)
  load(paste0(pred.dir,pred.files.pa[i]))
  row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
  all.sel.ra <- c(all.sel.ra, sel.not.na[sel[row.numbers]])
  all.predY.mean <- rbind(all.predY.mean, predY.pa.mean)
  all.predY.median <- rbind(all.predY.median, predY.pa.median)
  all.predY.se <- rbind(all.predY.se, predY.pa.se)
  all.predY.5 <- rbind(all.predY.5, predY.pa.5)
  all.predY.95 <- rbind(all.predY.95, predY.pa.95)
}

## now fill raster cells with values from the prediction file, save mean and se output for each species
colnames(all.predY.mean) <- gsub("[/()]", "", colnames(all.predY.mean))
pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,1]
pred.ra[[2]][all.sel.ra]  <- all.predY.median[,1]
pred.ra[[3]][all.sel.ra]  <- all.predY.se[,1]
pred.ra[[4]][all.sel.ra]  <- all.predY.5[,1]
pred.ra[[5]][all.sel.ra]  <- all.predY.95[,1]
sp.nam <- colnames(all.predY.mean)[1]
names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
writeRaster(pred.ra, file=paste0(pa.base.str,sp.nam,".tif"), overwrite=TRUE)
for(j in 2:ncol(all.predY.mean)){
  print(j)
  pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
  pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,j]
  pred.ra[[2]][all.sel.ra]  <- all.predY.median[,j]
  pred.ra[[3]][all.sel.ra]  <- all.predY.se[,j]
  pred.ra[[4]][all.sel.ra]  <- all.predY.5[,j]
  pred.ra[[5]][all.sel.ra]  <- all.predY.95[,j]
  print("calculations done, saving now")
  sp.nam <- colnames(all.predY.mean)[j]
  names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
  writeRaster(pred.ra, file=paste0(pa.base.str,sp.nam,".tif"), overwrite=TRUE)
}
# library(foreach)
# library(doParallel)
# ## Set up parallel backend
# numCores <- detectCores() - 1  # use one less than the number of available cores
# cl <- makeCluster(numCores)
# registerDoParallel(cl)
# ## Parallel processing
# foreach(j = 2:ncol(all.predY.mean)) %dopar% {
#   print(j)
#   pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
#   pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,j]
#   pred.ra[[2]][all.sel.ra]  <- all.predY.median[,j]
#   pred.ra[[3]][all.sel.ra]  <- all.predY.se[,j]
#   pred.ra[[4]][all.sel.ra]  <- all.predY.5[,j]
#   pred.ra[[5]][all.sel.ra]  <- all.predY.95[,j]
#   print("calculations done, saving now")
#   sp.nam <- colnames(all.predY.mean)[j]
#   names(pred.ra) <- c(paste0(sp.nam, "_mean"), paste0(sp.nam, "_median"), paste0(sp.nam, "_se"), paste0(sp.nam, "_5"), paste0(sp.nam, "_95"))
#   writeRaster(pred.ra, file=paste0(pa.base.str, sp.nam, ".tif"), overwrite=TRUE)
# }
# ## Stop the cluster
# stopCluster(cl)



###### ABUNDANCE DATA:
#### load all prediction files into single files, starting with the first file and then loop
#### THIS CAN BE DONE MUCH FASTER BY SPLITTING THE JOBS RATHER THAN ONE LONG RUN (3h)
load(paste0(pred.dir,pred.files.ab[1]))
## which cells do we need to fill with data
row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
all.sel.ra <- sel.not.na[sel[row.numbers]]
## fill objects for each value
all.predY.mean <- predY.mean
all.predY.median <- predY.median
all.predY.se <- predY.se
all.predY.5 <- predY.5
all.predY.95 <- predY.95
## loop through all other species
for(i in 2:length(pred.files.ab)){
  print(i)
  load(paste0(pred.dir,pred.files.ab[i]))
  row.numbers <- as.numeric(rownames(xy.grid[sel.loop,]))
  all.sel.ra <- c(all.sel.ra, sel.not.na[sel[row.numbers]])
  all.predY.mean <- rbind(all.predY.mean, predY.mean)
  all.predY.median <- rbind(all.predY.median, predY.median)
  all.predY.se <- rbind(all.predY.se, predY.se)
  all.predY.5 <- rbind(all.predY.5, predY.5)
  all.predY.95 <- rbind(all.predY.95, predY.95)
}

## now fill raster cells with values from the prediction file, save mean and sd output for each species
colnames(all.predY.mean) <- gsub("[/()]", "", colnames(all.predY.mean))
pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,1]
pred.ra[[2]][all.sel.ra]  <- all.predY.median[,1]
pred.ra[[3]][all.sel.ra]  <- all.predY.se[,1]
pred.ra[[4]][all.sel.ra]  <- all.predY.5[,1]
pred.ra[[5]][all.sel.ra]  <- all.predY.95[,1]
sp.nam <- colnames(all.predY.mean)[1]
names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
writeRaster(pred.ra, file=paste0(ab.base.str,sp.nam,".tif"))
for(j in 2:ncol(all.predY.mean)){
  print(j)
  pred.ra <- c(empty.ra, empty.ra, empty.ra, empty.ra, empty.ra)
  pred.ra[[1]][all.sel.ra]  <- all.predY.mean[,j]
  pred.ra[[2]][all.sel.ra]  <- all.predY.median[,j]
  pred.ra[[3]][all.sel.ra]  <- all.predY.se[,j]
  pred.ra[[4]][all.sel.ra]  <- all.predY.5[,j]
  pred.ra[[5]][all.sel.ra]  <- all.predY.95[,j]
  sp.nam <- colnames(all.predY.mean)[j]
  print("calculations done, saving now")
  names(pred.ra) <- c(paste0(sp.nam,"_mean"), paste0(sp.nam,"_median"), paste0(sp.nam,"_se"), paste0(sp.nam,"_5"), paste0(sp.nam,"_95"))
  writeRaster(pred.ra, file=paste0(ab.base.str,sp.nam,".tif"))
}



#############################################
## calculate richness and total abundance based on individual species' predictions

ra.list <- list.files(dir4)
ra.list.pa <- ra.list[grep("1_pa_",ra.list)]
ra.list.pa <- ra.list.pa[grep("Biota",ra.list.pa)]
ra.list.ab <- ra.list[grep("1_abundcondpres_",ra.list)][-79]

## Presence-Absence and Richness
hmsc.maps.pa <- rast(paste0(dir4,ra.list.pa))
hmsc.maps.pa.mean   <- subset(hmsc.maps.pa, seq(1,nlyr(hmsc.maps.pa),by=5))
hmsc.maps.pa.median <- subset(hmsc.maps.pa, seq(2,nlyr(hmsc.maps.pa),by=5))
hmsc.maps.pa.se     <- subset(hmsc.maps.pa, seq(3,nlyr(hmsc.maps.pa),by=5))
hmsc.maps.pa.5      <- subset(hmsc.maps.pa, seq(4,nlyr(hmsc.maps.pa),by=5))
hmsc.maps.pa.95     <- subset(hmsc.maps.pa, seq(5,nlyr(hmsc.maps.pa),by=5))

hmsc.maps.richness.median <- sum(hmsc.maps.pa.median, na.rm=TRUE)
writeRaster(hmsc.maps.richness.median, file=paste0(dir4,substr(ra.list.pa[[1]],1,58),"richness_median.tif"))
hmsc.maps.richness.mean <- sum(hmsc.maps.pa.mean, na.rm=TRUE)
writeRaster(hmsc.maps.richness.mean, file=paste0(dir4,substr(ra.list.pa[[1]],1,58),"richness_mean.tif"))
hmsc.maps.richness.se <- sum(hmsc.maps.pa.se, na.rm=TRUE)
writeRaster(hmsc.maps.richness.se, file=paste0(dir4,substr(ra.list.pa[[1]],1,58),"richness_se.tif"))
hmsc.maps.richness.5 <- sum(hmsc.maps.pa.5, na.rm=TRUE)
writeRaster(hmsc.maps.richness.5, file=paste0(dir4,substr(ra.list.pa[[1]],1,58),"richness_lci.tif"))
hmsc.maps.richness.95 <- sum(hmsc.maps.pa.95, na.rm=TRUE)
writeRaster(hmsc.maps.richness.95, file=paste0(dir4,substr(ra.list.pa[[1]],1,58),"richness_uci.tif"))

## Abundance
hmsc.maps.ab <- rast(paste0(dir4,ra.list.ab))
hmsc.maps.ab.mean   <- subset(hmsc.maps.ab, seq(1,nlyr(hmsc.maps.ab),by=5))
hmsc.maps.ab.median <- subset(hmsc.maps.ab, seq(2,nlyr(hmsc.maps.ab),by=5))
hmsc.maps.ab.se     <- subset(hmsc.maps.ab, seq(3,nlyr(hmsc.maps.ab),by=5))
hmsc.maps.ab.5      <- subset(hmsc.maps.ab, seq(4,nlyr(hmsc.maps.ab),by=5))
hmsc.maps.ab.95     <- subset(hmsc.maps.ab, seq(5,nlyr(hmsc.maps.ab),by=5))
# hmsc.maps.mean.corrected_1 <- hmsc.maps.ab.mean[[1:25]]
# hmsc.maps.mean.corrected_2 <- hmsc.maps.ab.mean[[26:50]]
# hmsc.maps.mean.corrected_3 <- hmsc.maps.ab.mean[[51:78]]
# for(i in 1:25){
#   print(i)
#   sel.ra.cells1 <- which(hmsc.maps.mean.corrected_1[[i]][]>540)
#   sel.ra.cells2 <- which(hmsc.maps.mean.corrected_2[[i]][]>540)
#   hmsc.maps.mean.corrected_1[[i]][sel.ra.cells1] <- 540
#   hmsc.maps.mean.corrected_2[[i]][sel.ra.cells2] <- 540
# }
# for(i in 1:28){
#   print(i)
#   sel.ra.cells3 <- which(hmsc.maps.mean.corrected_3[[i]][]>540)
#   hmsc.maps.mean.corrected_3[[i]][sel.ra.cells3] <- 540
# }
# hmsc.maps.mean.corrected <- c(hmsc.maps.mean.corrected_1, hmsc.maps.mean.corrected_2, hmsc.maps.mean.corrected_3)
# rm(hmsc.maps.mean.corrected_1, hmsc.maps.mean.corrected_2, hmsc.maps.mean.corrected_3)
# writeRaster(hmsc.maps.mean.corrected, file=paste0(dir4,substr(ra.list.ab[[1]],1,69),"zMeanWithOverpredictionsReducedTo100%Cover.tif"))
hmsc.maps.mean.corrected <- rast(paste0(dir4,substr(ra.list.ab[[1]],1,69),"zMeanWithOverpredictionsReducedTo100%Cover.tif"))
hmsc.maps.totalabundcondpres.mean <- sum(hmsc.maps.mean.corrected, na.rm=TRUE)/5.4
writeRaster(hmsc.maps.totalabundcondpres.mean, file=paste0(dir4,substr(ra.list.ab[[1]],1,69),"totalabundance_corrected.tif"))


##################
### Bootstrapping the median probabilities to get standard errors for richness
species_richness <- function(data) {
  richness <- rowSums(data, na.rm=TRUE)
  return(richness)
}
# Number of bootstrap samples (running 10 times 100 because it crashes otherwise)
n_bootstrap <- 100
for(k in 1:10){
  message(k)
# empty raster to save the output
  richness_bootstrapped <- rast(hmsc.maps.pa.median, nlyrs=100)
  for(i in 1:n_bootstrap){
    print(i)
    # bootstrap species
    sample_indices <- sample(1:nlyr(hmsc.maps.pa.median), replace = TRUE)
    # run bootstrap calculation of richness
    richness_bootstrapped[[i]] <- app(hmsc.maps.pa.median[[sample_indices]], fun=species_richness)
  }  
  writeRaster(richness_bootstrapped, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_bootstrapped_",k,".tif"))
}
# compile outputs and create raster
boots.files <- list.files(dir4, pattern="bootstrapped")
rich.ra.bootstrapped <- rast(paste0(dir4,boots.files))
rich.sd.bootstrapped <- app(rich.ra.bootstrapped, fun="sd")
rich.median.bootstrapped <- app(rich.ra.bootstrapped, fun="median")

writeRaster(rich.sd.bootstrapped, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_bootstrapped_sd.tif"))
writeRaster(rich.median.bootstrapped, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_bootstrapped_median.tif"))





###############################
###### calculate beta-diversity
library(adespatial)

# #########
# #### testing on a large region
# test.ra <- crop(hmsc.maps.pa.median, ext(-500000,500000,-2100000,-1200000))
# plot(test.ra[[1]])
# 
# dat <- values(test.ra)
# na.sel <- which(is.na(rowSums(dat)))
# dat.clean <- dat[-na.sel,]
# 
# ## 20 runs with every 20th cell
# #dat.beta.lcbd <- rep(NA, nrow(dat.clean))
# step_size <- 20
# sample.size <- floor(nrow(dat.clean)/step_size)
# sample.size1 <- nrow(dat.clean)-(step_size-1)*sample.size
# ## get samples right:
# selected_indices_list <- list()
# for(i in 1:step_size){
#   print(i)
#   if(i == 1){
#     selected_indices_list[[i]] <- sample(1:nrow(dat.clean), sample.size1, replace = FALSE)
#     selected_indices_full <<- selected_indices_list[[i]]
#   }else{
#     # Generate a random subset of indices for the current loop, excluding already selected indices
#     available_indices <- setdiff(1:nrow(dat.clean), selected_indices_full)
#     selected_indices_list[[i]] <- sample(available_indices, sample.size, replace = FALSE)
#     # Update the selected indices
#     selected_indices_full <<- c(selected_indices_full, selected_indices_list[[i]])
#   }
# }
# 
# for(i in 1:step_size){
#   print(i)
#   ##
#   dat.clean.loop <- dat.clean[selected_indices_list[[i]],]
#   dat.beta.loop <- adespatial::beta.div(dat.clean.loop)
#   if(i==1){
#     save(dat.beta.loop, selected_indices_list, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_TestSample",sprintf("%03d", i),".Rdata"))
#   }else{
#     save(dat.beta.loop, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_TestSample",sprintf("%03d", i),".Rdata"))
#   }
# }
# 
# dat.beta.lcbd <- rep(NA, length(selected_indices_full))
# for(i in 1:step_size){
#   print(i)
#   load(paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_TestSample",sprintf("%03d", i),".Rdata"))
#   dat.beta.lcbd[selected_indices_list[[i]]] <- dat.beta.loop$LCBD
# }
# 
# ross.ra.beta <- rast(test.ra[[1]])
# values(ross.ra.beta)[-na.sel] <- dat.beta.lcbd
# plot(ross.ra.beta)


#########
#### full run
dat <- values(hmsc.maps.pa.median)
na.sel <- which(is.na(rowSums(dat)))
dat.clean <- dat[-na.sel,]

## 50 runs with every 20th cell
#dat.beta.lcbd <- rep(NA, nrow(dat.clean))
step_size <- 100
sample.size <- floor(nrow(dat.clean)/step_size)
sample.size1 <- nrow(dat.clean)-(step_size-1)*sample.size
## get samples right:
selected_indices_list <- list()
for(i in 1:step_size){
  print(i)
  if(i == 1){
    selected_indices_list[[i]] <- sample(1:nrow(dat.clean), sample.size1, replace = FALSE)
    selected_indices_full <<- selected_indices_list[[i]]
  }else{
    # Generate a random subset of indices for the current loop, excluding already selected indices
    available_indices <- setdiff(1:nrow(dat.clean), selected_indices_full)
    selected_indices_list[[i]] <- sample(available_indices, sample.size, replace = FALSE)
    # Update the selected indices
    selected_indices_full <<- c(selected_indices_full, selected_indices_list[[i]])
  }
}

for(i in 1:step_size){
  print(i)
  ##
  dat.clean.loop <- dat.clean[selected_indices_list[[i]],]
  dat.beta.loop <- adespatial::beta.div(dat.clean.loop)
  if(i==1){
    save(dat.beta.loop, selected_indices_list, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_",sprintf("%03d", i),".Rdata"))
  }else{
    save(dat.beta.loop, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_",sprintf("%03d", i),".Rdata"))
  }
}
##
load(paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_",sprintf("%03d", 1),".Rdata"))
dat.beta.lcbd <- rep(NA, length(c(unlist(selected_indices_list))))
dat.beta.lcbd[selected_indices_list[[1]]] <- dat.beta.loop$LCBD
for(i in 2:step_size){
  print(i)
  load(paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_",sprintf("%03d", i),".Rdata"))
  dat.beta.lcbd[selected_indices_list[[i]]] <- dat.beta.loop$LCBD
}

ra.beta <- rast(hmsc.maps.pa.median[[1]])
values(ra.beta)[-na.sel] <- dat.beta.lcbd
plot(ra.beta)

writeRaster(ra.beta, file=paste0(dir4,"2km_model_cells_1_pa_chains_4_thin_10_samples_800_envonly_median_BetaDiversity_LCBD.tif"))

















# Initialize an array to store bootstrap estimates
bootstrap_estimates <- array(NA, dim = c(n_bootstrap, length(sel_na)))




for (chunk in 1:num_chunks) {
  message(chunk)
  start_index <- (chunk - 1) * chunk_size + 1
  end_index <- min(chunk * chunk_size, length(sel_na))
  
  # Perform bootstrap sampling for the current chunk
  for (i in 1:n_bootstrap) {
    if(i%in%seq(50,n_bootstrap, by=50)) print(i)
    print(paste("Chunk:", chunk, "Bootstrap sample:", i))
    sample_indices <- sample(1:nlyr(r), replace = TRUE)
    bootstrap_estimates[i, start_index:end_index] <- species_richness(r.dat.list[[chunk]], sample_indices)
  }
}





species_richness(data=r.crop[1:2,])


library(boot)
library(terra)
# Identify cells that are not NA
sel_na <- which(!is.na(r[[1]][]))
# Number of bootstrap samples
n_bootstrap <- 1000
# Initialize an array to store bootstrap estimates
bootstrap_estimates <- array(NA, dim = c(n_bootstrap, length(sel_na)))

# Process the raster in chunks
chunk_size <- 1000  # Adjust this size based on your system's memory capacity
num_chunks <- ceiling(length(sel_na) / chunk_size)

r.dat.list <- list()
for(chunk in 1:num_chunks) {
  message(chunk)
  start_index <- (chunk - 1) * chunk_size + 1
  end_index <- min(chunk * chunk_size, length(sel_na))
  sel_na_chunk <- sel_na[start_index:end_index]
  r.dat.list[[chunk]] <- r[sel_na_chunk,]
}

for (chunk in 1:num_chunks) {
  message(chunk)
  start_index <- (chunk - 1) * chunk_size + 1
  end_index <- min(chunk * chunk_size, length(sel_na))

  # Perform bootstrap sampling for the current chunk
  for (i in 1:n_bootstrap) {
    if(i%in%seq(50,n_bootstrap, by=50)) print(i)
    print(paste("Chunk:", chunk, "Bootstrap sample:", i))
    sample_indices <- sample(1:nlyr(r), replace = TRUE)
    bootstrap_estimates[i, start_index:end_index] <- species_richness(r.dat.list[[chunk]], sample_indices)
  }
}

# Calculate the standard error for each selected raster cell
standard_errors <- apply(bootstrap_estimates, 2, function(x) sd(x, na.rm = TRUE))

# Create a new raster for standard errors
se_raster <- raster(r)
values(se_raster)[sel_na] <- standard_errors



## Define the function to calculate the mean probability
calc_mean_prob <- function(data, indices) {
  sample_data <- data[indices]
  return(mean(sample_data, na.rm = TRUE))
}

# Load your raster data
# Assuming your raster is loaded as 'hmsc.maps.pa.env.median'
sel.na <- which(!is.na(values(hmsc.maps.pa.median[[1]])))
raster_values <- values(hmsc.maps.pa.median)[sel.na,]

# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data = raster_values, statistic = calc_mean_prob, R = 1000)

# Calculate the standard error
standard_error <- sd(bootstrap_results$t)
print(standard_error)

hmsc.maps.richness.env.bootstrapped.se






