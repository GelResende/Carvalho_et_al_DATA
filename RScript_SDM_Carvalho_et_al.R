############################################################################################################################
### Supporting Information to ###
setwd("E://")

# Install and load R packages needed to run the analysis:
needed_packages<-c("devtools", "raster", "magrittr", "dplyr", "stringr", "dismo", "stats", "plyr",
                   "foreach", "randomForest", "splitTools", "mapview", "ggspatial", "ggplot2",
                   "kernlab", "maxnet", "mgcv", "glmnet", "parallel", "doParallel", "purrr",
                   "rgdal", "sf", "sp", "data.table")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(needed_packages, require, character.only = TRUE)
removeTmpFiles(0.00001) # remove temporary files
rm(needed_packages, new.packages)

#####

# STEP 1 - PREPARE THE RASTER FILES ON THE BACKGROUND AREA
##########################################################################################################################
# STEP 1 - PREPARE THE RASTER FILES ON THE BACKGROUND AREA
rm(list=ls()); gc()

# Load environmental variables based on orthogonal PCA axes (PCA variables):
envT <- raster::stack(list.files("Rasters/", full.names=T))
names(envT)
names(envT)<-c("DryMonths", "LandTenure", "NearbtDND", "ResearchEduc", "TravelTime")

# Reclassify raster values and set the 65535 value to NA:
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE)

# Convert raster to points:
envT_df<-rasterToPoints(x=envT, spatial=FALSE)
envT_df[is.nan(envT_df)]<-NA # convert NaN values to NA
envT_df<-as.data.table(envT_df)
  
# Rename columns (predictors) to ease interpretability:
names(envT_df)<-c("x", "y", "DryMonths", "LandTenure", "NearbtDND", "ResearchEduc", "TravelTime")
summary(envT_df)

# Remove pixels outside legal Amazon:
envT_df<-envT_df[!is.na(envT_df$TravelTime),]

# Set 'LandTenure' to factor:
#envT_df$LandTenure<-as.factor(envT_df$LandTenure)

# Check the multicollinearity of the background area:
usdm::vif(envT_df[,3:ncol(envT_df)])

# Create the shapefile representing the grid cells of your study area:
envT_df$CompleteObs<-rowSums(envT_df[,3:7], na.rm=F)
length(which(is.na(envT_df$CompleteObs)))/nrow(envT_df)
grid_cells<-envT_df[,c(1:2,8)]
coordinates(grid_cells)<- ~ x + y # convert to spatial points
gridded(grid_cells)<-TRUE  # griddify your set of points

# Convert the SpatialPixelsDataframe into a raster object and assign coordinate reference system:
sampleable_cells <- raster(grid_cells) 
raster::crs(sampleable_cells) <- crs(envT)
sampleable_cells<-raster::reclassify(sampleable_cells, cbind(0, +Inf, 1), right=TRUE)
save(sampleable_cells, file="RData/SampleableCells.RData")

# Mask all predictors to consider only the sampleable cells:
cl <- parallel::makePSOCKcluster(parallel::detectCores()-1) # number of cores in computer
doParallel::registerDoParallel(cl)
getDoParWorkers()

# Rebuild projection raster to match the exact extent of the sampleable cells raster:
# Extract each raster layer in parallel: 
load("RData/SampleableCells.RData")
foreach(b = 1:nlayers(envT), .packages = c("raster", "sp", "rgdal")) %dopar% {
    
  # Load one layer for a future projection:
  SingleLayer<-raster::subset(envT, subset=b)
  
  # Crop raster layer to the extent of the 'sampleable_cells':
  SingleLayer<-raster::crop(SingleLayer, sampleable_cells)
  
  # Create an empty raster to receive the valuesof SingleLayer:
  EmptyRaster<-raster(ext=extent(sampleable_cells),
                      crs=crs(sampleable_cells),
                      nrows=dim(sampleable_cells)[1],
                      ncols=dim(sampleable_cells)[2])
  
  # Rewrite values over the extent of original predictors:
  RevaluedLayer <- raster::raster(vals=values(SingleLayer),
                                  ext=extent(EmptyRaster),
                                  crs=crs(EmptyRaster),
                                  nrows=dim(EmptyRaster)[1],
                                  ncols=dim(EmptyRaster)[2])
  
  # Write the raster:
  return(writeRaster(RevaluedLayer, filename=paste0("Predictors/", names(envT)[b], ".tif"), format="GTiff", overwrite=TRUE))
    
  } # end of foreach

parallel::stopCluster(cl)

#####

# STEP 2 - GET THE TRAINING AND TESTING DATA FOLDS
##########################################################################################################################
# STEP 2 - GET THE TRAINING AND TESTING DATA FOLDS
rm(list=ls()); gc()

# Read txt file with occurrence data:
all_occ <- data.table::fread("Datasets/wetland_points.csv", h = T, sep = ',', stringsAsFactors = F)
names(all_occ)<-c("group", "x", "y")

# Load raster data on sampleable cells:
load("RData/SampleableCells.RData")

# Identify which occurrence records will not result in NA predictors:
NotNA_occ<-data.frame(raster::extract(sampleable_cells, all_occ[,2:3], na.rm=T))
NotNA_occ<-which(!is.na(NotNA_occ)) # only two occurrence records are NA
occ<-all_occ[NotNA_occ,]
fwrite(occ, "Datasets/CleannedOcc.txt", sep="\t")

# Load environmental variables based on orthogonal PCA axes (PCA variables):
envT <- raster::stack(list.files("Predictors/", full.names=T))
occ <- data.table::fread("Datasets/CleannedOcc.txt", h = T, stringsAsFactors = F)

# Create a new column indicating the presence value: 
occ$Presence<-1

# Create a list of occurrence data for each group:
occ_xy <- split(occ[,-1], f=occ$group)
spN <- names(occ_xy) # get group names

# Function to create random points (to use next):
Fast_Random_Points <- function(r, n, p) {
  v <- raster::getValues(r) # get values from the raster
  v <- which(!is.na(v)) # remove NA
  v <- v[!v%in%raster::cellFromXY(r,p)] # select pixels not occupied by an occurrence
  v <- sample(v, n) # randomly sample n pixels from v
  v <- raster::xyFromCell(r, v) # return as coordinates
  return(v)
}

# Get the set of 10K random points for each group:
absences<-list()
for(i in 1:length(spN)){
  
  # Create the same number of presences as pseudoabsences:
  new_absences <- Fast_Random_Points(r=sampleable_cells, # mask for pseudoabsences
                                     n=nrow(occ_xy[[i]]), # number of pseudoabsences (same as presences)
                                     p=occ_xy[[i]] # list of occurrences (presences)
                                     )
  
  # Create a new column informing that the presence is false:
  new_absences<-as.data.frame(new_absences)
  new_absences$Presence<-0
  
  # Store in a list:
  absences[[i]]<-new_absences
  rm(new_absences)
  
}
names(absences)<-spN

# For each group, create the training and testing data partitions for use in cross-validation:
TrainingData<-list()
TestingData<-list()
for(i in 1:length(spN)){
  
  # Set seed to allow reproducibility:
  set.seed(123)
  
  # Create empty lists to store the training and testing data of the group 'i':
  TR_Data<-list()
  TS_Data<-list()
  
  # Get a vector of row numbers used in each partition for the group 'i':
  k_index_1 <- caret::createFolds(y = 1:nrow(occ_xy[[i]]), k = 5, list=TRUE, returnTrain=TRUE) # presences
  k_index_0 <- caret::createFolds(y = 1:nrow(absences[[i]]), k = 5, list=TRUE, returnTrain=TRUE) # absences
  
  for(k in 1:length(k_index_1)){
    
    # Store the training fold 'k ' for the group 'i':
    TR_presences<-occ_xy[[i]][k_index_1[[k]],]
    TR_absences<-absences[[i]][k_index_0[[k]],]
    
    # Same as above, but for testing data:
    TS_presences<-occ_xy[[i]][-k_index_1[[k]],]
    TS_absences<-absences[[i]][-k_index_0[[k]],]
    
    # Store each k partition:
    TR_Data[[k]]<-data.frame(rbind(TR_presences, TR_absences), Kfold=k)
    TS_Data[[k]]<-data.frame(rbind(TS_presences, TS_absences), Kfold=k)
    
    # Remove unnecessary objects before next iteration:
    rm(TR_presences, TR_absences, TS_presences, TS_absences)
    
  }
  
  # Name lists TR_Data and TS_Data:
  names(TR_Data)<-c("Kfold1", "Kfold2", "Kfold3","Kfold4", "Kfold5")
  names(TS_Data)<-c("Kfold1", "Kfold2", "Kfold3","Kfold4", "Kfold5")
  
  # Store all data partitions for the group 'i':
  TrainingData[[i]]<-TR_Data
  TestingData[[i]]<-TS_Data
  
  # Remove unnecessary objects before next iteration:
  rm(TR_Data, TS_Data, k_index_1, k_index_0)
  
}
names(TrainingData)<-spN
names(TestingData)<-spN

# Export to disk:
save(TestingData, TrainingData, file="RData/TR_and_TS_Data.RData")

#####

# STEP 3 - BUILD RANDOM FOREST MODELS
##########################################################################################################################
# STEP 3 - BUILD RANDOM FOREST MODELS
rm(list=ls()); gc()

# Load environmental variables based on orthogonal PCA axes (PCA variables):
envT <- raster::stack(list.files("Predictors/", full.names=T))
envT[[2]]<-as.factor(envT[[2]]) # convert 'LandTenure' to factor
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) # reclassify raster values and set the 65535 value to NA

# Load training and testing data:
load("RData/TR_and_TS_Data.RData")

# Create an empty list to receive model outputs:
RF_TrainingOutput<-list()
RF_TestingOutput<-list()

# For each group 'i', :
a<-Sys.time()
for(i in 1:length(TrainingData)){
  
  # Get the environmental data for each k-training and k-testing fold of the group 'i':
  TR_envData<-lapply(TrainingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  TS_envData<-lapply(TestingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  
  # Create an empty list to receive model outputs:
  RF_training_k<-list()
  RF_testing_k<-list()
  
  for(k in 1:length(TS_envData)){
  
    # Bind the spatial and environmental datasets:
    MyTRData<-cbind(TrainingData[[i]][[k]], TR_envData[[k]])
    MyTSData<-cbind(TestingData[[i]][[k]], TS_envData[[k]])
    MyTRData<-MyTRData[complete.cases(MyTRData),] # remove NA's (if any)
    MyTSData<-MyTSData[complete.cases(MyTSData),] # remove NA's (if any)
    
    # Convert 'LandTenure' to factor and certify all levels appear in both training and testing datasets:
    MyTRData$LandTenure<-as.factor(MyTRData$LandTenure)
    MyTSData$LandTenure<-as.factor(MyTSData$LandTenure)
    MyTRData$LandTenure<-factor(MyTRData$LandTenure, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    MyTSData$LandTenure<-factor(MyTSData$LandTenure, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    
    # Build and store model outputs for the partition 'k' and group 'i':
    RF_training_k[[k]]<-randomForest::tuneRF(x = MyTRData[,c(5:9)], # predictors
                                             y = MyTRData[,3], # response
                                             trace = F,
                                             stepFactor = 1,
                                             ntreeTry = 1000,
                                             doBest = T, # repeat model using tuning parameters
                                             plot = F,
                                             importance = T, # store variable importance
                                             )
    
    # Add a new list element holding model outputs for each data fold:
    RF_training_k[[k]][[ (length(RF_training_k[[k]])+1) ]]<-data.frame(tree=1:length(RF_training_k[[k]]$mse),
                                                                       group=names(TrainingData)[i],
                                                                       kfold=k,
                                                                       mse=RF_training_k[[k]]$mse,
                                                                       rsq=RF_training_k[[k]]$rsq)
    names(RF_training_k[[k]])[length(RF_training_k[[k]])]<-"RF_output_k" # name the last element of the list
    
    # Store variable importance computed for each k-fold:
    RF_training_k[[k]][[ (length(RF_training_k[[k]])+1) ]]<-as.data.frame(RF_training_k[[k]]$importance)
    names(RF_training_k[[k]])[length(RF_training_k[[k]])]<-"importance_k" # name the last element of the list
    
    # Get predictions for the testing data:
    RF_testing_k[[k]]<-data.frame(group=names(TrainingData)[i],
                                  MyTSData, 
                                  Predicted=predict(object=RF_training_k[[k]], newdata=MyTSData[,c(5:9)], type="response"))
    
  } # end of k for loop
  
  # Store outputs for the group 'i':
  RF_TrainingOutput[[i]]<-randomForest::combine(RF_training_k[[1]], # get the original lists
                                                RF_training_k[[2]],
                                                RF_training_k[[3]], 
                                                RF_training_k[[4]],
                                                RF_training_k[[5]])
  
  # Add a new list element holding model outputs for all data folds:
  RF_TrainingOutput[[i]] [[ (length(RF_TrainingOutput[[i]])-1) ]] <-data.frame(rbind(
                   RF_training_k[[1]][[(length(RF_training_k[[1]])-1)]],
                   RF_training_k[[2]][[(length(RF_training_k[[2]])-1)]],
                   RF_training_k[[3]][[(length(RF_training_k[[3]])-1)]],
                   RF_training_k[[4]][[(length(RF_training_k[[4]])-1)]],
                   RF_training_k[[5]][[(length(RF_training_k[[5]])-1)]]))
  
  # Add a new list element holding variable importance values for all data folds:
  RF_TrainingOutput[[i]] [[ length(RF_TrainingOutput[[i]]) ]] <-data.frame(rbind(
    RF_training_k[[1]][[length(RF_training_k[[1]])]],
    RF_training_k[[2]][[length(RF_training_k[[2]])]],
    RF_training_k[[3]][[length(RF_training_k[[3]])]],
    RF_training_k[[4]][[length(RF_training_k[[4]])]],
    RF_training_k[[5]][[length(RF_training_k[[5]])]]))
  
  RF_TestingOutput[[i]]<-rbindlist(RF_testing_k)

  # Remove unnecessary objects before next iteration:
  rm(RF_training_k, RF_testing_k, TR_envData, TS_envData, MyTRData, MyTSData)
  
} # end of i for loop
b<-Sys.time()
b-a

# Export to disk:
save(RF_TrainingOutput, RF_TestingOutput, file="RData/TR_and_TS_ModelOutputs.RData")

#####

# STEP 4 - ASSESS VARIABLE IMPORTANCE AND METRICS OF MODEL PERFORMANCE AND EVALUATION
##########################################################################################################################
# STEP 4 - ASSESS VARIABLE IMPORTANCE AND METRICS OF MODEL PERFORMANCE AND EVALUATION
rm(list=ls()); gc()

# Load model output for training and testing data:
load("RData/TR_and_TS_ModelOutputs.RData")

# Combine all testing outputs in a single dataframe:
TestingOutput<-rbindlist(RF_TestingOutput)

# Extract outputs of model performance for each training fold:
ModelPerformance <- lapply(RF_TrainingOutput, purrr::pluck, "RF_output_k")

# Perform model evaluation for the complete set of predicted values using the testing data:
ModelEvaluation<-list()
TestingOutput$group<-as.factor(TestingOutput$group)
for(i in 1:nlevels(TestingOutput$group)){
  
  ModelEvaluation[[i]]<-ENMTML:::Eval_Jac_Sor_TMLA(
    p = TestingOutput[TestingOutput$group==levels(TestingOutput$group)[i] & TestingOutput$Presence==1,]$Predicted,
    a = TestingOutput[TestingOutput$group==levels(TestingOutput$group)[i] & TestingOutput$Presence==0,]$Predicted,
    thr = NULL)
  
  names(ModelEvaluation)[i]<-levels(TestingOutput$group)[i]

}

ModelEvaluation$ant$SorensenTHR

##report the threshold of models
ModelEvaluation$ant$SorensenTHR
ModelEvaluation$beetle$SorensenTHR
ModelEvaluation$tree$SorensenTHR
ModelEvaluation$tree$SorensenTHR

###report  sorensenof the models
max(ModelEvaluation$ant$Sorensen, na.rm=T)
max(ModelEvaluation$beetle$Sorensen, na.rm=T)
max(ModelEvaluation$bird$Sorensen, na.rm=T)
max(ModelEvaluation$tree$Sorensen, na.rm=T)

#TestingOutput

######eexporta matrix
as.matrix(ModelEvaluation)->ModelEvaluation
write.table (ModelEvaluation, 'ModelEvaluation.csv')
save(ModelEvaluation, file="ModelEvaluation.RData")


#ModelEvaluation<-rbindlist(ModelEvaluation)
#ModelEvaluation
#ModelEvaluation<-as.data.table(ModelEvaluation)

#######################################################
ModelEvaluation<-rbindlist(ModelEvaluation)
ModelEvaluation
ModelEvaluation<-ModelEvaluation[, .(.N,
                                       avg_Sorensen = mean(Sorensen, na.rm=T),
                                       min_Sorensen = min(Sorensen, na.rm=T),
                                       max_Sorensen = max(Sorensen, na.rm=T),
                                       sd_Sorensen = sd(Sorensen, na.rm=T),
                                       avg_Jaccard= mean(Jaccard, na.rm=T),
                                       min_Jaccard = min(Jaccard, na.rm=T),
                                       max_Jaccard = max(Jaccard, na.rm=T),
                                       sd_Jaccard = sd(Jaccard, na.rm=T)),
                                   by=.(presences,absences)] 
ModelEvaluation

##na verdade o sorensen de acurácia, é o máximo da função acima

# Compute the average metrics of model performance for each group across each iteration (tree) and folds (k):
ModelPerformance<-rbindlist(ModelPerformance)
ModelPerformance
ModelPerformance<-ModelPerformance[, .(.N,
                                        avg_mse = mean(mse, na.rm=T),
                                        min_mse = min(mse, na.rm=T),
                                        max_mse = max(mse, na.rm=T),
                                        sd_mse = sd(mse, na.rm=T),
                                        avg_rsq = mean(rsq, na.rm=T),
                                        min_rsq = min(rsq, na.rm=T),
                                        max_rsq = max(rsq, na.rm=T),
                                        sd_rsq = sd(rsq, na.rm=T)),
                                    by=.(group, tree)] 

# Filter the model performance output for the last tree grown:
ModelPerformance<-ModelPerformance[ModelPerformance$tree==500,]
ModelPerformance

# Extract outputs of model performance for each training fold:
VarImportance <- lapply(RF_TrainingOutput, purrr::pluck, "importance_k")
names(VarImportance)<-ModelPerformance$group # register the name of the group as list name

# Compute normalized metrics of variable importance per group:
VarImp_Normalized<-data.frame()
for(i in 1:length(VarImportance)){
  
  VarImportance[[i]]$Kfold<-c(rep(1,5), rep(2,5), rep(3,5), rep(4,5), rep(5,5))
  VarImportance[[i]]$Variables<-rep(row.names(RF_TrainingOutput[[1]]$importance), 5)
  
  for(k in 1:5){
    
    MyData<-VarImportance[[i]][VarImportance[[i]]$Kfold==k,]
    MyData$INP_Normalized<-(MyData$IncNodePurity)/(sum(MyData$IncNodePurity))
    MyData$IMSE_Normalized<-(MyData$X.IncMSE)/(sum(MyData$X.IncMSE))
    MyData$group<-ModelPerformance$group[i]
    VarImp_Normalized<-rbind(VarImp_Normalized, MyData)
  }
}
rm(VarImportance)

# Prepare data.frame for plotting:
VarImp_Normalized<-as.data.table(VarImp_Normalized)
VarImportance<-VarImp_Normalized[, .(avg_INP = mean(INP_Normalized, na.rm=T),
                                     min_INP = min(INP_Normalized, na.rm=T),
                                     max_INP = max(INP_Normalized, na.rm=T),
                                     sd_INP = sd(INP_Normalized),
                                     avg_IMSE = mean(IMSE_Normalized, na.rm=T),
                                     min_IMSE = min(IMSE_Normalized, na.rm=T),
                                     max_IMSE = max(IMSE_Normalized, na.rm=T),
                                     sd_IMSE = sd(IMSE_Normalized, na.rm=T)),
                                   by=.(group, Variables)]

VarImportance
# Relabel levels to improve visualization:
VarImportance$group1 <- factor(VarImportance$group, 
                            levels=c("tree", "bird", "beetle", "ant"),
                            labels=c("Trees", "Birds", "Beetles", "Ants"))
VarImportance$Variables1 <- factor(VarImportance$Variables, 
                                   levels=rev(c("LandTenure", "DryMonths", "NearbtDND", "ResearchEduc", "TravelTime")),
                                   labels=rev(c("Land\nTenure", "Dry season\nLength", "Degradation", "Education", "Accessibility")))

# Plot the normalized IncNodePurity for each variable and group:
MyPlot1<-ggplot(data=VarImportance, aes(x=group1, y=avg_INP, shape=group1)) +
  
  # Add geopoints:
  geom_point(shape=16, aes(colour=group1, shape=group1, fill=group1), size=2) + 
  
  # Add error bars representing min and max values of variable importance:
  geom_errorbar(aes(ymin=min_INP, ymax=max_INP, col=group1), width=0.5, cex=0.65) +
  
  # Set the colour ramp:
  scale_color_manual(values = c("#4F724B", "#4F724B", "#4F724B", "#4F724B"), labels = c("Ants", "Beetles", "Birds", "Trees")) +
  
  # Split geom_points according to levels of 'Variables':
  facet_wrap(~Variables1, strip.position="left", nrow=6, scales = "free_y") +
  
  # Axes titles:
  xlab("") + ylab('Mean Decrease Gini (Normalized)') +
  
  # Other aesthetics:
  coord_flip() +  # flip coordinates (puts labels on y axis)
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)),
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)),
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        strip.placement = "outside")

MyPlot1

# Same as above, but for the normalized %IncMSE (Proportional Increase in Mean Squared Error):
MyPlot2<-ggplot(data=VarImportance, aes(x=group1, y=avg_IMSE)) +
  geom_point(shape=16, aes(colour=group1, fill=group1), size=2) + 
  geom_errorbar(aes(ymin=min_IMSE, ymax=max_IMSE, col=group1), width=0.5, cex=0.65) +
  scale_color_manual(values = c("#C51B7D", "#B35806", "#542788", "#4D9221"),labels = c("Ants", "Beetles", "Birds", "Trees")) +
  facet_wrap(~Variables1, strip.position="left", nrow=6, scales = "free_y") +
  xlab("") + ylab('% Increase in MSE (Normalized)') +
  coord_flip() +  
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.line = element_line(colour="black"), 
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)),
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)),
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        strip.placement = "outside")

# Build the multipanel plot:
MultipanelPlot<-ggpubr::ggarrange(MyPlot1 + theme(legend.justification = c(0, 0),
                                                  legend.position=c(0.75, 0.05),
                                                  legend.background=element_rect(fill="white", size=0.3, linetype="solid", colour="black"),
                                                  legend.title=element_text(size=0, hjust=0),
                                                  legend.spacing.x = unit(0.25, 'cm'),
                                                  legend.text=element_text(size=10),
                                                  legend.key.width=unit(0.4,"cm"),
                                                  legend.key.height=unit(0.2,"cm"),
                                                  legend.key = element_blank(),
                                                  legend.margin = margin(0.12, 0.12, 0.12, 0.12, "cm")) +
                                             guides(shape="none", fill="none", color=guide_legend(ncol=1, byrow=F, reverse=F)),
                                  MyPlot2,
                                  labels=c(" A", " B"), 
                                  font.label=list(size=12, color = "black"), ncol=2, nrow=1)

# Export to disk:
ggsave("Figures/VariableImportance.png", plot=MultipanelPlot, width=10, height=5, units="in", bg="transparent")
ggsave("Figures/VariableImportance.pdf", plot=MultipanelPlot, width=10, height=5, units="in", bg="transparent")

#####

# STEP 5 - PREDICT THE GEOGRAPHICAL PATTERNS OF RESEARCH SUITABILITY
##########################################################################################################################
# STEP 5 - PREDICT THE GEOGRAPHICAL PATTERNS OF RESEARCH SUITABILITY
rm(list=ls()); gc()

# Load model output for training and testing data:
load("RData/TR_and_TS_ModelOutputs.RData")

# Load raster data:
envT <- raster::stack(list.files("Predictors/", full.names=T))
envT[[2]]<-as.factor(envT[[2]]) # convert 'LandTenure' to factor
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) # reclassify raster values and set the 65535 value to NA

# Get raster cell values as a dataframe:
envT_df<-rasterToPoints(x=envT, spatial=FALSE)
envT_df[is.nan(envT_df)]<-NA # convert NaN values to NA
envT_df<-as.data.table(envT_df) # convert to data.table
envT_df<-envT_df[!is.na(envT_df$TravelTime),] # remove pixels outside legal Amazon
envT_df$CompleteObs<-rowSums(envT_df[,3:7], na.rm=F)
envT_df$LandTenure<-as.factor(envT_df$LandTenure) # set 'LandTenure' as factor
envT_df<-envT_df[!is.na(envT_df$CompleteObs), c(1:7)]

# Predict 'research suitability' for each group in parallel:
cl<-makePSOCKcluster(detectCores()-1, type="SOCK")
registerDoParallel(cl)
getDoParWorkers()
MyPredictions<-foreach(i = 1:length(RF_TrainingOutput), 
                   .export = 'c', 
                   .packages = c("raster", "sp", "randomForest")) %dopar% {
                     predict(object = RF_TrainingOutput[[i]], newdata = envT_df,  type = "response")
                    }
parallel::stopCluster(cl)

# Convert the prediction tables in a SpatialPixelDataFrame:
MyPixelDF<-(data.frame(envT_df[,1:2], 
                       Ants=MyPredictions[[1]],
                       Beetles=MyPredictions[[2]],
                       Birds=MyPredictions[[3]],
                       Trees=MyPredictions[[4]]))
coordinates(MyPixelDF)<- ~ x + y # convert to spatial points
gridded(MyPixelDF)<-TRUE  # griddify your set of points

# Convert the SpatialPixelDataFrame into raster layers:
MyRasters<-list()
for(i in 1:ncol(MyPixelDF)){
  MyRasters[[i]] <- raster(MyPixelDF, layer=i) 
  }

# Convert the list of rasters to raster.stack:
MyRasters<-do.call(raster::stack, MyRasters)
raster::crs(MyRasters) <- raster::crs(envT)

# Export rasters to disk:
dir.create("PredictedMaps/", showWarnings = F)
for(i in 1:4){writeRaster(MyRasters[[i]], filename=paste0("PredictedMaps/", names(MyRasters)[i], ".tif"), format="GTiff", overwrite=TRUE)}

#####

# STEP 6: PLOT THE MAPS
########################################################################################################################
# STEP 6: PLOT THE MAPS
rm(list=ls()); gc()

### PLOTS:
# Let's add polygon boundaries for the background area:
BackgroundArea<-rgdal::readOGR(dsn="Shapefiles", layer='Amz_biome') # change directory as needed
BackgroundArea<-sp::spTransform(BackgroundArea, crs(raster::stack(list.files("Rasters/", full.names=T))))
BackgroundArea_sf<-sf::st_as_sf(BackgroundArea)

# Load raster files on each group. To visualize a raster in R using gplot2, we need to convert it to a dataframe: 
RasterData <- raster::stack(list.files("PredictedMaps/", full.names=T))
RasterData <- as.data.frame(RasterData, xy=TRUE)
group_names<-names(RasterData)

# Load spatial points:
OccRecords <- data.table::fread("Datasets/CleannedOcc.txt", h = T, stringsAsFactors = F)
OccRecords$group <- as.factor(OccRecords$group)
OccRecords<-SpatialPointsDataFrame(coords=OccRecords[,c("x", "y")], data=OccRecords[,.(group, x, y)])
wgs84<-"+proj=longlat +datum=WGS84 +no_defs"
crs(OccRecords)<-wgs84 # set the CRS
OccRecords_sf<-sf::st_as_sf(OccRecords) # convert to sf object

# Create an empty list to receive the plots:
MyColours<-list()
MyColours[[1]]<-scale_fill_gradientn(colours=c("#F7F7F7", "#FDE0EF", "#F1B6DA", "#DE77AE", "#C51B7D", "#8E0152"),
                                     values=c(0, 0.05, 0.1, 0.3, 0.5, 1), na.value = NA) # shades of pink
MyColours[[2]]<-scale_fill_gradientn(colours=c("#F7F7F7", "#FEE0B6", "#FDB863", "#E08214", "#B35806", "#7F3B08"), 
                                     values=c(0, 0.05, 0.1, 0.3, 0.5, 1), na.value = NA) # shades of orange
MyColours[[3]]<-scale_fill_gradientn(colours=c("#F7F7F7", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B"), 
                                     values=c(0, 0.05, 0.1, 0.3, 0.5, 1), na.value = NA) # shades of purple
MyColours[[4]]<-scale_fill_gradientn(colours=c("#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#276419"),
                                     values=c(0, 0.05, 0.1, 0.3, 0.5, 1), na.value = NA) # shades of green
MyColours[[5]]<-scale_fill_gradientn(colours=c("#F7F7F7", "#E6F5D0", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
                                     values=c(0, 0.05, 0.1, 0.3, 0.5, 1), na.value = NA) # shades of red

# Plot maps in a loop:
MyMaps<-list()
for(i in 1:4){
  
  # Get the group 'i' represented by the column 'group':
  RasterData$group<-RasterData[,(i+2)]
    
  # Map the research suitability:
  MyMaps[[i]] <- ggplot2::ggplot() +
    
    # Add the raster layer in the background:
    geom_raster(data = RasterData, aes(x=x, y=y, fill=group), na.rm=T) + 
    
    # Add polygon boundaries for the study area:
    geom_sf(data=BackgroundArea_sf, fill=NA, colour="black", size=0.3) +
    
    # Add occurrence records:
    geom_sf(data=OccRecords_sf[OccRecords_sf$group==levels(OccRecords_sf$group)[i],], colour="black", fill=NA, size=1, shape=3) +
    
    # Set the colour ramp:
    MyColours[[i]] +
    
    # Specify other aesthetics:
    theme(axis.line=element_blank(), # no axis line
          axis.text=element_blank(), # no axis text
          axis.ticks=element_blank(), # no tick marks
          axis.title=element_blank(), # no axis titles
          panel.grid.minor=element_blank(), # no minor grids (subdivisions of plotting space)
          panel.grid.major=element_blank(), # no major grids
          panel.background=element_blank(),
          plot.background=element_blank(), 
          plot.margin=unit(c(0, 0, 0, 0), "cm"),  # top, right, bottom, left
          panel.spacing=unit(c(0, 0, 0, 0), "cm"),  # top, right, bottom, left
          panel.border=element_rect(fill=NA, colour="black"),
          legend.position=c(0.9, 0.25), # x and y axis
          legend.direction="vertical", # legend orientation
          legend.title=element_text(size=10, hjust=0, face="bold"),
          legend.text=element_text(size=rel(0.8), hjust=0),
          legend.background=element_blank()
    ) +
    
    # Restrict the map to study area extent and define projection:
    coord_sf(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
             xlim = (raster::extent(BackgroundArea_sf)[1:2]+c(-1,1)),
             ylim = (raster::extent(BackgroundArea_sf)[3:4]+c(-1,1)), expand = FALSE) +
    
    # Adjust colour bar for the legend:
    guides(fill=guide_colourbar(title="Research\nPresence", label=T, nbin=15, barwidth=1, 
                                barheight=8, draw.ulim=T, draw.llim=T, frame.colour="black", 
                                frame.linewidth=1, ticks.linewidth=1, ticks=T, ticks.colour="black"))
}

# Plot the average result across groups:
GroupedRaster <- raster::stack(list.files("PredictedMaps/", full.names=T))
GroupedRaster<-raster::calc(x=GroupedRaster, fun=sum, na.rm=T)/4
GroupedRaster<-as.data.frame(GroupedRaster, xy=TRUE)
names(GroupedRaster)<-c("x", "y", "group")
MyMaps[[5]] <- ggplot2::ggplot() +
  
  # Add the raster layer in the background:
  geom_raster(data = GroupedRaster, aes(x=x, y=y, fill=group), na.rm=T) + # raster file
  
  # Add polygon boundaries for the study area:
  geom_sf(data=BackgroundArea_sf, fill=NA, colour="black", size=0.3) +
  
  # Add occurrence records:
  geom_sf(data=OccRecords_sf, colour="black", fill=NA, size=1, shape=3) +
  
  # Set the colour ramp:
  MyColours[[5]] +
  
  # Specify other aesthetics:
  theme(axis.line=element_blank(), # no axis line
        axis.text=element_blank(), # no axis text
        axis.ticks=element_blank(), # no tick marks
        axis.title=element_blank(), # no axis titles
        panel.grid.minor=element_blank(), # no minor grids (subdivisions of plotting space)
        panel.grid.major=element_blank(), # no major grids
        panel.background=element_blank(),
        plot.background=element_blank(), 
        plot.margin=unit(c(0, 0, 0, 0), "cm"),  # top, right, bottom, left
        panel.spacing=unit(c(0, 0, 0, 0), "cm"),  # top, right, bottom, left
        panel.border=element_rect(fill=NA, colour="black"),
        legend.position=c(0.9, 0.25), # x and y axis
        legend.direction="vertical", # legend orientation
        legend.title=element_text(size=10, hjust=0, face="bold"),
        legend.text=element_text(size=rel(0.8), hjust=0),
        legend.background=element_blank()
  ) +
  
  # Restrict the map to study area extent and define projection:
  coord_sf(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
           xlim = (raster::extent(BackgroundArea_sf)[1:2]+c(-1,1)),
           ylim = (raster::extent(BackgroundArea_sf)[3:4]+c(-1,1)), expand = FALSE) +
  
  # Adjust colour bar for the legend:
  guides(fill=guide_colourbar(title="Research\nPresence", label=T, nbin=15, barwidth=1, 
                              barheight=8, draw.ulim=T, draw.llim=T, frame.colour="black", 
                              frame.linewidth=1, ticks.linewidth=1, ticks=T, ticks.colour="black"))

# Build a multipanel plot:
MultipanelPlot<-ggpubr::ggarrange(MyMaps[[1]], MyMaps[[2]], MyMaps[[3]], MyMaps[[4]],
                                  labels=c("  A", "  B", "  C", "  D"), 
                                  font.label=list(size=12, color = "black"), ncol=2, nrow=2) +
  annotate("text", label="Ants", x = 0.025, y = 0.95, size = 4, colour = "black", hjust=0, angle=0, fontface=2) +
  annotate("text", label="Beetles", x = 0.525, y = 0.95, size = 4, colour = "black", hjust=0, angle=0, fontface=2) +
  annotate("text", label="Birds", x = 0.025, y = 0.45, size = 4, colour = "black", hjust=0, angle=0, fontface=2) +
  annotate("text", label="Trees", x = 0.525, y = 0.45, size = 4, colour = "black", hjust=0, angle=0, fontface=2)

ggsave("Figures/MultipanelMaps.pdf", plot=MultipanelPlot, width=12, height=8.5, units="in", bg="transparent")
ggsave("Figures/MultipanelMaps.png", plot=MultipanelPlot, width=12, height=8.5, units="in", bg="transparent")
ggsave("Figures/AverageMap.png", plot=MyMaps[[5]], width=6, height=5, units="in", bg="transparent")

#####

# STEP 7 - BUILD PARTIAL PLOTS
##########################################################################################################################
# STEP 7 - BUILD PARTIAL PLOTS
rm(list=ls()); gc()

# Load model training and testing data:
load("RData/TR_and_TS_Data.RData")

# Load model output for training and testing data:
load("RData/TR_and_TS_ModelOutputs.RData")

# Load environmental variables based on orthogonal PCA axes (PCA variables):
envT <- raster::stack(list.files("Predictors/", full.names=T))
envT[[2]]<-as.factor(envT[[2]]) # convert 'LandTenure' to factor
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) # reclassify raster values and set the 65535 value to NA

# For each group 'i': 1 = ants, 2 = beetles, 3 = birds, 4 = trees
Partial_Plots<-list()
for(i in 1:length(TestingData)){ 
  
  # Extract the environmental data for all presences/pseudoabsences records:  
  All_envData<-lapply(TestingData[[i]], function(x)
  raster::extract(envT, x[,1:2], na.rm=T))
  All_envData<-lapply(All_envData, as.data.frame)
  
  # Combined all data folds used for trainning:
  MyData<-cbind(rbindlist(TestingData[[i]]), rbindlist(All_envData))
  MyData<-MyData[complete.cases(MyData),] # remove NA's (if any)
  MyData$LandTenure<-as.factor(MyData$LandTenure)
  MyData$LandTenure<-factor(MyData$LandTenure, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  
  # Separate the randomForest build for te group 'i':
  SelectedRF_Model<-RF_TrainingOutput[[i]]
  
  # Create an empty list to receive model outputs:
  myPlots<-list()
  PredictorNames<-names(MyData)[5:9]
  # OBS: PredictorNames[2] refers to the categorical variable 'LandTenure'
    
  # Build Partial plots for all variables:
    for(v in 1:5){ # v number of variables
    
      # Plots for continuous predictors:
      if(v!=2){
        
        myPlots[[v]]<-partialPlot(x = SelectedRF_Model, # object of class randomForest
                                  pred.data = MyData, # training data for the respective randomForest
                                  x.var = PredictorNames[v], # variable name to be examined
                                  which.class = NULL, # for classification data, define which class
                                  plot = TRUE, 
                                  add = FALSE,
                                  rug = TRUE, 
                                  ylab="Research Presence",
                                  xlab=PredictorNames[v]
                                  )
        
        myPlots[[v]]$Group<-names(TestingData)[i]
        myPlots[[v]]$Predictor<-PredictorNames[v]
        
      }
      
      # Plots for categorical predictors:
      if(v==2){
        
        myPlots[[v]]<-partialPlot(x = SelectedRF_Model, # object of class randomForest
                                  pred.data = MyData, # training data for the respective randomForest
                                  x.var = PredictorNames[v], # variable name to be examined
                                  which.class = NULL, # for classification data, define which class
                                  plot = TRUE, 
                                  add = FALSE,
                                  rug = TRUE, 
                                  ylab="Research Presence",
                                  xlab=PredictorNames[v]
        )
        
        myPlots[[v]]$Group<-names(TestingData)[i]
        myPlots[[v]]$Predictor<-PredictorNames[v]
      }
      
    }
  
  Partial_Plots[[i]]<-myPlots
  
}
Partial_Plots[[1]]
Partial_Plots[[2]]
Partial_Plots[[3]]
Partial_Plots[[4]]

# Partial_Plots is a set of 4 lists, one for each group:
# Partial_Plots[[1]] contains five lists, one for each predictor variable
# Each list, has the x-value and y-values of the partialplot for the informed variable

