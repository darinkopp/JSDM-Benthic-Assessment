# Fitting Random Forest Models
rm(list = ls())
library(sf)
library(terra)
library(quantregForest)
library(ggplot2)
#citation(package = "quantregForest")


# Random forest models were used to relate present-day physiochemical conditions 
# to natural and anthropogenic variables and predict values the represent a 
# removal or reduction in anthropogenic disturbance.

tmp <- read.csv("Present_Day_Environmental_Data.csv")

#Define variables used as predictors
######

# These variables were selected because of their hypothesized relationship 
# with the physiochemical response variables. Note that Legacy P, an important 
# predictor from Sabo et al. 2023, is excluded from list of candidate predictors 
# because it is strongly correlated with P_Inputs. 

# All other candidate variables had correlations <0.7 (analysis below).

varsNTL <- c("WsAreaSqKm", "ElevWs", "RunoffWs", 
             "P2O5Ws","NWs", "SWs", "KffactWs", 
             "BFIWs", "PctHbWetXXXXWs", "PermWs", 
             "PctCropXXXXWs", "RdDensWs", "N_inputs", 
             "mean.n_tw.2018", "TMEAN_S_XXXX_PT",
             "PSUMPY_XXXX_PT")

varsPTL <-  c("WsAreaSqKm", "ElevWs", "RunoffWs", "P2O5Ws", 
              "SWs", "KffactWs", "BFIWs", "PctHbWetXXXXWs", 
              "PermWs","ClayWs", "PctCropXXXXWs", "RdDensWs", 
              "P_inputs","mean.pdep2013","TMEAN_S_XXXX_PT",
              "PSUMPY_XXXX_PT")

varsCL <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "NWs", 
            "P2O5Ws", "SWs", "KffactWs", "BFIWs", 
            "PctHbWetXXXXWs", "PermWs", "ClayWs", 
            "PctCropXXXXWs", "RdDensWs", "CoalMineDensWs", 
            "MineDensWs", "mean.cl_tw.2018", "mean.s_tw.2018", 
            "TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT")

varsSO4 <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "NWs", 
             "P2O5Ws", "SWs", "KffactWs", "BFIWs", 
             "PctHbWetXXXXWs", "PermWs", "ClayWs", 
             "PctCropXXXXWs", "RdDensWs", "CoalMineDensWs", 
             "MineDensWs", "mean.s_tw.2018", "TMEAN_S_XXXX_PT", 
             "PSUMPY_XXXX_PT")

varsLSUB <- c("SandWs", "PctAlluvCoastWs", "NABD_NrmStorWs",  
              "RdDensWs", "PctCropXXXXWsRp100", "PctFstWs", 
              "CoalMineDensWs", "MineDensWs", "W1_HAG", 
              "KffactWs","WsAreaSqKm", "StreamPower","RunoffWs")


a<-cor(tmp[complete.cases(tmp[, varsSO4]),varsSO4], 
       method = "pearson")
library(HH)
for (i in list(varsLSUB,varsSO4,varsCL,varsPTL,varsNTL)){
  #i<-varsNTL
  vf<-vif(tmp[complete.cases(tmp[, i]),i])  
  print(vf[which.max(vf)])
}


apply(a,1,function(x) sum(abs(x)>0.7))
################################

#Define Model formulas
######
formulas <- list()
for (i in c("NTL", "PTL", "CL", "SO4", "LSUB_DMM")){
  if(i == "NTL"){
    formulas[[i]] <- formula(paste(i, "~", paste(c(varsNTL), collapse = " + ")))  
  } 
  if(i == "PTL"){
    formulas[[i]] <- formula(paste(i, "~", paste(c(varsPTL), collapse = " + "))) 
  }
  if(i %in% c("CL")){
    formulas[[i]] <- formula(paste(i, "~", paste(c(varsCL), collapse = " + "))) 
  }
  if(i %in% c("SO4")){
    formulas[[i]] <- formula(paste(i, "~", paste(c(varsSO4), collapse = " + "))) 
  }
  if(i %in% c("LSUB_DMM")){
    formulas[[i]] <- formula(paste(i, "~", paste(c(varsLSUB), collapse = " + ")))
  }
}

################################


# Fit Random Forest Models 
########
#We used the quantile random forest modeling because they produced prediction 
# intervals that we initially thought could aid in the assessment of the 
# potential impacts of athropogneic activities. We since favored a more 
# straightforward approach using RMSE. The predictions generated from quantile 
# RF should be analogous to the RandomForest r package commonly used. Each 
# gradient was run separately using data from all ecoregions combined. 
# 
# A 80/20 training/testing split of the data was performed prior to modeling. 

library(vegan)

allvars <- unique(unlist(lapply(formulas, function(x) all.vars(x))))
tmp <- tmp[complete.cases(tmp[, allvars]),]

# sample training and testing data consistently for all models
allrows <- 1:nrow(tmp)
trainrows <- sample(allrows, 
                    ceiling(0.80 * nrow(tmp)), 
                    replace = F)
testrows <- allrows[!allrows %in% trainrows]

# Units for response variables: NTL and PTL = ln(x+1) ug/L; CL = ln(x+1) 
# mg/L; SO4 = ln(x) mg/L; W1_HALL = ln(x+1); LSUB_DMM = log10(x) mm; 
# TMEAN_S_XXXX_PT = degC; PSUMPY_XXXX_PT = ln(x) mm


for (p in 1:length(formulas)){
  # p <- 5
  vars <- all.vars(formulas[[p]])
  response <- vars[1]
  
  tmp.cc <- tmp[, c("UID", "ECOWSA9", vars)]
 
  traindat <- tmp.cc[trainrows,]
  testdat <- tmp.cc[testrows,]
  
  # Quantile random forest
  qrf <- quantregForest(x = traindat[,vars[-1]], 
                        y = traindat[,vars[1]], 
                        ntree = 2500, 
                        importance = T)
  
  # create model fitted object as a list include training and testing data
  rf_model_list <- list(qrf = qrf, 
                        train = traindat, 
                        test = testdat)
  
  # final Models were fit and saved 4/11/2024 - Refit for revision. Ensured each 
  # model has the same calibration and validation data
  #save(rf_model_list, file = paste0("Model_Fit_RF_Final/", response, "_RF_optimize_2b"))
}
################################

# Delete this
#####
df <- t(apply(df, 1, scale))
#rf_model_list$train
rownames(imp)
# make sure the datasets have complete.cases with respect the the 
# predictor variables
newdata.id <- as.character(rf_model_list$train$UID)
train.id <- rownames(rf_model_list$train)
df <- rbind(rf_model_list$train[,rownames(imp)],
            nd[newdata.id, rownames(imp)])
df <- t(apply(df, 1, scale))


#calculate dissimilarity between each site and its hindcast condition
d <- vegdist(df, method = "euclidean")
d <- as.matrix(d)
trainMean <- mean(d[train.id, train.id])
round(min(d[train.id, site])/trainMean,3)


#extract values from the dissimilarity matrix
diff.tbl <- data.frame(matrix(NA, 0, 2))
for (i in 1:length(newdata.id)){
  #i<-1
  site <- newdata.id[i]
  #if(nchar(site) == 7){
    diff.tbl[i, 1] <- site
    diff.tbl[i, 2] <- round(min(d[train.id, site])/trainMean,3)
    #diff.tbl[i, 2] <- d[site, paste0(site, "1")] 
  }
#}
windows()
quantile(diff.tbl$X2, probs = c(0.025, 0.975))

#add ecoregions to stratify validation
dissimilar <- data.frame(diff.tbl, AG_ECO9 = pd[diff.tbl[,1], "AG_ECO9"])
names(dissimilar)
dissimilar
validation.2 <- split(dissimilar, dissimilar$AG_ECO9)


validation.2 <- lapply(validation.2, function(x) x[order(x[,2])[1:10],])
validation.2 <- lapply(validation.2, function(x) 
  x[sample(1:nrow(x), size = 5, replace = F),])
validation.2 <- do.call(rbind, validation.2)
validation.2 <- tmp[tmp$UID%in%validation.2[,1],]

# remove these variables from the dataset
tmp <- tmp[!tmp$UID%in%validation.2[,1],]
#####################################

# Delete This too
#######
library("quantregForest")
dim(getTree(rf_model_list$qrf, k=1, labelVar = T))
rf_mod <- rf_model_list$qrf
NewData <- rf_model_list[["validation"]]
NewData <- NewData[complete.cases(NewData),]
test.x <- predict(rf_mod, 
                  newdata = NewData, 
                  what = c(0.05, 0.5, 0.95))
test.y <- NewData[,"PTL"]

cor(test.y,test.x)^2
# diagnostics for testing data
vldRsq <- round(1 - (sum((test.y - test.x[,2])^2)/
                       sum((test.y - mean(test.y))^2)), 2) 
vldRMSE <- round(sqrt(mean((test.y - test.x[,2])^2)),2)

plot(test.y~test.x[,2])
summary(lm(test.y~test.x[,2]))
# find representative tree
library(reprtree)
load(paste0("Model_Fit_RF_Final/", "PTL", "_RF_optimize_4"))
rf <- rf_model_list$qrf
nd <- rf_model_list$train

# code would not work with quantile RF. Pulled from ReprTree
preds <- reprtree:::predict2.randomForest(object=rf, newdata = nd, predict.all = T)
d <- reprtree:::dist.fn(t(preds$individual), method = "euclidean")
D <- colMeans(d)
index <- which(D == min(D))
trees <- lapply(as.list(index), function(i) getTree(rf, i, labelVar = TRUE))
names(trees) <- as.character(index)

#remove spaces in a word
collapse <- function(x){
  x<-sub(" ","_",x)
  
  return(x)
}
names(trees[[1]]) <- sapply(names(trees[[1]]),collapse)

# functions to get rules from the tree
getConds2 <- function(tree){
  #store all conditions into a list
  #conds <- list()
  edge_data <- data.frame()
  #start by the terminal nodes and find previous conditions
  id.leafs <- which(tree$status == -1)
  j <- 0
  for(i in id.leafs){
    #i<-id.leafs[1]
    j <- j+1
    prevConds <- prevCond(tree,i)
    conds[[j]] <- prevConds$cond
    edge_data[j, "to"] <- i
    edge_data[j, "from"] <- prevConds$id
    edge_data[j, "split"] <- "terminal"
    edge_data[j, "path"] <- i
    edge_data[j, "prediction"]<-tree$prediction[i]
    while(prevConds$id > 1){
      j <- j + 1
      edge_data[j,"to"] <- prevConds$id
      prevConds <- prevCond(tree, prevConds$id)
      #conds[[j]] <- paste(conds[[j]]," & ", prevConds$cond)
      edge_data[j,"from"] <- prevConds$id 
      edge_data[j,"split"] <- prevConds$cond
      edge_data[j, "path"] <- i
      edge_data[j, "prediction"] <- tree$prediction[i]
    }
    if(prevConds$id == 1){
      #conds[[j]]<-paste(conds[[j]]," => ",tree$prediction[i])
    }
    
  }
  
  return(edge_data)
  #return(conds)
}
prevCond <- function(tree,i){
  if(i %in% tree$right_daughter){
    # GREATER than split point--- This will be the case for most anthropogenic 
    # variables
    id <- which(tree$right_daughter == i)
    cond <- paste(tree$split_var[id], ">", round(tree$split_point[id],2))
  }
  
  if(i %in% tree$left_daughter){
    # LESS than the split point --- this might only be relevant for natural 
    # vegetation cover
    id<-which(tree$left_daughter==i)
    cond<-paste(tree$split_var[id],"<",round(tree$split_point[id],2))
  }
  
  return(list(cond=cond,id=id))
}


edge_data <- getConds2(trees[[1]])
edge.criteria <- do.call(rbind, strsplit(edge_data$split, " "))
colnames(edge.criteria)<-c("Variable", "Sign", "Value")
edge_data <- data.frame(edge_data, edge.criteria)
edge_data[edge_data$Sign == "terminal", "Value"] <- edge_data$prediction[edge_data$Sign == "terminal"]
edge_data$Value <- as.numeric(edge_data$Value)
edge.split <- split(edge_data, edge_data$path)

edge.split[[34]]
edge_data$Value[edge_data$Sign == "terminal"]

# remove edges manage.vars
rm.edge <- lapply(edge.split, function(a) {
  #a <- edge.split[["1017"]]
  any(a$Variable %in% c("PctCropXXXXWs", 
                        "RdDensWs", "CoalMineDensWs",
                        "MineDensWs", "W1_HAG",
                        "W1_HNOAG", "W1H_CROP",
                        "PctCropXXXXWsRp100") & 
        a$Sign == ">" & a$Value > 0 | 
        a$Variable == "PctFstWs" & a$Sign == "<" & a$Value > 100)
  })

edge_data[edge_data$to==42,]

edge.split[unlist(rm.edge)]               
edge.split <- do.call(rbind, edge.split[!unlist(rm.edge)])
edge.split[!unlist(rm.edge)]
manage.vars$manage.vars
lapply(edge.split[unlist(rm.edge)],3)
library(igraph)
library(ggraph)
edge.split$Value

edge_data <- unique(edge.split)
edge_data[edge_data$Variable%in%c("PctCropXXXXWs", 
              "RdDensWs", "CoalMineDensWs",
              "MineDensWs", "W1_HAG",
              "W1_HNOAG", "W1H_CROP",
              "PctCropXXXXWsRp100","PctFstWs"), "Disturbance"]<- 0.05
edge_data[edge_data$Disturbance==0.025,]
edge_data[!edge_data$Variable%in%c("PctCropXXXXWs", 
                                  "RdDensWs", "CoalMineDensWs",
                                  "MineDensWs", "W1_HAG",
                                  "W1_HNOAG", "W1H_CROP",
                                  "PctCropXXXXWsRp100", 
                                  "PctFstWs"), "Disturbance"]<- 0.025
edge_data$Disturbance<-as.factor(edge_data$Disturbance)
edge_data[edge_data$Sign!="terminal", "Value"] <- NA
edge_data[edge_data$Sign!="terminal", "W"] <- 0.25
edge_data[edge_data$Sign=="terminal", "W"] <- 1

ed <- edge_data[,c("from", "to", "path", "split", "Value", "Disturbance")]
ed <- edge_data[,c("from", "to", "Disturbance")]
names(ed)
ed[ed$to==30,]
ed[ed$path==30,]
mygraph <- ?graph_from_data_frame(unique(ed), directed = T)
z<-as_data_frame(mygraph, what=c("both"))
z$vertices[1:6,]
z$edges[1:6,]
as_data_frame(mygraph, what="edges")
V(mygraph)$type <- bipartite_mapping(mygraph)$type
labs <- V(mygraph)
E(mygraph)
g1 <- graph_from_adjacency_matrix(mygraph, add.colnames='path')
ggraph::
windows()
ss <- ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  #geom_edge_diagonal() +
  #geom_edge_parallel() +
  geom_edge_elbow(aes(color = Disturbance, width=Disturbance))+
  scale_edge_width_manual(values = c(0.25, 0.5))+
  geom_node_point(aes(filter = leaf == T))+
  geom_node_label(aes(label = labs), label.size = 0.001)

  geom_node_text(aes(label = ifelse(leaf, path, NA)), size = 3, repel = TRUE)
  theme_void()

# so this tree represents the predicted values 
print(ss)
#########################################

################################################################################

# Predict Nutrients, Salinity, PHAB and Climate under new conditions

################################################################################

# Initialize datasets 
#########
# Present day dataset need to calculate difference of hindcasted predictions 
pd <- read.csv("Present_Day_Environmental_Data.csv")
pd <- reshape2::melt(pd[,c("UID","ECOWSA9","LAT_DD","LON_DD",
                           "NTL","PTL","SO4","CL",
                           "LSUB_DMM","W1_HALL",
                           "PSUMPY_XXXX_PT","TMEAN_S_XXXX_PT")], 
                     value.name = "PD", 
                     id.vars = c("UID","ECOWSA9","LAT_DD","LON_DD"))

# data with anthropogenic factors removed
NewPredDat <- read.csv("Model_Fit_RF_Final/NoHumanPrediction.csv")
rownames(NewPredDat) <- NewPredDat$UID

#######################

# Random forest models names
files <- list.files("Model_Fit_RF_Final/", 
                    full.names = T)
optimize_rf <- grep("RF_optimize_2b", files, value = T)
######################################

# make predictions of data
######

# initialize new prediction dataset in long format
newJSDMpreds <- data.frame() 

for (i in optimize_rf){
  #i <- optimize_rf[2]
  load(i)
  r <- unlist(lapply(strsplit(i,"/"), "[[", 2))
  r <- unlist(lapply(strsplit(r, "_"), "[[",1))
  
  r <- ifelse(r=="LSUB", paste0(r,"_DMM"),r)
  
  # take rownames from fitted model
  PredDat <- rbind(rf_model_list[["train"]], rf_model_list[["test"]])
  nd <- NewPredDat[as.character(PredDat$UID), c(r, rownames(rf_model_list[["qrf"]]$importance))]
  nd <- nd[complete.cases(nd),]
  
  # Make predictons with new data
  p <- predict(rf_model_list[["qrf"]], 
               newdata = nd, 
               what = c(0.025, 0.5, 0.975))
  p <- setNames(data.frame(p), c("PI_Lwr", "HC", "PI_Upr"))
  rmse <- sqrt(tail(rf_model_list[["qrf"]]$mse, 1))
  
  newJSDMpreds <- rbind(newJSDMpreds, 
                        data.frame(variable = r,
                                   UID = rownames(nd),                                 
                                   HC = p$HC,
                                   rmse))
}

# Temp from 1940-1950 already calculated in NewD
# extract standard deviation from 10yr average
TMEAN_S_0050_sd <- terra::rast("Data_Raw/TMEAN_S_0050_sd.grd")
q <- NewPredDat %>%
  st_as_sf(coords = c("LON_DD","LAT_DD"), 
           crs = 4269) %>% 
  st_transform(crs = crs(TMEAN_S_0050_sd)) %>%
  st_coordinates() %>%
  data.frame() %>%
  terra::extract(TMEAN_S_0050_sd, y = .)

newJSDMpreds <- rbind(newJSDMpreds, 
                      data.frame(UID = NewPredDat$UID,
                                 variable = "TMEAN_S_XXXX_PT",
                                 HC = NewPredDat$TMEAN_S_XXXX_PT,
                                 rmse = q[,2]))


# precipitation from 1900-1950
# extract standard deviation from 10yr average
PSUM_0050_sd <- terra::rast("Data_Raw/PSUM_0050_sd.grd")
q <- NewPredDat %>%
  st_as_sf(coords = c("LON_DD","LAT_DD"), 
           crs = 4269) %>% 
  st_transform(crs = crs(PSUM_0050_sd)) %>%
  st_coordinates() %>%
  data.frame() %>%
  terra::extract(PSUM_0050_sd, y = .) 

#the SD is in regular units
newJSDMpreds <- rbind(newJSDMpreds, 
                      data.frame(UID = NewPredDat$UID,
                                 variable = "PSUMPY_XXXX_PT",
                                 HC = NewPredDat$PSUMPY_XXXX_PT,
                                 rmse = q[,2]))

# W1_HALL
# 0.33 is the threshold used for riparian disturbance in NRSA Phab
newJSDMpreds <- rbind(newJSDMpreds,
                      data.frame(UID = NewPredDat$UID,
                                 variable = "W1_HALL",
                                 HC = log(0.33 + 1),
                                 rmse = NA))

newJSDMpreds <- merge(pd, newJSDMpreds, by = c("UID", "variable"), all.x = TRUE)
######################################


# Identify sites that exceed to standard deviation thresholds
#######


# Initialize column with present day value
# Values that are not 
newJSDMpreds$newValue <- newJSDMpreds$PD
#calculate difference in units of standard deviations
newJSDMpreds$diff_sd <- (newJSDMpreds$HC - newJSDMpreds$PD)/newJSDMpreds$rmse

#need to calculate Precip with backtransformed values such that SD units are the same
pptIndex <- newJSDMpreds$variable=="PSUMPY_XXXX_PT"
newJSDMpreds[pptIndex, "diff_sd"] <-
  (exp(newJSDMpreds$HC[pptIndex]) - exp(newJSDMpreds$PD[pptIndex]))/newJSDMpreds$rmse[pptIndex]

# Identify values that exceed 2 rmse -- W1_HALL does not have RMSE
# for chemistry variables only sites with present-day values higher the hindcasted
newJSDMpreds$newValue <- ifelse(newJSDMpreds$variable %in% c("NTL", "PTL", "SO4", "CL")&
                                  !is.na(newJSDMpreds$rmse) & 
                                  newJSDMpreds$diff_sd <= -2 ,
                                newJSDMpreds$HC, newJSDMpreds$newValue)

# Identify substrate values that have differences >2 sd
# substrate can be coarser or finer, climate warmer/cooler wetter/drier
# newJSDMpreds$newValue <- ifelse(newJSDMpreds$variable %in% c("LSUB_DMM","TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT") & 
#                                   !is.na(newJSDMpreds$rmse) &
#                                   abs(newJSDMpreds$diff_sd) >= 2,
#                                 newJSDMpreds$HC, newJSDMpreds$newValue)

newJSDMpreds$newValue <- ifelse(newJSDMpreds$variable %in% c("LSUB_DMM") & 
                                  !is.na(newJSDMpreds$rmse) &
                                  abs(newJSDMpreds$diff_sd) >= 2,
                                newJSDMpreds$HC, newJSDMpreds$newValue)

newJSDMpreds$newValue <- ifelse(newJSDMpreds$variable %in% c("TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT") & 
                                  !is.na(newJSDMpreds$rmse) &
                                  abs(newJSDMpreds$diff_sd) >= 2,
                                newJSDMpreds$HC, newJSDMpreds$newValue)

# Identify W1_Hall values that exceed threshold
newJSDMpreds$newValue <- ifelse(newJSDMpreds$variable == "W1_HALL" &
                                  newJSDMpreds$PD >= newJSDMpreds$HC,
                                newJSDMpreds$HC, newJSDMpreds$newValue)

# Add indicator of exceedance
newJSDMpreds$Exceed <- ifelse(newJSDMpreds$PD == newJSDMpreds$newValue, "N", "Y")

#order and rename levelsfor plots 
newJSDMpreds[,"variable"] <- factor(newJSDMpreds[,"variable"],
                                    levels = c("NTL","PTL", "CL","SO4",
                                               "LSUB_DMM", "W1_HALL", 
                                               "TMEAN_S_XXXX_PT", 
                                               "PSUMPY_XXXX_PT"))

newJSDMpreds$titles <- newJSDMpreds[,"variable"]
levels(newJSDMpreds$titles) <- c("NTL","PTL", "CL","SO4", 
                                 "SUBD","RPDI","MSAT","TPRCP") 
######################################

# Add back transformations
######
# add backtransformed data field for plotting 
newJSDMpreds$HC.bktrans <- NA
newJSDMpreds$PD.bktrans <- NA

unique(newJSDMpreds$variable)
# units for SO4 "mg/L" and total precip "mm" were ln(x) transformed for JSDM
# Precipitation does not have RMSE
newJSDMpreds[newJSDMpreds$variable %in% c("SO4", "PSUMPY_XXXX_PT"), "HC.bktrans"] <- 
  exp(newJSDMpreds$HC[newJSDMpreds$variable %in% c("SO4", "PSUMPY_XXXX_PT")])
newJSDMpreds[newJSDMpreds$variable %in% c("SO4", "PSUMPY_XXXX_PT"), "PD.bktrans"] <- 
  exp(newJSDMpreds$PD[newJSDMpreds$variable %in% c("SO4", "PSUMPY_XXXX_PT")])

# Units for NTL "ug/L", PTL "ug/L", CL "mg/L", and W1_HALL "unitless" were ln(x+1) transformed for JSDM
newJSDMpreds[newJSDMpreds$variable %in% c("PTL", "NTL", "CL"), "HC.bktrans"] <- 
  exp(newJSDMpreds$HC[newJSDMpreds$variable %in% c("PTL", "NTL", "CL")]) - 1
newJSDMpreds[newJSDMpreds$variable %in% c("PTL", "NTL", "CL"), "PD.bktrans"] <- 
  exp(newJSDMpreds$PD[newJSDMpreds$variable %in% c("PTL", "NTL", "CL")]) - 1 

# W1_HALL does not have RMSE but uses the 0.33 threshold
newJSDMpreds[newJSDMpreds$variable %in% c("W1_HALL"), "HC.bktrans"] <- 
  exp(newJSDMpreds$newValue[newJSDMpreds$variable %in% c("W1_HALL")]) - 1
newJSDMpreds[newJSDMpreds$variable %in% c("W1_HALL"), "PD.bktrans"] <- 
  exp(newJSDMpreds$PD[newJSDMpreds$variable %in% c("W1_HALL")]) - 1 


# Units for LSUB_DMM "mm" were Log10(x) transformed
newJSDMpreds[newJSDMpreds$variable %in% c("LSUB_DMM"), "HC.bktrans"] <- 
  10^(newJSDMpreds$HC[newJSDMpreds$variable %in% c("LSUB_DMM")]) 
newJSDMpreds[newJSDMpreds$variable %in% c("LSUB_DMM"), "PD.bktrans"] <- 
  10^(newJSDMpreds$PD[newJSDMpreds$variable %in% c("LSUB_DMM")]) 

# Units for MSAT "degC" were not transformed
newJSDMpreds[newJSDMpreds$variable %in% c("TMEAN_S_XXXX_PT"), "HC.bktrans"] <- 
  newJSDMpreds$HC[newJSDMpreds$variable %in% c("TMEAN_S_XXXX_PT")] 
newJSDMpreds[newJSDMpreds$variable %in% c("TMEAN_S_XXXX_PT"), "PD.bktrans"] <- 
  newJSDMpreds$PD[newJSDMpreds$variable %in% c("TMEAN_S_XXXX_PT")] 

##############################
names(newJSDMpreds)[3] <- "Ecoregion"

write.csv(newJSDMpreds, "Hindcast_Environmental_Data_r1.csv")

# delete 
#######
z<-read.csv("Hindcast_Environmental_Data_r1.csv")
quantile(z[z$variable=="TMEAN_S_XXXX_PT","HC"]-
            z[z$variable=="TMEAN_S_XXXX_PT","PD"], probs = c(0.05,0.95))
z[z$variable=="PSUMPY_XXXX_PT","HC"]
#####################

# this should have been moved to the main analysis document... 
# check and remove variables to avoid potential risks associated with 
# extrapolation
######################
is_contained <- function(small_vec, large_vec) {
  all(small_vec >= min(large_vec) & small_vec <= max(large_vec))
  }

# Example usage:
vec1 <- c(2, 5, 7,9)
vec2 <- c(1, 3, 6, 8)


# Adjust for extrapolation, regional models should not be extrapolated, 
# beyond the conditions observed under present conditions. 
# list of fitted model names (see Kopp et al 2023 for details)

FittedModelPath <- "Model_Fit_JSDM_Final"
modelnames <- grep("1819", list.files(FittedModelPath, full.names = F), value = T)

Dat <- read.csv("Hindcast_Environmental_Data.csv")

#create wide format (i.e. site x environmental variable)
NewPredDat <- reshape2::dcast(UID + Ecoregion + LAT_DD + LON_DD ~ variable, 
                              data = Dat, value.var = "newValue")
rownames(NewPredDat) <- NewPredDat$UID


out <- data.frame()
for (mn in modelnames){
  # mn <- modelnames[1]
  r <- unlist(lapply(strsplit(mn, "_"), "[", 3))
  load(paste0(FittedModelPath, "/",  mn))
  tmp <- setNames(data.frame(matrix(NA, 
                                    nrow = 1, 
                                    ncol = dim(model.env$XData)[2])),
                  names(model.env$XData))
  
  for (var in names(model.env$XData)){
    #var <- names(model.env$XData)[1]
    
    vec1 <- model.env$XData[,var]
    #hindcast Conditions 
    vec2 <- NewPredDat[rownames(model.env$XData),var]
    
    # ask whether the range of hindcast values are 
    # within the range of the present day values 
    TF <- is_contained(vec2, vec1)
    tmp[1, var] <- TF
    rownames(tmp) <- r
    
    if(!TF){
      print(sum(vec2 < min(vec1)|vec2 > max(vec1)))
      #vec2[which(vec2 < min(vec1))]
      #vex2[which(vec2 > max(vec1))]
      
      NewPredDat[rownames(model.env$XData)[which(vec2 < min(vec1))],var] <- min(vec1)
      NewPredDat[rownames(model.env$XData)[which(vec2 > max(vec1))],var] <- min(vec1)
      
      #stop()
    }
    
  }
  out <- rbind(out,tmp)
}

write.csv(NewPredDat, "Hindcast_Environmental_Data_NoExtrapolation.csv")

#tmean summer
21/nrow(NewPredDat)
