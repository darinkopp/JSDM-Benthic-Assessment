#compare parameter estimates between models
library(overlapping)

# functon for entire posterior distribution
getPost <- function(hM, parName){
  bind0 = function(...) {abind::abind(..., along = 3)}
  mpost1 <- poolMcmcChains(hM$postList)
  valList = lapply(mpost1, function(a) a[[parName]])
  val = do.call(bind0, valList)
  if (parName == "Beta") {
    colnames(val) = hM$spNames
    rownames(val) = hM$covNames
  }
  return(val)
}

#### NAP ####
# I wanted to retain NTL and CL in the model so removed SO4 and PTL
# for the reduced model. 

model.dir <- "Model_Fit_JSDM_Final/NAP/Output/Fit_Models/"

load(paste0(model.dir,"HMSC_Model_NAP_10A_1819"))
full <- model.env
load(paste0(model.dir, "HMSC_Model_NAP_10A_1819_rmSO4_PTL"))
reduced <- model.env

# evaluate whether there is a difference in model performance
#dWAIC <- abs(computeWAIC(full) - computeWAIC(reduced))

# Extract the full posterior
postBeta_Full <- getPost(full, parName = "Beta")
postBeta_Reduced <- getPost(reduced, parName = "Beta")

similarity <- data.frame(matrix(NA, nrow = dim(postBeta_Full)[2], ncol=0))
rownames(similarity) <- colnames(postBeta_Full)

supportLevel = 0.95
var <- "CL"
for (spName in colnames(postBeta_Full)){
  # spName<-colnames(postBeta_Full)[1]
  tmp <- overlapping::overlap(list(postBeta_Full[var, spName,], 
                                   postBeta_Reduced[var, spName,]))  
  
  mbeta <- mean(postBeta_Full[var, spName,])
  betaP <- mean(postBeta_Full[var, spName,] > 0)
  signFull = sign(mbeta)
  signFull = signFull * ((betaP > supportLevel) + 
                           (betaP < (1 - supportLevel)) > 0)
  
  
  mbeta <- mean(postBeta_Reduced[var, spName,])
  betaP <- mean(postBeta_Reduced[var, spName,] > 0)
  signReduced = sign(mbeta)
  signReduced = signReduced * ((betaP > supportLevel) + 
                                 (betaP < (1 - supportLevel)) > 0)
  
  similarity[spName, "OV"] <- tmp$OV
  similarity[spName, "signReduced"] <- signReduced
  similarity[spName, "signFull"] <- signFull
}
##########################
dWAIC
round(quantile(similarity$OV,probs=c(0.025,0.5,0.975)),2)
similarity[similarity$OV < 0.70,]

#### UMW ##########
# removed MSAT for reduced model and evaluated whether CL was different
model.dir <- "Model_Fit_JSDM_Final/UMW/Output/Fit_Models/"

load(paste0(model.dir,"HMSC_Model_UMW_100A_1819"))
full <- model.env

load(paste0(model.dir, "HMSC_Model_UMW_10B_1819_rmMSAT"))
reduced <- model.env

#dWAIC <- abs(computeWAIC(reduced)-computeWAIC(full))

# Extract the full posterior
postBeta_Full <- getPost(full, parName = "Beta")
postBeta_Reduced <- getPost(reduced, parName = "Beta")

similarity <- data.frame(matrix(NA, nrow = dim(postBeta_Full)[2], ncol=0))
rownames(similarity) <- colnames(postBeta_Full)

supportLevel=0.95
var <- "CL"
for (spName in colnames(postBeta_Full)){
  # spName<-colnames(postBeta_Full)[1]
  tmp <- overlapping::overlap(list(postBeta_Full[var, spName,], 
                                   postBeta_Reduced[var, spName,]))  
  
  mbeta <- mean(postBeta_Full[var, spName,])
  betaP <- mean(postBeta_Full[var, spName,] > 0)
  signFull = sign(mbeta)
  signFull = signFull * ((betaP > supportLevel) + 
                           (betaP < (1 - supportLevel)) > 0)
  
  
  mbeta <- mean(postBeta_Reduced[var, spName,])
  betaP <- mean(postBeta_Reduced[var, spName,] > 0)
  signReduced = sign(mbeta)
  signReduced = signReduced * ((betaP > supportLevel) + 
                                 (betaP < (1 - supportLevel)) > 0)
  
  similarity[spName, "OV"] <- tmp$OV
  similarity[spName, "signReduced"] <- signReduced
  similarity[spName, "signFull"] <- signFull
}
###########################
dWAIC
round(quantile(similarity$OV,probs=c(0.025,0.5,0.975)),2)
similarity[similarity$OV < 0.80,]

##### XER ############
# removed SO4 for reduced model, assessed whether CL estimates were different

# model directory
model.dir <- "Model_Fit_JSDM_Final/XER/Output/Fit_Models/"

load(paste0(model.dir,"HMSC_Model_XER_10A_1819"))
full <- model.env

load(paste0(model.dir,"HMSC_Model_XER_10A_1819_rmSO4"))
reduced <- model.env

# evaluate whether there is a difference in model performance
dWAIC <- abs(computeWAIC(full) - computeWAIC(reduced))

# Ext
postBeta_Full <- getPost(full, parName = "Beta")
postBeta_Reduced <- getPost(reduced, parName = "Beta")


similarity <- data.frame(matrix(NA, nrow = dim(postBeta_Full)[2], ncol=0))
rownames(similarity) <- colnames(postBeta_Full)

supportLevel=0.95
var <- "CL"
for (spName in colnames(postBeta_Full)){
  # spName<-colnames(postBeta_Full)[1]
  tmp <- overlapping::overlap(list(postBeta_Full[var, spName,], 
                                   postBeta_Reduced[var, spName,]))  
  
  mbeta <- mean(postBeta_Full[var, spName,])
  betaP <- mean(postBeta_Full[var, spName,] > 0)
  signFull = sign(mbeta)
  signFull = signFull * ((betaP > supportLevel) + 
                           (betaP < (1 - supportLevel)) > 0)
  
  
  mbeta <- mean(postBeta_Reduced[var, spName,])
  betaP <- mean(postBeta_Reduced[var, spName,] > 0)
  signReduced = sign(mbeta)
  signReduced = signReduced * ((betaP > supportLevel) + 
                                 (betaP < (1 - supportLevel)) > 0)
  
  similarity[spName, "OV"] <- tmp$OV
  similarity[spName, "signReduced"] <- signReduced
  similarity[spName, "signFull"] <- signFull
}
##############################
dWAIC
round(quantile(similarity$OV,probs=c(0.025,0.5,0.975)),2)
similarity[similarity$OV < 0.7,]