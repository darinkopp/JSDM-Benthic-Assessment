# Prepare HMSC Data for each region as its own .CSV
rm(list=ls())
gc()

library(Hmsc)
library(reshape2)

#This is the directory that contains fitted JSDM models
FittedModelPath <- "Model_Fit_JSDM_Final"

#vector of environmental gradient names
gradients <- c("NTL", "PTL", "CL", "SO4",
               "TMEAN_S_XXXX_PT", "PSUMPY_XXXX_PT",
               "LSUB_DMM", "W1_HALL")

#Read Environmental Data 
#####
#gis data 
env.preds <- read.csv("Data_Raw/envgis.data_12022021.csv")
names(env.preds)
#Field data chemistry and substrate
nrsa.preds <- read.csv("Data_Raw/envchmphy.data_12022021.csv")
env.preds <- merge(env.preds, nrsa.preds, all.x = T, by = "UID")
rownames(env.preds) <- as.character(env.preds$UID)

########################

#create bug data matrix
#####
bug.dat <- read.csv("Data_Raw/bug.data.Full_12022021.csv")
sum(is.na(bug.dat[,"GENUS"]) & bug.dat$OTU_NRSA!="DELETE")/nrow(bug.dat)
bug.dat <- bug.dat[!is.na(bug.dat[,"GENUS"]),]
bug.mat <- acast(bug.dat, UID ~ GENUS, 
                 value.var = "TOTAL", 
                 fill = 0, fun.aggregate = sum)

#bug.mat <- acast(bug.dat, UID ~ GENUS, 
#                value.var = "TOTAL_300", 
#               fill = 0, fun.aggregate = sum)
bug.mat[bug.mat >= 1] <- 1
######################

#taxonomy
#####
taxonomy <- unique(bug.dat[,c("PHYLUM", "CLASS", "ORDER", 
                              "FAMILY", "GENUS")])
taxonomy <- taxonomy[!duplicated(taxonomy$GENUS),]
rownames(taxonomy) <- taxonomy$GENUS
######################

#traits
#####
traits <- read.csv("Data_Raw/NRSA_benthicTaxa-OTUs_CH_7272021_ThrmOpt.csv", 
                   na.strings=c("","NA"), stringsAsFactors = F)
traits <- traits[traits$GENUS %in% colnames(bug.mat), 
                 c("GENUS", "FFG", "HABIT", "PTV", "ThrmOptV")]
rownames(traits)<-traits$GENUS
traits[,"PTV"] <- as.numeric(traits[,"PTV"])
traits[,"ThrmOptV"] <- as.numeric(traits[,"ThrmOptV"])

q <- split(traits[rownames(taxonomy),], 
           taxonomy$FAMILY)

MeanMajVals <- lapply(q, function(x){ 
  #x<-q[[4]]
  FFG <- names(which.max(table(x$FFG)))
  FFG <- ifelse(is.null(FFG), NA, FFG)
  HABIT <- names(which.max(table(x$HABIT)))
  HABIT <- ifelse(is.null(HABIT), NA, HABIT)
  PTV <- mean(x$PTV, na.rm = T)
  ThrmOptV <- mean(x$ThrmOptV, na.rm = T)
  data.frame(FFG, HABIT, PTV, ThrmOptV)})

for (x in names(q)){
  #x <- names(q)[4]
  for (i in c("FFG", "HABIT", "PTV", "ThrmOptV")){
    #i <- "FFG"
    if(any(is.na(q[[x]][,i])) & i %in% c("HABIT","FFG")){
      q[[x]][,i][is.na(q[[x]][,i])] <- as.character(MeanMajVals[[x]][,i])
    }
    
    if(any(is.na(q[[x]][,i])) & !i %in% c("HABIT","FFG")){
      q[[x]][,i][is.na(q[[x]][,i])] <- MeanMajVals[[x]][,i]
    }
  }
}

x <- do.call(rbind, q)
x <- x[!is.na(x$GENUS),]
rownames(x) <- x$GENUS
traits<-x

sel <- apply(traits, 1, function(x) all(!is.na(x)))
traits <- traits[sel,]
traits$Clinger <- factor(ifelse(traits$HABIT=="CN", "Y", "N"))
traits$Scraper <- factor(ifelse(traits$FFG=="SC", "Y", "N"))
######################

# Apply Data transformations to Environmental data
#####
ind <- env.preds$TMEAN_S_XXXX_PT == max(env.preds$TMEAN_S_XXXX_PT)
env.preds$TMEAN_S_XXXX_PT[ind] <- NA
ind <-env.preds$PSUMPY_XXXX_PT == max(env.preds$PSUMPY_XXXX_PT)
env.preds$PSUMPY_XXXX_PT[ind] <- NA

env.preds <- env.preds[complete.cases(env.preds[, gradients]),]

env.preds$SO4 <- log(env.preds$SO4)
env.preds$NTL <- log(env.preds$NTL + 1)
env.preds$PTL <- log(env.preds$PTL + 1)
env.preds$CL <- log(env.preds$CL + 1)
env.preds$W1_HALL <- log(env.preds$W1_HALL + 1)
env.preds$PSUMPY_XXXX_PT <- log(env.preds$PSUMPY_XXXX_PT)

env.preds <- env.preds[complete.cases(env.preds[, c(gradients)]),]
######################

# select 2018-2019 Visit #1 survey data for calibration and
# visit 2 for validation
#####
RS <- env.preds[env.preds$YEAR %in% c(2018, 2019),]
RS <- split(RS, as.character(RS$UNIQUE_ID))
VAL <- do.call(rbind, lapply(RS, function(x) x[x$VISIT_NO != 1,]))
rownames(VAL) <- VAL$UID
VAL$RT_MASTER <- as.character(VAL$RT_MASTER)

RS <- do.call(rbind, lapply(RS, function(x) x[x$YEAR == max(x$YEAR) & x$VISIT_NO == 1,]))
rownames(RS) <- RS$UID
RS$RT_MASTER <- as.character(RS$RT_MASTER)

RS <- split(RS, as.character(RS$ECOWSA9))
VAL <-split(VAL, VAL$ECOWSA9)
######################

sum(unlist(lapply(taxon,function(x)dim(x))))
lapply(names(RS[[1]]))

# Identify sites that have complete data
######
#select bug data associated with RS
taxon <- lapply(RS, function(x){
  #records in RS that have bug data
  wbugs <- rownames(x)[(rownames(x) %in% row.names(bug.mat))]
  bug.mat[wbugs,]})

#check all taxa (columns) have trait values
taxon <- lapply(taxon, function(x) 
  x[,intersect(colnames(x), rownames(traits))])

#Minimum Prevalence
taxon <- lapply(taxon, function(x){ 
  #remove rare taxa with traits 
  rare <- 0.1 #remove rare taxa (>10% Prevalence w/in an ecoregion)
  nonrare <- colnames(x)[colSums(x)/nrow(x) >= rare]
  z <- x[,nonrare]
  z[rowSums(z) > 0, ]})

#Select sites w/ covariates 
#select sites in bug data w/covatiates
sites <- lapply(taxon, function(x) env.preds[rownames(x),])
traits <- lapply(taxon, function(x) traits[colnames(x),])
taxonomy <- lapply(taxon, function(x) taxonomy[colnames(x),])

#visit 2 validation data
taxon.val <- lapply(VAL, function(x){
  #records in RS that have bug data
  wbugs <- rownames(x)[(rownames(x) %in% row.names(bug.mat))]
  bug.mat[wbugs,]})
#Select sites w/ covariates 
#select sites in bug data w/covatiates
sites.val <- lapply(taxon.val, function(x) env.preds[rownames(x),])
######################


#check for strong correlations between variables within regions
#######
library(HH)
#calculate VIF, all < 5 
vifTbl <- lapply(sites, function(x) round(vif(x[,gradients]),2))
vifTbl <- do.call(rbind, vifTbl)

#write.csv(vifTbl, "Figures/TableS3_JSDM_VIF.csv")
# inital submission considered values >5 to indicate strong 
# potential for multicolinearity. After peer-review we wanted to
# to check wether a stricter criteria would change the results. 
# Three regions with VIF > 3 were re run after removing the variables 
# that we not altered by human activity. See Figure 4 in the manuscript
vif(sites[["NAP"]][,gradients[-c(2,4)]])
vif(sites[["UMW"]][,gradients[-5]])
vif(sites[["XER"]][,gradients[-4]])
######################

#check matching rownames
#####
sapply(names(taxon), function(x)
  all(all(rownames(taxon[[x]])==rownames(sites[[x]])),
      all(colnames(taxon[[x]])==rownames(traits[[x]])),
      all(colnames(taxon[[x]])==rownames(taxonomy[[x]]))))
######################

#Write data files to sub-directories 
#####
for (i in names(taxon)){
  #i<-"SAP"
  #print(dim(taxon[[i]]))
  write.csv(taxon[[i]], paste0("Model_Fit_JSDM_Final/", i, "/Model_Data/Taxa_", i, ".csv"))
  write.csv(sites[[i]], paste0("Model_Fit_JSDM_Final/", i, "/Model_Data/Sites_", i, ".csv"))
  write.csv(traits[[i]], paste0("Model_Fit_JSDM_Final/", i, "/Model_Data/Traits_", i, ".csv"))
  write.csv(taxonomy[[i]], paste0("Model_Fit_JSDM_Final/", i, "/Model_Data/Taxonomy_", i, ".csv"))
  write.csv(sites.val[[i]], paste0("Model_Fit_JSDM_Final/", i, "/Model_Data/ValSites_", i, ".csv"))
  write.csv(taxon.val[[i]], paste0("Model_Fit_JSDM_Final/", i, "/Model_Data/ValTaxa_", i, ".csv"))
}
######################
Regions <- c("NAP", "SAP", "CPL", "UMW", "NPL", 
             "TPL","SPL", "XER", "WMT")



# after model fitting. 
# copy fitted models to common directory 
######
modelnames <- grep("1819", 
                   list.files(FittedModelPath, full.names = F, recursive = T), 
                   value = T)
# confirm that file names represent the model that converged for each regions subdirectory
# UMW was only model that needed to increase thinning
modelnames <- modelnames[c(1:7, 10:11)]

#copy fitted model to parent directory
for(i in modelnames){
  #i<-modelnames[1]
  nm <- unlist(lapply(strsplit(i,"/"), "[[", 4))
  file.copy(from = paste0(FittedModelPath, "/", i), 
            to = paste0(FittedModelPath, "/", nm))
}
############################




