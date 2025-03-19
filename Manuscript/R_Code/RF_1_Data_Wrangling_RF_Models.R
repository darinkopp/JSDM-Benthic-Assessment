# Data wrangling for random forest models 
# Creates present-day data set. 


#This is the directory that contains fitted JSDM models
FittedModelPath <- "Model_Fit_JSDM_Final"

# Reads and combines data from JSDM models
##########
#create master site table. 
SiteDataFiles <- grep("Sites", 
                      list.files(FittedModelPath, full.names = F, recursive = T), 
                      value = T)

#copy fitted model to parent directory
SiteData <- data.frame()
for(i in SiteDataFiles){
  # i <- SiteDataFiles[1]
  tmp <- read.csv(paste0(FittedModelPath, "/", i),row.names = "X")
  SiteData <- rbind(SiteData, tmp)
}
################################

# Sum natural cover types 
########
FstType <- c("PctConifXXXXWsRp100", "PctDecidXXXXWsRp100", "PctMxFstXXXXWsRp100",
             "PctGrsXXXXWsRp100", "PctShrbXXXXWsRp100")
SiteData$PctFstWs <- apply(SiteData[,FstType], 1, sum)
################################

# Atmospheric Deposition
######
NADP <- read.csv("Data_Raw/Atmospheric_DepositionNRSA1819.csv")
NADP <- NADP[!duplicated(NADP$UNIQUE_ID),]
SiteData <- merge(SiteData, NADP, by = "UNIQUE_ID", all.x = T)
################################

# Nutrient inventory data 
##########
# Remove duplicated COMIDs for merging
NutInvDat <- read.csv("Data_Raw/Downscaled_Final_Cat_Ws_AddComids.csv", header=T, row.names = "X")
NutInvDatWs <- NutInvDat[!duplicated(NutInvDat$COMID), ]
NutInvDatWs <- NutInvDatWs[, grep("COMID|Ws$", names(NutInvDatWs))]

SiteData <- merge(SiteData, NutInvDatWs, by = "COMID", all.x = T)

SiteData$N_inputs <- apply(SiteData[,c("N_Fert_FarmWs", "N_Fert_UrbanWs", 
                                       "N_Human_WasteWs", "N_Livestock_WasteWs")], 1, sum)

SiteData$P_inputs <- apply(SiteData[,c("P_f_fertilizerWs",
                                       "P_human_wasteWs", 
                                       "P_livestock_WasteWs", 
                                       "P_nf_fertilizerWs")], 1, sum)
SiteData$Legacy_PWs <- ifelse(SiteData$Legacy_PWs < 0, 0, SiteData$Legacy_PWs)
################################

# Slope values from NHD; NHDplusV2 Data stored locally.
########
# #iterate through each vpu
# to extract value for a COMID

nhd.dir <- "L:/Public/dkopp/NHDPLUSV2"
SiteData$VPU <- ifelse(nchar(SiteData$VPU) == 1, paste0("0", SiteData$VPU), SiteData$VPU)
q <- split(SiteData, SiteData$VPU)

out <- data.frame()
for (i in 1:length(q)){
  #i = 12
  VPU = unique(q[[i]]$VPU) 
  
  #read NHDPLUS Slope values
  slope <- paste0("NHDPlus", VPU,"/NHDPlusAttributes/elevslope.dbf")
  slope <- grep(slope, list.files(nhd.dir, 
                                  full.names = T, 
                                  recursive=T),
                value = T)
  slope = foreign::read.dbf(slope)
  
  out <- rbind(out, slope[slope$COMID%in%q[[i]]$COMID, c("COMID","SLOPE")])
}

SiteData <- merge(SiteData, out, by = "COMID", all.x = T)
################################

# Additional phab metrics proviced by PRK
##############
Phab_Thresholds <- read.csv("Data_Raw/Phab_OE1819_PK_04062022.csv")
LSD <- read.csv("Data_Raw/NRSA18_19_FINAL_TABLE.csv")

vnames <- names(Phab_Thresholds)[!names(Phab_Thresholds) %in% names(SiteData)]
SiteData <- merge(SiteData, 
                  Phab_Thresholds[, c(vnames, "SITE_ID", "VISIT_NO", "YEAR")], 
                  by = c("SITE_ID", "VISIT_NO", "YEAR"), all.x = T)

#((LSD$PSUMPY_2019Ws - LSD$ETWs)/1000)*(1000000/31536000)*AREAWSkm2_use
LSD <- LSD[,c("UNIQUE_ID", "PSUMPY_2019Ws", "ETWs")]
LSD <- merge(SiteData[c("UNIQUE_ID", "WsAreaSqKm")], LSD, by = "UNIQUE_ID")
LSD$StreamPower <- (((LSD$PSUMPY_2019Ws - LSD$ETWs)/1000) * (1000000/31536000)*LSD$WsAreaSqKm)^0.5
SiteData <- merge(SiteData, LSD[,-2], by = "UNIQUE_ID")
SiteData$StreamPower <- SiteData$StreamPower*SiteData$SLOPE


p <- read.csv("Data_Raw/nrsa_1819_physical_habitat_larger_set_of_metrics_-_data.csv")
SiteData <- merge(SiteData, p[,c("UID", "LDMB_BW5", "LRBS_BW5")], all.x = T, by = "UID")

#some missing values for critical diameter with simple algebra. Cleared w/ PK
Dcbf_cl <- log10((10^(xdata$LSUB_DMM))/(10^xdata$LRBS_use))
SiteData$LDCBF_G08 <- ifelse(is.na(SiteData$LDCBF_G08), Dcbf_cl, SiteData$LDCBF_G08)
# a few metrics for
SiteData[is.na(SiteData$LDMB_BW5), "LDMB_BW5"] <- SiteData$LDCBF_G08[is.na(SiteData$LDMB_BW5)]
################################

# Add ecoregions to sample sites
#############
eco_l3 <- read_sf("Data_Raw/us_eco_l3.shp")
# Assign lower ecoregions
pts <- SiteData %>%
  st_as_sf(coords = c("LON_DD","LAT_DD"), crs = 4269) %>%
  st_transform(crs = st_crs(eco_l3))

z <- st_join(eco_l3, pts, join = st_contains)
z <- unique(z[,c("UID", "US_L3CODE", "US_L3NAME")])
z <- st_set_geometry((z[!is.na(z$UID),]),NULL)

SiteData <- merge(SiteData, z, by = "UID", all.x = T)
SiteData$PCT_FNSN <- SiteData$PCT_FN + SiteData$PCT_SA
###############################

# add weights 
################
weights <- read.csv("Data_Raw/nrsa_1819_site_information_-_data.csv", stringsAsFactors = F)
weights <- weights[weights$WGT_TP_CORE != 0 & weights$VISIT_NO == 1,]


SiteData <- merge(weights[,c("UID", "AG_ECO9", "WGT_TP_CORE")], SiteData, by = "UID", all.y = T)
###############################

#write.csv(SiteData, paste0("Present_Day_Environmental_Data.csv"), row.names = F)



################################################################################

# Creating dataset with anthropogenic disturbance removed/reduced

################################################################################
pd <- read.csv("Present_Day_Environmental_Data.csv")
NewPredDat <- read.csv("Present_Day_Environmental_Data.csv")

# manually identify management variables
#########
# unique(c(varsNTL, varsPTL, varsCL, varsSO4, varsLSUB))
manage.vars <- c("TMEAN_S_XXXX_PT", "PSUMPY_XXXX_PT",
                 "PctCropXXXXWs", "RdDensWs", 
                 "N_inputs", "P_inputs",
                 "CoalMineDensWs", "MineDensWs", 
                 "mean.n_tw.2018", "mean.s_tw.2018",
                 "W1_HAG", "W1_HNOAG", "W1H_CROP",
                 "PctCropXXXXWsRp100","PctFstWs") 

manage.vars.names <- c("MSAT", "TPRCP", "PctCrop", "RdDens",
                       "N_input", "P_input", "CoalMineDen", 
                       "MineDen", "N_dep", "S_dep","W1_HAG",
                       "W1_HNOAG","W1H_CROP", 
                       "PctCropRP", "PctNatRP")

manage.vars <- data.frame(manage.vars, manage.vars.names)
rownames(manage.vars) <- manage.vars[,1]
######################################

rownames(NewPredDat) <- NewPredDat$UID
mv <- manage.vars[,1][-c(1,2)]


#remove/reduce anthropogenic disturbance 
#######
# variables are not set to zero because they are not direct measures of disturbance
mv <- mv[!mv %in% c("mean.s_tw.2018", "mean.n_tw.2018", 
                    "TMEAN_S_XXXX_PT", "PSUMPY_XXXX_PT",
                    "PctFstWs")]

# climate variable 1900-1950
TMEAN_S_0050 <- terra::rast("Data_Raw/TMEAN_S_0050.grd")
q <- NewPredDat %>%
  st_as_sf(coords = c("LON_DD","LAT_DD"), 
           crs = 4269) %>% 
  st_transform(crs = crs(TMEAN_S_0050)) %>%
  st_coordinates() %>%
  data.frame() %>%
  terra::extract(TMEAN_S_0050, y = .)

NewPredDat[, "TMEAN_S_XXXX_PT"] <- q[,2]


PSUM_0050 <- terra::rast("Data_Raw/PSUM_0050.grd")
q <- NewPredDat %>%
  st_as_sf(coords = c("LON_DD","LAT_DD"), 
           crs = 4269) %>% 
  st_transform(crs = crs(PSUM_0050)) %>%
  st_coordinates() %>%
  data.frame() %>%
  terra::extract(PSUM_0050, y = .) 

NewPredDat$PSUMPY_XXXX_PT <- log(q[,2])

# Nutrients -- constant deposition
# https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/eap.1703 used 0.4kgN/ha; 
# No sites are below this level. In discussion of paper, Clark states that before 1900 could be roughly 3-5KgN/ha in the east
# and this is also consistent with some critical loads. . 0.1kgN/ha was used by Clark et al for S; 
# No sites sampled in this study are below this level. Other efforts have estimated pre-industrial S deposition 
# at 0.32–2.98 kg S·ha−1·yr−1 (Granat et al. 1976, Fakhraei et al. 2016) - used middle number
NewPredDat[,"mean.n_tw.2018"] <- ifelse(NewPredDat[,"mean.n_tw.2018"] > 5, 
                                        5, NewPredDat[,"mean.n_tw.2018"])

NewPredDat[,"mean.s_tw.2018"]<- ifelse(NewPredDat[,"mean.s_tw.2018"] > 1.65, 
                                       1.65, NewPredDat[,"mean.s_tw.2018"])

# Natural vegetation 
NewPredDat[,"PctFstWs"] <- 100

# variables directly related to human disturbance
NewPredDat[,names(NewPredDat) %in% mv] <- 0
######################################

write.csv(NewPredDat, "Model_Fit_RF_Final/NoHumanPrediction_r1.csv")



##################
# Extra
##################

#a start to removing excess rows

dim(SiteData)

# Define variable names for each response variable, 
#####
varsNTL <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "P2O5Ws","NWs", "SWs", 
             "KffactWs", "BFIWs", "PctHbWetXXXXWs", "PermWs",
             "PctCropXXXXWs", "RdDensWs", "N_inputs", "mean.n_tw.2018",
             "TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT")

varsPTL <-  c("WsAreaSqKm", "ElevWs", "RunoffWs", "P2O5Ws", "SWs", "KffactWs",
              "BFIWs", "PctHbWetXXXXWs", "PermWs","ClayWs", "PctCropXXXXWs", 
              "RdDensWs", "P_inputs","mean.pdep2013","TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT")

varsCL <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "NWs", "P2O5Ws", "SWs", "KffactWs", 
            "BFIWs", "PctHbWetXXXXWs", "PermWs", "ClayWs", "PctCropXXXXWs", "RdDensWs",
            "CoalMineDensWs", "MineDensWs", "mean.cl_tw.2018", "mean.s_tw.2018",
            "TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT")

varsSO4 <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "NWs", "P2O5Ws", "SWs", "KffactWs", 
             "BFIWs", "PctHbWetXXXXWs", "PermWs", "ClayWs", "PctCropXXXXWs", "RdDensWs",
             "CoalMineDensWs", "MineDensWs", "mean.s_tw.2018", "TMEAN_S_XXXX_PT", "PSUMPY_XXXX_PT")

#varsLSUB <- c("PSUMPY_XXXX_PT","TMEAN_S_XXXX_PT", "SandWs","PctAlluvCoastWs",
#             "NABD_NrmStorWs", "PctFstWs", "RdDensWs","PctCropXXXXWs",
#            "CoalMineDensWs", "MineDensWs", "W1_HAG", "W1_HNOAG", "W1H_CROP",
#           "KffactWs","WsAreaSqKm", "ElevWs", "RunoffWs","RckDepWs", "SLOPE","strPow")

varsLSUB <- c("SandWs", "PctAlluvCoastWs", "NABD_NrmStorWs",  "RdDensWs",
              "PctCropXXXXWsRp100", "PctFstWs",
              "CoalMineDensWs", "MineDensWs", "W1_HAG", 
              "KffactWs","WsAreaSqKm", "ElevWs", "RckDepWs", "strPow","RunoffWs")
#"PctFstWs",
#xdata$NABD_NrmStorWs2<-xdata$NABD_NrmStorWs/xdata$RunoffWs
#varsLDMB <- c("PSUMPY_XXXX_PT","TMEAN_S_XXXX_PT", "SandWs","PctAlluvCoastWs",
#             "NABD_NrmStorWs", "PctFstWs", "RdDensWs","PctCropXXXXWs",
#            "CoalMineDensWs", "MineDensWs", "W1_HAG", "W1_HNOAG", "W1H_CROP",
#            "KffactWs","WsAreaSqKm", "ElevWs", "SLOPE", "RunoffWs","RckDepWs")

#varsFNSA <- c("PSUMPY_XXXX_PT","TMEAN_S_XXXX_PT", "SandWs","PctAlluvCoastWs",
#             "NABD_NrmStorWs", "PctFstWs", "RdDensWs","PctCropXXXXWs",
#            "CoalMineDensWs", "MineDensWs", "W1_HAG", "W1_HNOAG", "W1H_CROP",
#           "KffactWs","WsAreaSqKm", "ElevWs", "SLOPE", "RunoffWs","RckDepWs")

#add PCT_FN to match NRSA thresholds
varsALL <- unique(c("UID", "UNIQUE_ID", "SITE_ID", "VISIT_NO", "YEAR", "COMID",
                    "LAT_DD", "LON_DD", "ECOWSA9","NTL", "PTL", "CL", "SO4", 
                    "LSUB_DMM","W1_HALL", varsLSUB, varsCL, varsSO4, varsPTL, 
                    varsNTL))
#############################################

xdata <- xdata[,varsALL]







dim(SiteData)

# Define variable names for each response variable, 
#####
varsNTL <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "P2O5Ws","NWs", "SWs", 
             "KffactWs", "BFIWs", "PctHbWetXXXXWs", "PermWs",
             "PctCropXXXXWs", "RdDensWs", "N_inputs", "mean.n_tw.2018",
             "TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT")

varsPTL <-  c("WsAreaSqKm", "ElevWs", "RunoffWs", "P2O5Ws", "SWs", "KffactWs",
              "BFIWs", "PctHbWetXXXXWs", "PermWs","ClayWs", "PctCropXXXXWs", 
              "RdDensWs", "P_inputs","mean.pdep2013","TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT")

varsCL <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "NWs", "P2O5Ws", "SWs", "KffactWs", 
            "BFIWs", "PctHbWetXXXXWs", "PermWs", "ClayWs", "PctCropXXXXWs", "RdDensWs",
            "CoalMineDensWs", "MineDensWs", "mean.cl_tw.2018", "mean.s_tw.2018",
            "TMEAN_S_XXXX_PT","PSUMPY_XXXX_PT")

varsSO4 <- c("WsAreaSqKm", "ElevWs", "RunoffWs", "NWs", "P2O5Ws", "SWs", "KffactWs", 
             "BFIWs", "PctHbWetXXXXWs", "PermWs", "ClayWs", "PctCropXXXXWs", "RdDensWs",
             "CoalMineDensWs", "MineDensWs", "mean.s_tw.2018", "TMEAN_S_XXXX_PT", "PSUMPY_XXXX_PT")

#varsLSUB <- c("PSUMPY_XXXX_PT","TMEAN_S_XXXX_PT", "SandWs","PctAlluvCoastWs",
#             "NABD_NrmStorWs", "PctFstWs", "RdDensWs","PctCropXXXXWs",
#            "CoalMineDensWs", "MineDensWs", "W1_HAG", "W1_HNOAG", "W1H_CROP",
#           "KffactWs","WsAreaSqKm", "ElevWs", "RunoffWs","RckDepWs", "SLOPE","strPow")

varsLSUB <- c("SandWs", "PctAlluvCoastWs", "NABD_NrmStorWs",  "RdDensWs",
              "PctCropXXXXWsRp100", "PctFstWs",
              "CoalMineDensWs", "MineDensWs", "W1_HAG", 
              "KffactWs","WsAreaSqKm", "ElevWs", "RckDepWs", "strPow","RunoffWs")
#"PctFstWs",
#xdata$NABD_NrmStorWs2<-xdata$NABD_NrmStorWs/xdata$RunoffWs
#varsLDMB <- c("PSUMPY_XXXX_PT","TMEAN_S_XXXX_PT", "SandWs","PctAlluvCoastWs",
#             "NABD_NrmStorWs", "PctFstWs", "RdDensWs","PctCropXXXXWs",
#            "CoalMineDensWs", "MineDensWs", "W1_HAG", "W1_HNOAG", "W1H_CROP",
#            "KffactWs","WsAreaSqKm", "ElevWs", "SLOPE", "RunoffWs","RckDepWs")

#varsFNSA <- c("PSUMPY_XXXX_PT","TMEAN_S_XXXX_PT", "SandWs","PctAlluvCoastWs",
#             "NABD_NrmStorWs", "PctFstWs", "RdDensWs","PctCropXXXXWs",
#            "CoalMineDensWs", "MineDensWs", "W1_HAG", "W1_HNOAG", "W1H_CROP",
#           "KffactWs","WsAreaSqKm", "ElevWs", "SLOPE", "RunoffWs","RckDepWs")

#add PCT_FN to match NRSA thresholds
varsALL <- unique(c("UID", "UNIQUE_ID", "SITE_ID", "VISIT_NO", "YEAR", "COMID",
                    "LAT_DD", "LON_DD", "ECOWSA9","NTL", "PTL", "CL", "SO4", 
                    "LSUB_DMM","W1_HALL", varsLSUB, varsCL, varsSO4, varsPTL, 
                    varsNTL))
#############################################

xdata <- xdata[,varsALL]





