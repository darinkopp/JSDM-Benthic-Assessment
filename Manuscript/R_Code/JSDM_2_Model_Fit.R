library(Hmsc)
library(ape)
#run this once before starting the script
old.wd <- "C:/Users/DKOPP/OneDrive - Environmental Protection Agency (EPA)/HMSC/Manuscript_JSDM_OE_Assessment/Data_Analysis_Rscripts"

#run each regions on different R sessions
# c("XER","TPL","NPL","SAP","NAP","CPL","UMW","SPL","WMT")

# Region of Interest - regions were run in separate R sessions to use more 
# computational resources
Regions <- c("CPL","UMW")

timereport <- data.frame()
# Outter loop: iterate through thinning allows fitted models to be written without 
# waiting on models that need longer to converge
for (q in c("1A", "10A", "10B", "100A", "100B", "1000A", "1000B")){
  # q <- "1A" 
  print(q)
  
  #model parameters
  nChains = 3 #number of chains
  
  #take a total of 1000 samples by recording 5th observation 
  #for 1000 samples with thin = 5 is 5000 interations
  # model parameters
  thinning <- as.numeric(substr(q, start = 1, (nchar(q)-1)))
  thin = thinning #number of iterations per sample
  samples = 1000 #number of samples 
  transient  = 50 * thin #additional iterations used for warmup 
  verbose = 50 * thin #print progress
  
  setwd(old.wd)
  for (ROI in Regions){
    # ROI <- "XER" 
    roi <- ROI
    print(ROI)
    
    #point directory to files 
    setwd(paste0("Model_Fit_JSDM_Final/", ROI))
    
    candvar <- c("LAT_DD", "LON_DD", 
                 "NTL", "PTL", "CL", "SO4", 
                 "TMEAN_S_XXXX_PT", "PSUMPY_XXXX_PT",
                 "W1_HALL", "LSUB_DMM")
    
    #Taxa matrix
    Y <- read.csv(paste0("Model_Data/Taxa_", ROI, ".csv"), row.names = "X")
    
    #environmental Data 
    XData <- read.csv(paste0("Model_Data/Sites_", ROI, ".csv"), row.names = "X")
    XData <- XData[, candvar[-c(1,2)]]
    
    #phylogeny data 
    phyloTree <- read.csv(paste0("Model_Data/Taxonomy_",ROI,".csv"), row.names = "X", stringsAsFactors = TRUE)
    phyloTree <- as.phylo(~PHYLUM/CLASS/ORDER/FAMILY/GENUS, data = phyloTree, collapse = F)
    phyloTree$edge.length = rep(1, length(phyloTree$edge))
    
    #sample design
    sample.id = as.factor(rownames(Y))
    studyDesign = data.frame(sample = sample.id)
    rL = HmscRandomLevel(units = studyDesign$sample)
    rL$nfMax=2
    
    #formulas 
    XFormula = ~NTL + PTL + CL + SO4 + 
      TMEAN_S_XXXX_PT + PSUMPY_XXXX_PT + 
      W1_HALL + LSUB_DMM
    
    #specify model 
    model.env <- Hmsc(Y = 1 * (Y > 0), 
                      XData = XData, 
                      XFormula = XFormula,
                      phyloTree = phyloTree, 
                      studyDesign = studyDesign, 
                      ranLevels = list("sample" = rL), 
                      dist = "probit")
    
    #fit model
    T1 <- Sys.time()
    print(T1)
    
    model.env <- sampleMcmc(model.env, 
                            thin = thin, 
                            samples = samples, 
                            transient = transient,
                            nChains = nChains, 
                            verbose = verbose, 
                            initPar = "fixed effects", #sets initial condition for MCMC using ML (pg 94)
                            nParallel = nChains) #number of cores cannot be greater than the chains   
    T2 <- Sys.time()
    print(T2 - T1)
    
    save(model.env, file = paste0("Output/Fit_Models/HMSC_Model_", ROI,"_", q, "_1819"))
    
    
    #if A has already run, add it to B to increase sample size
    if(length(grep("B", q)) > 0){
      model.env_A <- model.env
      
      initRun <- paste0("Output/Fit_Models/HMSC_Model_", ROI, "_", paste0(thin, "A"),"_1819")
      load(initRun)
      model.env_B<-model.env
      
      model.env <- c(model.env_A, model.env_B)
    }
    
    # Check Convergence of model, if satisfactory break loop
    # potential scale reduction factor: goal for PSRF < 1.1. This is slightly high 
    # but acceptable given model complexity. 
    #####
    mpost <- convertToCodaObject(model.env, 
                                 spNamesNumbers = c(T, F), 
                                 covNamesNumbers = c(T, F)) 
    
    # Calculate convergence diagnostics (psrf)
    psrf.Beta <- gelman.diag(mpost$Beta, multivariate = F)$psrf
    
    # Phylo parms
    psrf.Rho <- gelman.diag(mpost$Rho, multivariate = F)$psrf
    
    # Assoc. parms: Random sub-sample to avoid excessive computation (Pg 305) 
    tmp <- mpost$Omega[[1]]
    z <- dim(tmp[[1]])[2]
    if(z > 1000){ sel <- sample(1:z, size = 1000)} else {sel <- 1:z}
    for(j in 1:length(tmp)){tmp[[j]] = tmp[[j]][,sel]}
    psrf.Omega <- gelman.diag(tmp, multivariate = F)$psrf
    
    
    # check diagnostics
    Beta.psrf.UCI = quantile(psrf.Beta[,"Upper C.I."], 0.99) < 1.1
    Omega.psrf.UCI = quantile(psrf.Omega[,"Upper C.I."], 0.99) < 1.1
    Rho.psrf.UCI = psrf.Rho[,"Upper C.I."] < 1.1
    
    print(c(beta = Beta.psrf.UCI,
            rho = Rho.psrf.UCI, 
            omega = Omega.psrf.UCI))
    
    timereport <- rbind(timereport,data.frame(ROI, q, time = T2 - T1,
                                              beta = Beta.psrf.UCI,                      
                                              rho = Rho.psrf.UCI, 
                                              omega = Omega.psrf.UCI))
    
    # break loop if successfully converged
    # update region vector so it is not re-run with more 
    # iterations
    if(all(Beta.psrf.UCI, 
           Rho.psrf.UCI, 
           Omega.psrf.UCI)){
      
      Regions <- Regions[Regions != ROI]
      #break
    }
    
    # reset directory
    setwd(old.wd)
  }
}
