# Validation of JSDM using Visit 2
# analysis and explanation provided in supplement. 

library(Hmsc)
library(ggplot2)
###############################################
FittedModelPath <- "Model_Fit_JSDM_Final"
modelnames <- grep("1819", 
                   list.files(FittedModelPath, 
                              full.names = F, 
                              recursive = F), 
                   value = T)
# thresholds used for model fitting to remove taxa with 
# rare predicted occurrence probabilities
P.threshold <- read.csv("Figures/Table3_jsdmPerformance.csv", row.names = "Row.names")

par(mfrow = c(3,3))
post.check <- data.frame()
df.plot <- data.frame()
for (i in 1:9){
  # i <- 6
  mn <- modelnames[i]
  load(paste0(FittedModelPath, "/", mn)) # model.env
  roi <- unlist(lapply(strsplit(mn, "_"), "[", 3))
  
  sites.val <- read.csv(paste0("Model_Fit_JSDM_Final/", roi, "/Model_Data/ValSites_", roi, ".csv"))
  sites <- read.csv(paste0("Model_Fit_JSDM_Final/", roi, "/Model_Data/Sites_", roi, ".csv"))
  
  #Identifies visit 2 sites 
  newIDS <- sites[sites$UNIQUE_ID %in% sites.val$UNIQUE_ID, c("UNIQUE_ID","UID")]
  NewData <- sites.val
  NewData <-  merge(NewData, newIDS, by = "UNIQUE_ID")
  rownames(NewData) <- NewData$UID.y
  print(mn)
  print(nrow(NewData))

  #identifies visit 2 taxa
  taxon.val <- read.csv(paste0("Model_Fit_JSDM_Final/", roi, "/Model_Data/ValTaxa_", roi, ".csv"), 
                        row.names = "X")
  obs <- taxon.val[as.character(NewData$UID.x),]
  obs <- obs[,colnames(obs)[colSums(obs)>0]]
  
  # add columns for taxa that were not collected within the region during visit 2
  absentTaxa <- colnames(model.env$Y)[!colnames(model.env$Y) %in% colnames(obs)]
  if(length(absentTaxa)>0){
    absentTaxa <- setNames(data.frame(matrix(0,nrow=nrow(obs), ncol=length(absentTaxa))),
                           absentTaxa)
    obs <- cbind(obs, absentTaxa)
  }
  
  # select species from the regional species pool that were modeled
  obs <- obs[,colnames(model.env$Y)]
  yc <- obs
  
  # set all rownames for indexing/site matching
  yc <- yc[as.character(NewData$UID.x),]
  rownames(yc) <- NewData$UID.y
  rownames(obs) <- NewData$UID.y
  
  # between site visits, there are a number of pheonomena that could influence 
  # the presence or absence of a taxon beyond changes in environmental conditions, 
  # including emergence, colonization (e.g. mass effects) or differences in 
  # labratory sub-sampling. Because the model does not account for these differences, 
  # it is unreasonable to explect it to preform well at a given time period, especially
  # if environmental conditons have not changed drastically. Thus we evaluate performance 
  # based on taxa that were either present or absent in both sampling
  
  # if a taxon was collected during visit 1 but NOT visit 2 treat it as known absent 
  a <- yc[rownames(yc), colnames(model.env$Y)] < 
    model.env$Y[rownames(yc), colnames(model.env$Y)]
  a <- ifelse(a, 0, 1)
  
  # if a taxon was collected during visit 2 but NOT visit 1 treat it as known present
  p <- yc[rownames(yc), colnames(model.env$Y)] > 
    model.env$Y[rownames(yc), colnames(model.env$Y)]
  p <- ifelse(p, 1, 0)
  
  # we are interested in predicting the occurrence probabilities for taxa that 
  # were either present or absent during both visits
  taxa2pred <- yc[rownames(yc), colnames(model.env$Y)]==  
    model.env$Y[rownames(yc), colnames(model.env$Y)] 
  taxa2pred <- ifelse(taxa2pred, NA, 9999)
  
  d <- taxa2pred * a * p
  d[d == 9999] <- 1
  
  # r<-29
  # check<-data.frame(a[rownames(yc)[r],],
  #                   p[rownames(yc)[r],],
  #                   taxa2pred[rownames(yc)[r],],
  #                   d[rownames(yc)[r],],
  #                   t(yc[rownames(yc)[r],]),
  #                   model.env$Y[rownames(yc)[r],])
  # head(check, 20)
  yc <- d
  
  
  #rownames(NewData)%in%rownames(model.env$Y)
  NewData <- NewData[,colnames(model.env$XData)]
  sample.id = as.factor(rownames(NewData))
  studyDesign = data.frame(sample = sample.id)
  rL = HmscRandomLevel(units = studyDesign$sample)
  
  PVHC <- predict(model.env,
                  XData = NewData,
                  studyDesign = studyDesign, 
                  rL = rL, expected = T, Yc=yc)
  
  PVHC = abind::abind(PVHC, along = 3)
  
  # set absent taxa to zero so that they are not included in the prediction 
  #PVHC[,colnames(absentTaxa),] <- 0
  
  # so each site has a unique combination of NAs
  focal.taxa.richness.pred <- data.frame()
  for (r in rownames(yc)){
    #r <- rownames(yc)[1]
    #dim(model.env$Y)[]
    #sum(!is.na(yc[r,]) & yc[r,] == 1)
    #dim(PVHC[r, is.na(yc[r,]), 1])
    #x <- PVHC[r, is.na(yc[r,]), 1]
    #predFreq <- (rowMeans(colMeans(PVHC[, is.na(yc[r,]), ])))
    #predFreq <- colMeans(obs[,is.na(yc[r,])])
    #predFreq <- colMeans(model.env$Y[,is.na(yc[r,])])
    #meanPred <- apply(PVHC[, is.na(yc[r,]), ], c(1,2), mean))
    #s<-meanPred*model.env$Y[rownames(meanPred),is.na(yc[r,])]
    
    #sum(x[x>=predFreq[names(x)]])
    pth <- P.threshold[roi,"th"]
    #pth <- 1-mean(is.na(yc[r,]))
    #focal.taxa.richness.tmp <- apply(PVHC[r, , ], 2, function(x){sum(x)})
    focal.taxa.richness.tmp <- apply(PVHC[r, , ], 2, function(x){sum(x[x>pth])})
    
    #focal.taxa.richness.tmp <- apply(PVHC[r, is.na(yc[r,]), ], 2, function(x){sum(x)})
    #focal.taxa.richness.tmp <- apply(PVHC[r, is.na(yc[r,]), ], 2, function(x){sum(x[x>pth])})
    #focal.taxa.richness.tmp <- apply(PVHC[r, is.na(yc[r,]), ], 2, function(x){sum(sort(x, decreasing = T)[1:floor(sum(x))])})
    #focal.taxa.richness.tmp <- apply(PVHC[r, is.na(yc[r,]), ], 2, function(x){sum(x[x>=predFreq[names(x)]])})
    #focal.taxa.richness.tmp <- apply(PVHC[r, , ], 2, function(x){mean(x)*sum(is.na(yc[r,]))})
    #focal.taxa.richness.tmp <- apply(PVHC[r, is.na(yc[r,]), ], 2, function(x){sum(x[!names(x)%in%colnames(absentTaxa)])})
    
    focal.taxa.richness.tmp <- data.frame(t(focal.taxa.richness.tmp), row.names = r)
    focal.taxa.richness.pred <- rbind(focal.taxa.richness.pred, focal.taxa.richness.tmp)
  }
  
  focal.taxa.richness.obs <- data.frame()
  for (r in rownames(yc)){
    #r <- rownames(yc)[1]
    #dim(PVHC[r, is.na(yc[r,]), ])
    focal.taxa.richness.tmp <- data.frame(obs.sum = sum(obs[r, ]),row.names = r, P = sum(yc[r,],na.rm=T))
    #focal.taxa.richness.tmp <- data.frame(obs.sum = sum(obs[r, is.na(yc[r,])]), row.names = r)
    #focal.taxa.richness.tmp <- data.frame(obs.sum = sum(obs[r, ])-sum(yc[r,],na.rm=T),row.names = r)
    focal.taxa.richness.obs <- rbind(focal.taxa.richness.obs, focal.taxa.richness.tmp)
  }
  
  
  pred.mean <- apply(focal.taxa.richness.pred,1,mean)
  pred.lwr <- apply(focal.taxa.richness.pred,1,quantile, probs = 0.0001)
  pred.upr <- apply(focal.taxa.richness.pred,1,quantile, probs = 0.9999)
  df.tmp <- data.frame(roi, pred.mean, pred.lwr, pred.upr, 
                       obs.sum = focal.taxa.richness.obs[names(pred.mean),]) 
  
  df.plot <- rbind(df.plot, df.tmp)
  
  print(summary(lm(obs.sum.obs.sum~pred.mean, data = df.tmp)))
  
  outside <- NULL
  for (j in 1:dim(focal.taxa.richness.obs)[1]){
    #j<-1
    if(focal.taxa.richness.obs[j,"obs.sum"] < 
       min(focal.taxa.richness.pred[j,]) | 
       focal.taxa.richness.obs[j,"obs.sum"] > 
       max(focal.taxa.richness.pred[j,])){
      outside <- append(outside, 1,  length(outside))
    } else {
      outside <- append(outside, 0,  length(outside))
    }
  }
  
  post.check <- rbind(post.check, 
                      data.frame(roi, n=dim(focal.taxa.richness.obs)[1], 
                                 withinPost = round(1-mean(outside),2)))
}

# add color to points outside posterior predicted distribution
df.plot$outside<-
  ifelse(df.plot$obs.sum.obs.sum < df.plot$pred.lwr | 
           df.plot$obs.sum.obs.sum>df.plot$pred.upr,
         "Outside","")

Plot <- ggplot(df.plot, aes(x = pred.mean, 
                          y = obs.sum.obs.sum, 
                          color = outside))+
  geom_segment(aes(x=pred.mean, y = pred.upr, 
                   xend = pred.mean, yend = pred.lwr),
               color = "darkgrey")+
  geom_point(size = 2,
             aes(color = factor(outside)),
             shape = 21)+
  geom_abline(slope = 1, intercept = 0,lty = 2, lwd = 1.05)+
  #geom_smooth(method = lm, formula = y ~ x, se = FALSE)+
  facet_wrap("roi")+ 
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black","red"))+ 
  labs(x="posterior prediction of genus richness", y="observed genus richness")

#windows()
#Plot

ggsave(filename = "Figures/FigureS1_Visit2Validation_Plot.jpeg",
       Plot,width=5,height=5)
