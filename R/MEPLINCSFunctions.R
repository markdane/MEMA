#' Create pseudeoimages of 8 well MEMA data
#' @param DT A datatable with ArrayRow, ArrayColumn and pr columns
#' @param pr The name of the column of values to be shown in the images
#' @param prDisplay A name for the legend title
#' @return The ggplot grob to be displayed
#' @export
create8WellPseudoImage <- function(DT, pr, prDisplay){
  highThresh = .998
  #move outliers to maximum displayed value
  DT[[pr]][DT[[pr]]>=quantile(DT[[pr]],probs = highThresh,na.rm=TRUE)] <- as.integer(quantile(DT[[pr]],probs = highThresh,na.rm=TRUE))
  p <- ggplot(DT, aes_string(x="ArrayColumn", y="ArrayRow",colour=pr))+
    geom_point(size=rel(.3))+
    scale_y_reverse()+   scale_x_continuous(breaks= c(min(DT$ArrayColumn),round(mean(c(min(DT$ArrayColumn),max(DT$ArrayColumn)))),max(DT$ArrayColumn)))+
    scale_colour_gradient(low = "white", high = "red")+
    guides(colour = guide_legend(prDisplay, keywidth = .5, keyheight = .5))+
    ggtitle(paste("\n",prDisplay,"for plate",unique(DT$Barcode)))+
    xlab("")+ylab("")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.8)),
          axis.title.x = element_text(size=rel(.5)),
          plot.title = element_text(size = rel(.8)),
          strip.text = element_text(size = rel(.5)),
          legend.text=element_text(size = rel(.4)),legend.title=element_text(size = rel(.3)))+
    facet_wrap(~Well_Ligand, ncol=4)
}

#' Create histograms of 8 well MEMA data
#' @param DT A datatable with a pr column of values
#' @param pr The name of the column of values to be shown in the images
#' @param prDisplay A name for the legend title
#' @param binwidth Optional binwidth value for the histograms
#' @param upperProb Optional upper probablity to limit the display
#' @param ncol Optional number of columns to display in the ggplot facet
#' @return The ggplot grob to be displayed
create8WellHistograms <- function(DT, pr, prDisplay, binwidth = diff(quantile(DT[[pr]], probs = c(0,.98),na.rm=TRUE))/50, upperProb = .99, ncol = 4) {
  
  p <- ggplot(DT, aes_string(x=pr))+
    geom_histogram(binwidth = binwidth)+
    scale_x_continuous(limits = quantile(DT[[pr]],probs = c(0,upperProb),na.rm=TRUE))+
    ggtitle(paste("\n\n",prDisplay,"in",unique(DT$CellLine), "cells in plate",unique(DT$Barcode)))+
    xlab(prDisplay)+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.8)),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), 
          plot.title = element_text(size = rel(.5)),
          strip.text = element_text(size = rel(.5)),
          legend.text=element_text(size = rel(.3)),
          legend.title=element_text(size = rel(.3)))+
    facet_wrap(~Well, ncol=ncol)
}

# #' Display a heatmap from a datatable
# #' @param DT A datatable with MEP, Barcode and value columns
# #' @param title Optional title for the heatmap
# #' @param cols
# #' @export
# ubHeatmap <- function(DT, title = NULL, cols) {
#   #Remove MEP and Barcode columns and convert to a matrix
#   m <- as.matrix(DT[,grep("MEP|Barcode",colnames(DT),invert=TRUE), with=FALSE])
#   
#   #This assignment of names retains the order after matrix coercion
#   rownames(m) <- DT$MEP
#   
#   #Cluster the active MEPs, scaling the inputs
#   heatmap.2(m,scale="column", col = bluered, trace = "none", cexRow=.5, cexCol=.9, key=FALSE, main = paste(title), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1,5.0), lwid=c(.5,0.2,2.5,1.5),
#             mar=c(25,1),
#             RowSideColors=cols[DT$Barcode], colRow = cols[DT$Barcode], na.rm = TRUE)
#   return(m)
# }
# 
# #' @export
# #' 
# heatmapNoBC <- function(DT, title = NULL, cols = plateCol, activeThresh = .95) {
#   
#   activeFV <- DT
#   #Remove MEP and Barcode columns and convert to a matrix
#   activeFVM <- as.matrix(activeFV[,grep("MEP|Barcode",colnames(activeFV),invert=TRUE), with=FALSE])
#   
#   #This assignment of names retains the order after matrix coercion
#   rownames(activeFVM) <- names(DT$MEP)
#   
#   #Cluster the active MEPs, scaling the inputs
#   #plot.new()
#   heatmap.2(activeFVM,scale="column", col = bluered, trace = "none", cexRow=.5, cexCol=.9, key=FALSE, main = paste(title), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
#             lwid=c(1.5,0.2,2.5,2.5),mar=c(20,5), RowSideColors=cols[activeFV$Barcode], colRow = cols[activeFV$Barcode])
#   
#   return(activeFV)
# }


# #' @export
# #' 
# plotTotalDAPI <- function(l1, barcodes){
#   for (barcode in barcodes){
#     mDT <- l1[l1$Barcode == barcode]
#     mDT <- mDT[mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi > quantile(mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi, probs=.01, na.rm=TRUE) & mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi < quantile(mDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi,probs=.98, na.rm=TRUE)]
#     mDT <- mDT[,DNAThresh := min(Nuclei_CP_Intensity_IntegratedIntensity_Dapi[Nuclei_PA_Cycle_State==2]), by="Well"]
#     p <- ggplot(mDT, aes(x=Nuclei_CP_Intensity_IntegratedIntensity_Dapi))+geom_histogram(binwidth = 50000)+
#       geom_vline(data = mDT, aes(xintercept = DNAThresh), colour = "blue")+
#       facet_wrap(~Well_Ligand, nrow=2, scales="free_x")+
#       ggtitle(paste("\n\n","Total DAPI Signal,",barcode))+
#       ylab("Count")+xlab("Total Intensity DAPI")+
#       theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0.5, size=rel(.5)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.3)),strip.text.x = element_text(size = rel(.5)))
#     suppressWarnings(print(p))
#   }
# }
# 
# #' @export
# #' 
# plotTotalDAPIINCell <- function(l1, barcodes){
#   for (barcode in barcodes){
#     mDT <- l1[l1$Barcode == barcode]
#     mDT <- mDT[mDT$Nuclei_IxANuc > quantile(mDT$Nuclei_IxANuc, probs=.01, na.rm=TRUE) & mDT$Nuclei_IxANuc < quantile(mDT$Nuclei_IxANuc,probs=.98, na.rm=TRUE)]
#     mDT <- mDT[,DNAThresh := min(Nuclei_IxANuc[Nuclei_PA_Cycle_State==2]), by="Well"]
#     p <- ggplot(mDT, aes(x=Nuclei_IxANuc))+geom_histogram(binwidth = 1e+04)+
#       geom_vline(data = mDT, aes(xintercept = DNAThresh), colour = "blue")+
#       facet_wrap(~Well_Ligand, nrow=2, scales="free_x")+
#       #xlim(0,quantile(mDT$TotalIntensityDAPI,probs=.98, na.rm=TRUE))+
#       ggtitle(paste("\n\n","Total DAPI Signal,",barcode))+
#       ylab("Count")+xlab("Total Intensity DAPI")+
#       theme(strip.text = element_text(size = 5))
#     suppressWarnings(print(p))
#   }
# }

#' Print pseudoimages and histograms of the spot cell count values
#' @param l3 A datatable with Barcode, ECMp, Spot_PA_SpotCellCount, Well_Ligand and QA_score columns
#' @param barcodes The barcodes of the data to be displayed
#' @param lthresh The value used to gate the poor quuality spots.
#' @return none This function is called to print the images
#' @export
plotSCCHeatmapsQAHistograms <- function(l3, barcodes, lthresh){
  for (barcode in barcodes){
    DT <-l3[l3$Barcode==barcode,]
    #Remove the fiducial entries
    setkey(DT,ECMp)
    DT <- DT[!"fiducial"]
    DT <- DT[!"blank"]
    
    p <- create8WellPseudoImage(DT, pr = "Spot_PA_SpotCellCount",prDisplay = "Spot Cell Count")
    suppressWarnings(print(p))
    wellScores <- unique(DT[,list(Well_Ligand, QAScore=sprintf("%.2f",QAScore))])
    
    p <- ggplot(DT, aes(x=Spot_PA_LoessSCC))+
      geom_histogram(binwidth=.04)+
      geom_vline(xintercept=lthresh, colour="blue")+
      geom_text(data=wellScores, aes(label=paste0("QA\n",QAScore)), x = 2, y = 30, size = rel(3), colour="red")+
      ggtitle(paste("\n QA on Loess Model of Spot Cell Count\n for plate",unique(DT$Barcode)))+xlab("Loess Normalized Spot Cell Count")+xlim(0,3)+
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.5)),
            axis.title.x = element_text( size=rel(.5)),
            axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)),
            axis.title.y = element_text( size=rel(.5)),
            plot.title = element_text(size = rel(.8)),
            strip.text = element_text(size = rel(.5)),
            legend.text=element_text(size = rel(.3)),
            legend.title=element_text(size = rel(.3)))+
      facet_wrap(~Well_Ligand, ncol=4)
    suppressWarnings(print(p))
  }
}

# #' @export
# #' 
# filterl4 <- function(dt,lowQALigands){
#   #Remove failed QA wells
#   l4QA<- dt[!dt$Ligand %in% lowQALigands]
#   
#   setkey(l4QA, "ECMp")
#   l4QA <- l4QA[!"blank"]
#   l4QA <- l4QA[!"fiducial"]
#   l4QA <- l4QA[,grep("Center_X|Center_Y|Center_Theta",colnames(l4QA),value = TRUE, invert = TRUE), with = FALSE]
#   
#   #Define features for clustering
#   fv <- paste("^Barcode","MEP",
#               "Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker_Norm",
#               "Nuclei_CP_AreaShape_Area_Norm",
#               "Nuclei_CP_AreaShape_Eccentricity_Norm",
#               "Nuclei_CP_AreaShape_Perimeter_Norm",
#               "Nuclei_CP_Intensity_MedianIntensity_Dapi_Norm",
#               "Spot_PA_SpotCellCount_Norm",
#               "Nuclei_PA_AreaShape_Neighbors_Norm",
#               "Nuclei_PA_Cycle_DNA2NProportion_Norm$",
#               "Nuclei_CP_Intensity_MedianIntensity_Edu_Norm",
#               "Nuclei_PA_Gated_EduPositiveProportion_Norm_Norm",
#               "Cytoplasm_CP_Intensity_MedianIntensity_KRT19_Norm",
#               "Cytoplasm_CP_Intensity_MedianIntensity_KRT5_Norm",
#               "Cytoplasm_PA_Intensity_LineageRatio_Norm$",
#               sep="$|^")
#   
#   fv <- grep(fv, colnames(l4QA), value = TRUE)
#   #Create numeric feature vectors datatable
#   fvDT <- l4QA[,fv,with = FALSE]
#   return(fvDT)
# }
# 
# #' @export
# #' 
# filterl4RUV <- function(dt,lowQALigands){
#   #Remove failed QA wells
#   l4QA<- dt[!dt$Ligand %in% lowQALigands]
#   
#   setkey(l4QA, "ECMp")
#   l4QA <- l4QA[!"blank"]
#   l4QA <- l4QA[!"fiducial"]
#   l4QA <- l4QA[,grep("Center_X|Center_Y|Center_Theta",colnames(l4QA),value = TRUE, invert = TRUE), with = FALSE]
#   
#   #Define features for clustering
#   fv <- paste("^Barcode","MEP",
#               "Cytoplasm_CP_Intensity_MedianIntensity_MitoTrackerLog2RUVLoess",
#               "Nuclei_CP_AreaShape_AreaLog2RUVLoess",
#               "Nuclei_CP_AreaShape_EccentricityLog2RUVLoess",
#               "Nuclei_CP_AreaShape_PerimeterLog2RUVLoess",
#               "Nuclei_CP_Intensity_MedianIntensity_DapiLog2RUVLoess",
#               "Spot_PA_SpotCellCountLog2RUVLoess",
#               "Nuclei_PA_AreaShape_NeighborsLog2RUVLoess",
#               "Nuclei_PA_Cycle_DNA2NProportionLogitRUVLoess$",
#               "Nuclei_CP_Intensity_MedianIntensity_EdULog2RUVLoess",
#               "Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess",
#               "Cytoplasm_CP_Intensity_MedianIntensity_KRT19Log2RUVLoess",
#               "Cytoplasm_CP_Intensity_MedianIntensity_KRT5Log2RUVLoess",
#               "Cytoplasm_PA_Intensity_LineageRatioLog2RUVLoess$",
#               sep="$|^")
#   
#   fv <- grep(fv, colnames(l4QA), value = TRUE)
#   #Create numeric feature vectors datatable
#   fvDT <- l4QA[,fv,with = FALSE]
#   return(fvDT)
# }
# 
# #' @export
# #' 
# plotSCCRobustZScores <- function(dt, thresh = 3){
#   #Filter our FBS MEPs then plot spot cell count robust Z scores
#   #browser()
#   dt <- dt[!grepl("FBS",dt$MEP)]
#   p <- ggplot(dt, aes(x=Spot_PA_SpotCellCountLog2RUVLoess_RobustZ))+geom_histogram(binwidth = .1)+
#     geom_vline(xintercept = c(-thresh,thresh), colour = "blue")+
#     ggtitle(paste("\n\n","MEP Normalized Spot Cell Count Robust Z Scores Distribution"))+
#     ylab("Count")+xlab("Normalized Spot Cell Count Robust Z Scores")+
#     theme(strip.text = element_text(size = 5))
#   suppressWarnings(print(p))
#   
# }
# 
# 
# #' @export
# #' 
# combineSSs <- function(SSs){
#   #browser()
#   l4List <- lapply(SSs, function(ss){
#     l4 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine,"_",ss,"_Level4.txt"), showProgress = FALSE)
#     setkey(l4,"ECMp")
#     l4 <- l4[!"fiducial"]
#     l4 <- l4[!"blank"]
#     l4$SS <- ss
#     return(l4)
#   })
#   
#   l4SS1 <- l4List[[1]]
#   l4SS2 <- l4List[[2]]
#   l4SS3 <- l4List[[3]]
#   
#   setkey(l4SS1,"LigandAnnotID","ECMpAnnotID")
#   setkey(l4SS2,"LigandAnnotID","ECMpAnnotID")
#   setkey(l4SS3,"LigandAnnotID","ECMpAnnotID")
#   
#   #Bind the data
#   DT <- data.table(l4SS1, l4SS2, l4SS3, check.names = TRUE)
# }
# 
# 
# #' @export
# #' 
# integrateSSs <- function(SSs, cellLine = "PC3"){
#   #browser()
#   l4List <- lapply(SSs, function(ss){
#     l4 <- fread(paste0("./",cellLine,"/",ss,"/AnnotatedData/",cellLine,"_",ss,"_Level4.txt"), showProgress = FALSE)
#     setkey(l4,"ECMp")
#     l4 <- l4[!"fiducial"]
#     l4 <- l4[!"blank"]
#     setkey(l4, "MEP")
#     l4$SS <- ss
#     return(l4)
#   })
#   
#   l4SS1 <- l4List[[1]]
#   l4SS2 <- l4List[[2]]
#   l4SS3 <- l4List[[3]]
#   
#   setkey(l4SS1,"LigandAnnotID","ECMpAnnotID")
#   setkey(l4SS2,"LigandAnnotID","ECMpAnnotID")
#   setkey(l4SS3,"LigandAnnotID","ECMpAnnotID")
#   
#   #Bind the data using the common MEPs
#   DT <- data.table(l4SS1, l4SS2, l4SS3, check.names = TRUE)
#   
#   #Median summarize the FBS rows
#   setkey(DT,"MEP")
#   DTFBS <- DT[grepl("FBS", DT$MEP)]
#   #Get the medians of each numeric parameter
#   parms <- colnames(DTFBS)[unlist(lapply(DTFBS,class)) %in% c("numeric","integer")]
#   FBSMedians <- data.frame(t(as.matrix(apply(DTFBS[, parms,with=FALSE],2,median))),MEP="FBS", stringsAsFactors = FALSE)
#   
#   #Merge the metadata back in with the data
#   metadata <- colnames(DTFBS)[!unlist(lapply(DTFBS,class)) %in% c("numeric","integer")]
#   FBSMetadata <- DTFBS[, metadata, with = FALSE]
#   FBSMetadata$MEP <- "FBS"
#   FBSMetadata$MEP.1 <- "FBS"
#   FBSMetadata$MEP.2 <- "FBS"
#   FBSMetadata$ECMp <- NA
#   FBSMetadata$ECMp.1 <- NA
#   FBSMetadata$ECMp.2 <- NA
#   FBSMetadata$ECMpAnnotID <- NA
#   FBSMetadata$ECMpAnnotID.1 <- NA
#   FBSMetadata$ECMpAnnotID.2 <- NA
#   FBSMetadata$Well <- NA
#   FBSMetadata$Well.1 <- NA
#   FBSMetadata$Well.2 <- NA
#   FBSMetadata$Barcode <- NA
#   FBSMetadata$Barcode.1 <- NA
#   FBSMetadata$Barcode.2 <- NA
#   
#   FBSMetadata <- unique(FBSMetadata)
#   
#   FBSMisOrdered <- cbind(FBSMetadata[,MEP:=NULL],FBSMedians)
#   
#   #Replace all FBS rows with one row of medians as the last row
#   DT1FBS<- rbind(DT[!grepl("FBS", DT$MEP)],FBSMisOrdered,use.names=TRUE)
#   
# }
# 
# #' @export
# #' 
# createMEMATestRawData <- function(cellLine, ss, nrRows, nrCols){
#   library(XLConnect)
#   #Create a subset of raw data from a complete staining set
#   #Expects to find raw data in the ./cellLine/ss/Rawdata folder
#   #Will create or overwrite a folder named cellLine_TestData
#   
#   imageNumbers <- unlist(lapply(1:nrRows, function(x, nrCols){
#     xseq <- (x-1)*20+(1:nrCols)
#   }, nrCols = nrCols))
#   
#   
#   cellLineFiles <- dir(paste(".",cellLine,ss,"RawData","v1", sep = "/"),
#                        pattern = ".csv",full.names = TRUE)
#   for( x in cellLineFiles){
#     dt <-fread(x)
#     dtTest <- dt[dt$CP_ImageNumber %in% imageNumbers]
#     write.table(format(dtTest, digits = 4, trim=TRUE), file = sub("LI8X","LI8T",sub(cellLine,paste0(cellLine,"TestData"),x)),
#                 sep = ",", row.names = FALSE, quote=FALSE)
#   }
#   
#   metadataFiles <- dir(paste(".",cellLine,ss,"Metadata", sep = "/"),
#                        pattern = ".xlsx",full.names = TRUE)
#   for( x in metadataFiles){
#     wb<-loadWorkbook(x)
#     saveWorkbook(wb, file = sub("LI8X","LI8T",sub(cellLine,paste0(cellLine,"TestData"),x)))
#   }
# }

# #' @export
# #' 
# heatmapFromFBS <- function(DT, title = NULL, cols = plateCol, activeThresh = .95) {
#   #browser()
#   DT$Barcode <- as.factor(DT$Barcode)
#   #Get medians of high serum numeric features
#   fvDTHS <- DT[grepl("FBS", DT$MEP)]
#   hsMedians <- data.frame(t(as.matrix(apply(fvDTHS[,grep("MEP|Barcode", colnames(fvDTHS),invert=TRUE),with=FALSE],2,median, na.rm = TRUE))),MEP="FBS", Barcode = NA, stringsAsFactors = FALSE)
#   #Replace all FBS rows with one row of medians as the last row
#   DT<- rbind(DT[!grepl("FBS", DT$MEP)],hsMedians)
#   
#   #Calculate the dist matrix with euclidean method
#   dmm <- as.matrix(dist(DT[,grep("MEP|Barcode",colnames(DT),invert=TRUE), with=FALSE]), labels=TRUE)
#   #Extract the distance to the high serum medians
#   distHS <- dmm[which(DT$MEP == "FBS"),]
#   #Name the distance values
#   names(distHS) <- DT$MEP
#   
#   #Select the most active by distance from the median control fv
#   dmmThresh <- quantile(distHS,probs = activeThresh, na.rm = TRUE)
#   activeMEPs <- distHS[distHS>dmmThresh]
#   #browser()
#   #Create an active MEP subset matrix of the normalized data
#   activeFV <- DT[DT$MEP %in% names(activeMEPs)]
#   #Remove MEP and Barcode columns and convert to a matrix
#   activeFVM <- as.matrix(activeFV[,grep("MEP|Barcode",colnames(activeFV),invert=TRUE), with=FALSE])
#   
#   #This assignment of names retains the order after matrix coercion
#   rownames(activeFVM) <- activeFV$MEP
#   
#   #Cluster the active MEPs, scaling the inputs
#   heatmap.2(activeFVM,scale="column", col = bluered, trace = "none", cexRow=.5, cexCol=.9, key=FALSE, main = paste(title), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1,5.0), lwid=c(.5,0.2,2.5,1.5),
#             mar=c(25,1),
#             RowSideColors=cols[activeFV$Barcode], colRow = cols[activeFV$Barcode], na.rm = TRUE)
#   return(activeFV)
# }

#' @export
#' 
dfZoom <- function(x, min=.02, max=1){
  minMax <- quantile(unlist(x), probs=c(min,max), na.rm=TRUE)
  cl <- t(apply(x,1, function(c){
    c[c<minMax[1]] <- minMax[1]
    c[c>minMax[2]] <- minMax[2]
    return(c)
  }))
  return(data.frame(cl))
}
# 
# #' @export
# #' 
# fvZoom <- function(x, min=.02, max=1){
#   minMax <- quantile(unlist(x), probs=c(min,max), na.rm=TRUE)
#   x[x<minMax[1]] <- minMax[1]
#   x[x>minMax[2]] <- minMax[2]
#   return(x)
# }
# 
# #' @export
# #' 
# createl4KeepRaw <- function (l3) 
# {
#   #Select all signal names and the minimal metadata
#   l4Names <- grep("^Ligand$|^ECMp$|Barcode|MEP|^Well$|StainingSet|_CP_|_PA_",x = names(l3), value = TRUE)
#   #Remove the SE names
#   l4Names <- grep("_SE", l4Names, value = TRUE, invert = TRUE)
#   #Create a datatable of the signals and minimal metadata
#   l4Keep <- l3[, l4Names, with = FALSE]
#   l4DT <- l4Keep[, lapply(.SD, numericMedian), keyby = "Ligand,ECMp,MEP,Barcode,Well,StainingSet"]
#   return(l4DT)
# }
# 
# #' Add an ECMp and Ligand RLE value for each signal in a dataset
# #'
# #' 
# addRLEs <- function(dt){
#   #Reduce size of dt for debug
#   #dt <- dt[,grepl("StainingSet|^Well$|^MEP$|^ECMp$|^Ligand$|AreaShape_Area",colnames(dt)), with=FALSE]
#   if(!any(grepl("MEP",colnames(dt)))) stop ("There must be a MEP column to when adding RLE values")
#   if(!any(grepl("ECMp",colnames(dt)))) stop ("There must be an ECMp column to when adding RLE values")
#   if(!any(grepl("Ligand",colnames(dt)))) stop ("There must be a Ligand column to when adding RLE values")
#   dtRLEList <- apply(dt[,grepl("_CP_|_PA_|_IC_",colnames(dt)), with=FALSE], 2, function(signalValues){
#     DT <- data.table(Signal=signalValues,
#                      StainingSet=dt$StainingSet,
#                      ECMp=dt$ECMp,
#                      Ligand=dt$Ligand)
#     DT <- DT[,SignalLigandRLE := calcResidual(Signal), by="StainingSet,ECMp"]
#     DT <- DT[,SignalECMpRLE := calcResidual(Signal), by="StainingSet,Ligand"]
#   })
#   
#   #Copy signal name into column names
#   dtRLEList <- lapply(names(dtRLEList),function(dtName){
#     setnames(dtRLEList[[dtName]], "Signal", dtName)
#     setnames(dtRLEList[[dtName]], "SignalLigandRLE", paste0(dtName,"LigandRLE"))
#     setnames(dtRLEList[[dtName]], "SignalECMpRLE", paste0(dtName,"ECMpRLE"))
#     setkey(dtRLEList[[dtName]],StainingSet,ECMp,Ligand)
#     return(dtRLEList[[dtName]])
#   })
#   
#   #Combine list elements back into a common datatable
#   dtRLE <- Reduce(merge, dtRLEList)
#   
#   #Add the non-signal metadata columns from the original dt
#   mdColumns <- dt[,!grepl("_CP_|_PA_|_IC_",colnames(dt)), with=FALSE]
#   dtmdRLE <- merge(mdColumns,dtRLE,by=c("StainingSet","ECMp","Ligand"))
#   
#   return(dtmdRLE)
# }
# 
# 

# #'
# #' 
# numericMedianUniqueMetadata<-function(x){
#   if(is.numeric(x)){
#     as.numeric(median(x))
#   } else {
#     if(!length(unique(x))==1){
#       stop("metadata of summarized values must be identical")
#     } else{
#       unique(x)
#     }
#   }
# }
# 
# 

# #' @export
# #' 
# getMEMADataFileNames <- function(dataset){
#   library(stringr)
#   library(magrittr)
#   
#   #return a dataframe with columns: Barcode, Path, Type
#   platform <- regmatches(version$platform, regexpr("apple|linux", version$platform))
#   if(platform=="linux") {
#     path<- dataset[["Path"]]
#   } else if(platform=="apple"){
#     path<- dataset[["LocalPath"]]
#   } else stop("Code is running on unknown platform",platform)
#   
#   rdfnL <- lapply(dataset[grepl("Barcode",names(dataset))], function(barcode, version, path){
#     if(!length(dir(paste0(path,barcode,"/Analysis/",version),pattern = "csv"))==0){
#       Path=dir(paste0(path,barcode,"/Analysis/",version),pattern = "csv", full.names = TRUE)
#       Well <- str_extract(Path,"_[[:alpha:]][[:digit:]]{2}_") %>%
#         str_replace_all("_","")
#       Location <- str_extract(Path,"Image|Cells|Nuclei|Cytoplasm")
#       rdfns <- data.table(Barcode=barcode, Path=Path, Type="Raw", Well=Well, Location=Location)
#     }
#   }, version=dataset[["Version"]], path=path)
#   rdfns <- rbindlist(rdfnL)
#   
#   ifnL <- lapply(dataset[grepl("Barcode",names(dataset))], function(barcode, path){
#     if(!length(dir(paste0(path,barcode,"/Analysis"),pattern = "imageID"))==0){
#       Path <- dir(paste0(path,barcode,"/Analysis"),pattern = "imageID", full.names = TRUE)
#       rdfns <- data.table(Barcode=barcode, Path=Path, Type="imageID", Well=NA, Location=NA)
#       
#     }
#   }, path)
#   ifns <- rbindlist(ifnL)
#   
#   mdfnL <- lapply(dataset[grepl("Barcode",names(dataset))], function(barcode, path){
#     if(!length(dir(paste0(path,barcode,"/Analysis"),pattern = "json|xlsx|an2omero"))==0){
#       rdfns <- data.table(Barcode=barcode, Path=dir(paste0(path,barcode,"/Analysis"),pattern = "xlsx|an2omero", full.names = TRUE), Type="metadata", Well=NA, Location=NA)
#     }
#   }, path=path)
#   mdns <- rbindlist(mdfnL)
#   
#   gfnL <- lapply(dataset[grepl("Barcode",names(dataset))], function(barcode, path){
#     if(!length(dir(paste0(path,barcode,"/Analysis"),pattern = "gal"))==0){
#       rdfns <- data.table(Barcode=barcode, Path=dir(paste0(path,barcode,"/Analysis"),pattern = "gal", full.names = TRUE), Type="gal", Well=NA, Location=NA)
#     }
#   }, path=path)
#   gfns <- rbindlist(gfnL)
#   
#   xfnL <- lapply(dataset[grepl("Barcode",names(dataset))], function(barcode, path){
#     if(!length(dir(paste0(path,barcode,"/Analysis"),pattern = "xml"))==0){
#       rdfns <- data.table(Barcode=barcode, Path=dir(paste0(path,barcode,"/Analysis"),pattern = "xml", full.names = TRUE), Type="xml", Well=NA, Location=NA)
#     }
#   }, path=path)
#   xfns <- rbindlist(xfnL)
#   return(rbind(rdfns,ifns, mdns, gfns, xfns))
# }


#' Add the URLs for the Omero IDs on the OHSU Omero server
#' @param dt A datatable with an ImageID column
#' @return the dt datatable with OmeroDetailURL, OmeroThumbnailURL and OmeroImageURL columns added
#' @export
addOmeroIDs <- function(dt){
  if(any(grepl("^ImageID$",colnames(dt)))){
    dt <- dt[,OmeroDetailURL := paste0('<a href="https://meplincs.ohsu.edu/webclient/img_detail/',ImageID,'/"',' target="_blank">Omero</a>')]
    dt <- dt[,OmeroThumbnailURL := paste0('<a href="https://meplincs.ohsu.edu/webclient/render_thumbnail/',ImageID,'/"',' target="_blank">Omero</a>')]
    dt <- dt[,OmeroImageURL := paste0('<a href="https://meplincs.ohsu.edu/webclient/render_image/',ImageID,'/"',' target="_blank">Omero</a>')]
  }
  return(dt)
}

# #' Median summarize values
# #' @param dt A datatable with Ligand,ECMp,MEP, Well and signal value columns that contain CP and/or
# #' PA in their names
# #' @return 
# #' @export
# medianSummarizel4 <- function(dt){
#   #median summarize the MEP replicates in a datatable
#   #Select all signal names and the minimal metadata
#   dtNames <- grep("^Ligand$|^ECMp$|MEP|^Well$|_CP_|_PA_",x = names(dt), value = TRUE)
#   
#   #Create a datatable of the signals and minimal metadata
#   dtKeep <- dt[, dtNames, with = FALSE]
#   dtMedians <- dtKeep[, lapply(.SD, numericMedian), keyby = "Ligand,ECMp,MEP,Well"]
#   return(dtMedians)
# }
# 
# #' @export
# #' 
# meanSummarizel4 <- function(dt){
#   #median summarize the MEP replicates in a datatable
#   #Select all signal names and the minimal metadata
#   dtNames <- grep("^Ligand$|^ECMp$|MEP|^Well$|_CP_|_PA_",x = names(dt), value = TRUE)
#   
#   #Create a datatable of the signals and minimal metadata
#   dtKeep <- dt[, dtNames, with = FALSE]
#   dtMeans <- dtKeep[, lapply(.SD, numericMean), keyby = "Ligand,ECMp,MEP,Well"]
#   return(dtMeans)
# }
# 
# #' @export
# #' 
# meanSummarizel3 <- function(dt){
#   #median summarize the Spot replicates in a datatable
#   #Select all signal names and the minimal metadata
#   dtNames <- grep("^Ligand$|^ECMp$|MEP|^Well$|PrintSpot|_CP_|_PA_",x = names(dt), value = TRUE)
#   
#   #Create a datatable of the signals and minimal metadata
#   dtKeep <- dt[, dtNames, with = FALSE]
#   dtMeans <- dtKeep[, lapply(.SD, numericMean), keyby = "Ligand,ECMp,MEP,Well,PrintSpot"]
#   return(dtMeans)
# }
# 
# numericMean <- function(x){
#   as.numeric(mean(x, na.rm=TRUE))
# }
# 
#' Back Transform Logit values
#' @export
# btLogit <- function(x){
#   2^x/(1+2^x)
# }
