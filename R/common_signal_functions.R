#' Normalize the common signals across studies
#' 
#' @param paths A list of paths to level 2 data, each element is a different study.
#' @param k The number of factors to use in the RUV normalization
#' @param verbose A logical for progress messages
#' @return A datatable of raw and normalized spot level data with its metadata
#' @export
preprocessCommonSignalsLevel3 <- function(paths, k=256L, verbose=FALSE){
  #Read in the level 2 datasets
  l2List <- lapply(paths, getSpotLevelData)
  # Combine the common signals into a datatable
  l2 <- combineSignals(l2List)
  #Clean the common signal datatable
  l2 <- cleanCommonSignalDt(l2)
  
  signalsMinimalMetadata <- grep("_SE",grep("_CP_|_PA_|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^Drug1Conc$|^ArrayRow$|^ArrayColumn$|^CellLine$",colnames(l2), value=TRUE), value=TRUE, invert=TRUE)
  
  #RUVLoess normalize all signals
  if(!k==0){
    if(verbose)  message("Normalizing ", unlist(studyNameList),"\n")
    nDT <- normRUVLoessResiduals(l2[,signalsMinimalMetadata, with = FALSE], k)
    nDT$NormMethod <- "RUVLoessResiduals"
    l2$k <- k
    l3 <- merge(l2, nDT, by = c("BW","PrintSpot"))
  } else {
    l3 <- copy(l2)
    l3$NormMethod <- "none"
    l3$k <- k
  }
  
  #Add QA flags to the data
  l3 <- QASpotLevelData(l3, lowSpotCellCountThreshold=5,
                        lowRegionCellCountThreshold = 0.4,
                        lowWellQAThreshold = .7)
  
  return(l3)
}

#' Combine the common signals into a datatable
#' 
#' @param dtL a list of 1 to 3 datatables
#' @return a datatable with the columns that are common across all datatables in dtL
combineSignals <- function(dtL){
  #return a datatable with the common signals from 3 or fewer datatable
  if(length(dtL)==1) {
    dt <- dtL[[1]]
  } else if(length(dtL)==2) {
    commonSignals <- intersect(colnames(dtL[[1]]),colnames(dtL[[2]]))
    dt <- rbind(dtL[[1]][,commonSignals, with=FALSE],
                dtL[[2]][,commonSignals, with=FALSE])
  } else if(length(dtL)==3) {
    commonSignals <- intersect(intersect(colnames(dtL[[1]]),colnames(dtL[[2]])),colnames(dtL[[3]]))
    dt <- rbind(dtL[[1]][,commonSignals, with=FALSE],
                dtL[[2]][,commonSignals, with=FALSE],
                dtL[[3]][,commonSignals, with=FALSE])
  } else stop("Trying to find common signals for more than 3 studies")
  return(dt)
}

#' Clean the common signal datatable
#' 
#' Remove normalized signals and low quality data
#' @param dt A datatable
#' @return The input datatable without low quality data and normalized features
cleanCommonSignalDt <- function(dt){
  #Remove normalized and spot level signals
  dt <- dt[,grep("Norm|_SE",colnames(dt),invert=TRUE, value=TRUE), with=FALSE]
  #Remove blank spot data
  dt <- dt[!grepl("PBS|blank|Blank",dt$ECMp)]
  #Remove data from low quality wells
  # dt <-dt[!dt$QA_LowWellQA,]
  return(dt)
}