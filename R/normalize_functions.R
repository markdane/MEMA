#Normalization functions for processing MEMAs

#' Median normalize a vector
#' 
#' @param x numeric vector
#' 
#' @return A numeric vector that is \code{x} divided by its median
#' @export
#' 
medianNorm <- function(x){
  xMedian <- median(x, na.rm=TRUE)
  normedValues <- x/xMedian
}

# #' Normalize the proliferation ratio signal to the collagen 1 values
# #' @param x a dataframe or datatable with columns names ProliferatioRatio
# #' and ShortName. ShortName must include at least one entry of COL1 or COL I.
# #' @return The input dataframe of datatable with a normedProliferation column that has the ProliferationRatio values divided by the median collagen
# #' 1 proliferation value
# #' @export
# normProfToCol1 <- function(x){
#   col1Median <- median(x$ProliferationRatio[x$ShortName %in% c("COL1", "COL I")],na.rm = TRUE)
#   normedProliferation <- x$ProliferationRatio/col1Median
# }

# #' Normalize to a base MEP
# #'
# #' Normalizes one channel of values for all MEPs in a multi-well plate to one
# #' base MEP.
# #'
# #' @param DT A \code{data.table} that includes a numeric value column to be
# #'   normalized, a \code{ECMp} column that has the printed ECM names and a
# #'   \code{Growth.Factors} or \code{Ligand}column that has the growth factor names.
# #' @param value The name of the column of values to be normalized
# #' @param baseECM A regular expression for the name or names of the printed ECM(s) to be normalized against
# #' @param baseGF A regular expression for the name or names of the soluble growth factors to be normalized against
# #' @return A numeric vector of the normalized values
# #'
# #' @section Details: \code{normWellsWithinPlate} normalizes the value column of
# #'   all MEPs by dividing the median value of the replicates of the MEP that
# #'   is the pairing of baseECM  with baseGF.
# #'   @export
# normWellsWithinPlate <- function(DT, value, baseECM, baseGF) {
#   if(!c("ECMp") %in% colnames(DT)) stop(paste("DT must contain a ECMp column."))
#   if(!c(value) %in% colnames(DT)) stop(paste("DT must contain a", value, "column."))
#   if("Ligand" %in% colnames(DT)){
#     valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMp)  & grepl(baseGF,DT$Ligand)),value, with=FALSE]), na.rm = TRUE)
#   } else if (c("Growth.Factors") %in% colnames(DT)) {
#     valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMp)  & grepl(baseGF,DT$Growth.Factors)),value, with=FALSE]), na.rm = TRUE)
#   } else stop (paste("DT must contain a Growth.Factors or Ligand column."))
#   normedValues <- DT[,value,with=FALSE]/valueMedian
#   return(normedValues)
# }

# #' RZS Normalize a Column of Data
# #'
# #' \code{normRZSWellsWithinPlate} normalizes all elements of DT[[value]] by
# #' subtracting the median of DT[[value]] of all baseECM spots in the
# #' baseL wells, then divides the result by the MAD*1.48 of all baseECM spots in
# #' the baseL wells
# #'@param DT A datatable with value, baseECM and baseL, ECMp and
# #'Ligand columns
# #'@param value A single column name of the value to be normalized
# #'@param baseECM A single character string or a regular expression that selects
# #'the ECM(s) that are used as the base for normalization.
# #'@param baseL A single character string or a regular expression that selects
# #'the ligand used as the base for normalization.
# #'@return a vector of RZS normalized values
# #' @export
# #'
# normRZSWellsWithinPlate <- function(DT, value, baseECM, baseL) {
#   if(!"ECMp" %in% colnames(DT)) stop (paste("DT must contain an ECMp column."))
#   if(!"Ligand" %in% colnames(DT)) stop (paste("DT must contain a Ligand column."))
#   if(!c(value) %in% colnames(DT)) stop(paste("DT must contain a", value, "column."))
#   
#   valueMedian <- median(unlist(DT[(grepl(baseECM, DT$ECMp) & grepl(baseL,DT$Ligand)), value, with=FALSE]), na.rm = TRUE)
#   if (is.na(valueMedian)) stop(paste("Normalization calculated an NA median for",value, baseECM, baseL))
#   
#   valueMAD <- mad(unlist(DT[(grepl(baseECM, DT$ECMp)  & grepl(baseL,DT$Ligand)),value, with=FALSE]), na.rm = TRUE)
#   #Correct for 0 MAD values
#   valueMAD <- valueMAD+.01
#   normedValues <- (DT[,value,with=FALSE]-valueMedian)/valueMAD
#   return(normedValues)
# }

# #' Normalize selected values in a dataset on a plate basis
# #' 
# #' A wrapper function for \code{normRZSWellsWithinPlate} that selects the
# #' \code{_CP_|_QI_|_PA_|SpotCellCount|Lineage} columns of dt if they exist and 
# #' normalizes them on a plate basis
# #' @param dt A data.table with a \code{Barcode} column numeric values to be RZS normalized 
# #'  using all ECM proteins in the FBS well
# #' @return A datatable with the normalized values
# #' @export
# normRZSDataset <- function(dt){
#   parmNormedList <- lapply(grep("_CP_|_QI_|_PA_|SpotCellCount|Lineage",colnames(dt),value = TRUE), function(parm){
#     dt <- dt[,paste0(parm,"_RZSNorm") := normRZSWellsWithinPlate(.SD, value=parm, baseECM = ".*",baseL = "FBS"), by="Barcode"]
#     return(dt)
#   })
#   return(parmNormedList[[length(parmNormedList)]])
# }


#' Loess normalize an array using the spatial residuals
loessNorm <- function(Value,Residual,ArrayRow,ArrayColumn){
  dt <-data.table(Value=Value,Residual=Residual,ArrayRow=ArrayRow,ArrayColumn=ArrayColumn)
  lm <- loess(Residual~ArrayRow+ArrayColumn, dt, span=.7)
  dt$ResidualLoess<-predict(lm)
  dt <- dt[,ValueLoess := Value-ResidualLoess]
  return(ValueLoess = dt$ValueLoess)
}

#' Loess normalize values within an array
#'@export
loessNormArray <- function(dt){
  #Identify the Signal name
  signalName <- unique(dt$SignalName)
  setnames(dt,signalName,"Value")
  #Get the median of the replicates within the array
  dt <- dt[,mel := median(Value), by=c("BW","ECMp")]
  #Get the residuals from the spot median
  #disable loess normalization
  dt <- dt[,Residual := Value]
  #re-enable loess
  #dt <- dt[,Residual := Value - mel]
  #Subtract the loess model of each array's residuals from the signal
  #dt <- dt[, ValueLoess:= loessNorm(Value, Residual, ArrayRow, ArrayColumn), by="BW"]
  dt <- dt[, ValueLoess:= Value]
  setnames(dt,"ValueLoess", paste0(signalName,"Loess"))
  BW <- "BW"
  PrintSpot <- "PrintSpot"
  dt <- dt[,c(paste0(signalName,"Loess"),BW,PrintSpot), with=FALSE]
  setkey(dt,BW,PrintSpot)
  return(dt)
}
