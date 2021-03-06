#QA Functions for MEMAs

#' Calculate Coefficient of Variation (CV)
#'
#'Calculates the CV of a numeric vector as the standard deviation over the mean. All NA values are removed from the input vector.
#'@param x A numeric vector.
#'@return The CV of the non-na values in the x.
#'
#'  @export
CV <- function(x){
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
}

#' Calculate a QA score for MEMAs
#'
#' Return a QA score based on the number of spots missing or below a threshold
#' @param DT A datatable with MEMA data
#' @param threshold A numeric vector of length 1
#' @param maxNrSpot A numeric vector of length 1 that is the maximum number of spots printed on the MEMA
#' @param value A character vector of length 1 that is the column name of the values that determine the QA score
#' @return A single numeric value that is the proportion of spots below the threshold
#'
#' @export
calcQAScore <- function(DT, threshold, maxNrSpot=700, value){
  QAScore <- (nrow(DT)-sum(DT[,value,with=FALSE] < threshold))/maxNrSpot
  return (QAScore)
}

#' Add QA flags to the spot level data
#' 
#' @param dt A datatable to be QA'd
#' @param lowSpotCellCountThreshold Threshold for spots with not enough cells
#' @param lowRegionCellCountThreshold Threshold for loess regions with not enough cells
#' @param lowWellQAThreshold Threshold for low quality wells
#' @return The same datatable with the QA columns
#' @export
QASpotLevelData <- function(dt,lowSpotCellCountThreshold=5,
                            lowRegionCellCountThreshold = 0.4,
                            lowWellQAThreshold = .7){
  #Low cell count spots
  dt$QA_LowSpotCellCount <- dt$Spot_PA_SpotCellCount < lowSpotCellCountThreshold
  dt$QA_lowSpotCellCountThreshold <- lowSpotCellCountThreshold
  #Low quality DAPI
  dt$QA_LowDAPIQuality <- FALSE
  #Flag spots below automatically loess QA threshold
  dt$QA_LowRegionCellCount <- dt$Spot_PA_LoessSCC < lowRegionCellCountThreshold
  dt$QA_lowRegionCellCountThreshold <- lowRegionCellCountThreshold
  #Flag wells below automatically calculated QA threshold
  dt$QA_LowWellQA <- FALSE
  dt$QA_LowWellQA[dt$QAScore < lowWellQAThreshold] <- TRUE
  dt$QA_lowWellQAThreshold <- lowWellQAThreshold
  return(dt)
}