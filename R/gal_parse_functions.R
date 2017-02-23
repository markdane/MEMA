#Functions to parse gal files
#Author Mark Dane 3/2/2015

#'Read the spot metadata from a gal file
#'
#'\code{readSpotMetadata} returns a dataframe of the contents in a gal file and
#'adds array indices.
#'
#'@param galFile the name of a gal file to be read
#'@return A dataframe with the contents of each spot in the gal file. The spots
#'  are converted from print block, row and column space to array row and column
#'  space. The array space is aligned with the print space. That is,
#'  the 1,1,1 position in the print space is the first row and first column of the
#'  array space. The Name and ID fields of spots that were not printed have their
#'  dashes replaced with the word blank.
#' @export
readSpotMetadata <- function(galFile) {
  #Read the GAL file
  #browser()
  df <- limma::readGAL(galFile)
  df$Name <- gsub("^-$","blank",df$Name)
  df$ID <- gsub("^-$","blank",df$ID)
  layout <- limma::getLayout(df)
  nrCols <- layout$nspot.c*layout$ngrid.c
  nrRows <- layout$nspot.r*layout$ngrid.r
  colnames(df)[colnames(df) == "Grid"] <- "Block"
  df<-addArrayPositionNoRotate(df,gridsPerRow = layout$ngrid.c)
  DT<-data.table::data.table(df,key=c("Block","Row","Column"))
  return(DT)
}

#'Read metadata from a multi-sheet excel file
#'
#' \code{readMetadata} reads well level metadata from an excel file. Each sheet in the file represents a different data type. The name of the sheet will become the name of the column. The contents in the sheet will become the data values. Each sheet is foramtted to match a well plate with the A01 well in the upper left corner. First row and column of each sheet are left as labels. The values start in the second row and column.
#' @param xlsFile The name of the excel file.
#' @return a data frame with the data values from each sheet in a column with the name of the sheet.
#' @export
readMetadata<-function(xlsFile){
  #browser()
  sheetList<-sapply(gdata::sheetNames(path.expand(xlsFile)), gdata::read.xls, xls = path.expand(xlsFile), simplify = FALSE,stringsAsFactors=TRUE,check.names=FALSE,row.names="Row/Column")
  nrRows<-dim(sheetList[[1]])[1]
  nrCols<-max(as.numeric(colnames(sheetList[[1]])))
  nrWells=nrRows*nrCols
  sheetDF<-data.frame(lapply(sheetList,function(df,nrCols){
    #create a dataframe from all rows and columns of each sheet
    dfM<-matrix(t(df[,1:nrCols]),byrow=TRUE)
  }, nrCols=nrCols),WellIndex=1:nrWells,Well=wellAN(nrRows,nrCols),check.names=TRUE, stringsAsFactors=FALSE)
  return(sheetDF)
}
