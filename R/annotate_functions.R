#Functions to support MEMAs printed in 8 well plates

#' Rotate the metadata 180 degrees in Array space
#'
#'@param DT A data.table of metadata with Spot, ArrayRow and ArrayColumn columns.
#'@return The same data.table rotated 180 degrees in array space. The ArrayRow, Arraycolumn and Spot values are updated.
#'
#' @export
rotateMetadata <- function(DT){
  DT$ArrayRow <- max(DT$ArrayRow)+1-DT$ArrayRow
  DT$ArrayColumn <- max(DT$ArrayColumn)+1-DT$ArrayColumn
  DT$Spot <- as.integer(max(DT$Spot)+1-DT$Spot)
  return(DT)
}

#' Read in and parse an Aushon XML log file
#'
#' @param logFile An Aushon logfile
#' @return A datatable keyed by Row and Column with Depositions and
#'  PrintOrder columns.
#'
#' @export
readLogData<-function(logFile){
  #browser()
  data<-XML::xmlParse(logFile)
  dataList<-XML::xmlToList(data)
  #Only keep the sample attributes
  dataList<-dataList[names(dataList)=="Sample"]
  #Bind the XML data into a data table
  data<-data.table::rbindlist(dataList)
  #Create Row and Column data by shifting the values by 1
  data$Row<-as.integer(data$SpotRow)+1
  data$Column<-as.integer(data$SpotCol)+1
  #Convert deposition to an integer
  data$Depositions<-as.integer(data$Depositions)
  #Remove the 0 deposition entries
  data<-data[data$Depositions!=0,]
  #Create a print order column
  data$PrintOrder<-1:nrow(data)
  data.table::setkey(data,"PrintOrder")
  #Remove unneeded columns
  data <- data[,c("Row","Column","PrintOrder","Depositions"), with=FALSE]
  #Rotate by 90 degrees CW to match gal file orientation
  tmp <- data$Row
  data$Row <- data$Column
  data$Column <- 1+max(tmp)-tmp
  DT <- data.table::data.table(data,key="Row,Column")
  return(DT)
}

#' Convert column names in a data.table
#'
#' @param DT A data.table
#'
#' @return DT The same data.table with duplicated columns, invalid column name characters and trailing spaces removed.
#'
#' @export
convertColumnNames <- function (DT) {
  #Delete any duplicate names keeping the first instance
  DT <- DT[, unique(colnames(DT)), with = FALSE]
  #Replace invalid characters with a '.'
  data.table::setnames(DT, colnames(DT), make.names(colnames(DT)))
  #Remove all '.'s
  data.table::setnames(DT, colnames(DT), gsub("[.]", "", colnames(DT)))
}

#' Return the median of a vector as a numeric value
#'
#'\code{numericMedian} is a helper function for use within data.table that ensure the all medians are returned as numeric instead of numeric or integer values.
#' @param x integer or double vector
#'
#' @return The median of x as a numeric value
#'
#'@export
numericMedian <- function(x) as.numeric(median(x, na.rm = TRUE))


#' Process the metadata in an an2omero file
#' @param fileName The full name with path for the an2omero file
#' @return A datatable where each row is a unique spot in one well plate 
#'@export
processan2omero <- function (fileName) {
  #Process the file
  dt <- fread(fileName,header = TRUE)
  #Assign WellIndex values
  setkey(dt, Well)
  wi <- data.table(Well = unique(dt$Well), WellIndex = 1:length(unique(dt$Well)))
  dt <- merge(dt,wi)
  #Rename to preprocessing pipeline variable names
  setnames(dt,"OSpot","Spot")
  setnames(dt,"PlateID","Barcode")
  dt$EndpointDAPI <- dt[["395nm"]]
  dt$Endpoint488 <- dt[["488nm"]]
  dt$Endpoint555 <- dt[["555nm"]]
  dt$Endpoint647 <- dt[["640nm"]]
  dt$Endpoint750 <- dt[["750nm"]]
  #Shorten and combine Annot names
  dt$CellLine <- gsub("_.*","",dt$CellLine)
  dt$ECM1 <- compressHA(dt$ECM1)
  dt$ECM2 <- compressHA(dt$ECM2)
  dt$ECM3 <- compressHA(dt$ECM3)
  #Chain ECM proteins if the second one is not COL1
  dt$ECMp <-paste0(gsub("_.*","",dt$ECM1),"_",gsub("_.*","",dt$ECM2),"_",gsub("_.*","",dt$ECM3)) %>%
    gsub("_NA","",.) %>%
    gsub("_COL1|_$","",.)
  #Chain ligands
  dt$Ligand <-paste0(gsub("_.*","",dt$Ligand1),"_",gsub("_.*","",dt$Ligand2)) %>%
    gsub("_NA","",.)
  dt$MEP <- paste0(dt$ECMp,"_",dt$Ligand)
  dt$Drug <- gsub("_.*","",dt$Drug1)
  dt$MEP_Drug <-paste0(dt$MEP,"_",dt$Drug)
  dt$EndpointDAPI <-gsub("_.*","",dt$EndpointDAPI)
  dt$Endpoint488 <-gsub("_.*","",dt$Endpoint488)
  dt$Endpoint555 <-gsub("_.*","",dt$Endpoint555)
  dt$Endpoint647 <-gsub("_.*","",dt$Endpoint647)
  #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
  dt$PrintSpot <- dt$Spot
  nrArrayRows <- max(dt$ArrayRow)
  nrArrayColumns <- max(dt$ArrayColumn)
  dt$PrintSpot[grepl("B", dt$Well)] <- (nrArrayRows*nrArrayColumns+1)-dt$PrintSpot[grepl("B", dt$Well)]
  #Add a drug concentration in uM/L
  ConcFactor <- str_extract(dt$Drug1AnConcUnit,".*_per_.|.molar|volume_percent") %>%
    sapply(., function(x){
      switch(x,
             mmol_per_L = 1e3,
             umol_per_L = 1,
             nmol_per_L = 1e-3,
             pmol_per_L = 1e-6,
             mmolar = 1e6,
             umolar = 1e3,
             nmolar = 1,
             pmolar = 1e-3,
             0)
    })
  #Error handling for when there is no drug metadata
  dt$Drug1Conc <- dt$Drug1AnConc*ConcFactor
  dt$Drug1Conc[is.na(dt$Drug1Conc)] <- 0
  return(dt)
}

#' Read in and merge the Omero URLs
#' 
#' Adds Omero image IDs based on the WelIndex values
#' 
#' @param barcodePath The path to the file named barcode_imageIDs.tsv
#' @return a datatable with WellIndex, ArrayRow, ArrayColumn and ImageID columns
#' @export
getOmeroIDs <- function(barcodePath){
  barcode <- gsub(".*/","",barcodePath)
  dt <- fread(paste0(barcodePath,"/Analysis/",barcode,"_imageIDs.tsv"))[,list(WellName,Row,Column,ImageID)]
  #Extract well index and convert to alphanumeric label
  #8 and 96 well plates have different layouts in the WellName column
  if(grepl("Field",dt$WellName[1])){ #process as 96 well plate
    wellRow <- stringr::str_extract(dt$WellName,".-") %>%
      str_replace("-","") %>% 
      match(LETTERS)
    wellColumn <-  stringr::str_extract(dt$WellName,"-.*,") %>%
      str_replace("-","") %>% 
      str_replace(",","") %>% 
      as.integer()
    dt$WellIndex <- 12*(wellRow-1)+wellColumn
  } else { #process as 8 well plate
    dt[,WellIndex := as.integer(gsub(".*_Well","",WellName))]
    
  }
  setnames(dt,"Row","ArrayRow")
  setnames(dt,"Column","ArrayColumn")
  dt[,WellName := NULL]
  return(dt)
}