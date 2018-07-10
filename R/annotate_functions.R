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
  #Temp fix for repeated column
  repeatedColumn <- colnames(dt) %>%
str_which("Ecm-2-timeUnit")
  if(length(repeatedColumn==2))  colnames(dt)[repeatedColumn[2]] <- "Ecm-3-timeUnit"
  #Temp fix for no ligand metadata
  dt$`Ligand-1`[dt$`Ligand-1`==""] <- "PBS"
  ########
  #Assign WellIndex values
  setnames(dt, "iWell", "Well")
  setkey(dt, Well)
  wellRows <- dt$welllayout %>%
    str_remove("x.*") %>%
    unique() %>%
    as.integer()
  wellCols <- dt$welllayout %>%
    str_remove(".*x") %>%
    unique() %>%
    as.integer()
  dt$Well <- wellAN(wellRows,wellCols)[dt$Well]
  wi <- data.table(Well = unique(dt$Well), WellIndex = 1:length(unique(dt$Well)))
  dt <- merge(dt,wi)
  #Rename to preprocessing pipeline variable names
  #setnames(dt,"runid","Barcode")
  dt <- dt %>%
    rename(Barcode=runid,
           CellLine=`CellLine-1`,
           ECM1=`Ecm-1`,
           ECM2=`Ecm-2`,
           ECM3=`Ecm-3`,
           Ligand1=`Ligand-1`,
           Drug1=`Drug-1`,
           Drug1ConcUnit=`Drug-1-concUnit`,
           Drug1Conc=`Drug-1-conc`,
           Spot=iSpot) %>%
    mutate(Barcode=str_remove(Barcode,"^mema-"))
  
  #Format stain endpoint metadata by parsing the primary and secondary antibody assignments
  #Identify the stain numbers that are animals and therefore secondaries
  stainRecSetNames <- unique(dt$`Stain-recordSetCollapsed`) %>%
    str_split( fixed( " + ")) %>%
    unlist() %>%
    str_remove("_.*")
  #Identify which stains are secondaries 
  stainRecSetSecIndices <- stainRecSetNames %>%
    str_which("Mouse|Rat|Donkey|Rabbit")
  #Count how many stains are in the set
  stainRecSetIndices <- 1:(str_count(unique(dt$`Stain-recordSetCollapsed`),fixed(" + ")) + 1)
  #Identify the stains that are not secondaries
  stainRecSetStains <- stainRecSetNames[-stainRecSetSecIndices]
  stainRecSetStainIndices <- stainRecSetIndices[-stainRecSetSecIndices]
  #make a one row dataframe with Endpointxxx columns that have the biomarkers as values
  stainWavelength <- dt %>%
    select(paste0("Stain-",stainRecSetStainIndices,"-wavelengthNm")) %>%
    distinct() %>%
    gather(value = "Endpoint") %>%
    transmute(Stain=str_extract(key,"[[:digit:]]") %>%
             as.numeric(),
           Biomarker=stainRecSetNames[Stain],
           Endpoint=Endpoint) %>%
    select(-Stain) %>%
    spread(key = Endpoint, value = Biomarker, sep="") %>%
    rename(EndpointDAPI=Endpoint395)
  
  dt <- data.table(dt,stainWavelength)
  
  #Shorten and combine Annot names
  dt$CellLine <- gsub("_.*","",dt$CellLine)
  
  #Reorder the ECMs so that if COL1 is present and paired, it's listed as ECM2  
  swapECM1 <- str_which(dt$ECM2,"[^(COL1_go0005584)&^(COL4_go0005587)]")
  ECM1Hold <- dt$ECM1[swapECM1]
  dt$ECM1[swapECM1] <- dt$ECM2[swapECM1]
  dt$ECM2[swapECM1] <-ECM1Hold
  #Reorder the ECMs so that if NID1 is present and paired, it's listed as ECM1  
  swapECM <- str_which(dt$ECM3,"[(NID1|1_P14543|1)]")
  ECM1Hold <- dt$ECM1[swapECM]
  dt$ECM1[swapECM] <- dt$ECM3[swapECM]
  dt$ECM3[swapECM] <-ECM1Hold
  
  dt$ECM1 <- compressHA(dt$ECM1)
  dt$ECM2 <- compressHA(dt$ECM2)
  dt$ECM3 <- compressHA(dt$ECM3)
  
  dt$ECMp <-paste0(gsub("_.*","",dt$ECM1),"_",gsub("_.*","",dt$ECM2),"_",gsub("_.*","",dt$ECM3)) %>%
    str_remove("_$|__$")
  #Chain ligands
  dt$Ligand <- str_remove(dt$Ligand1,"_.*")
  dt$MEP <- paste0(dt$ECMp,"_",dt$Ligand)
  dt$Drug <- gsub("_.*","",dt$Drug1)
  dt$MEP_Drug <-paste0(dt$MEP,"_",dt$Drug)
  dt$EndpointDAPI <-gsub("_.*","",dt$EndpointDAPI)
  dt$Endpoint488 <-gsub("_.*","",dt$Endpoint488)
  dt$Endpoint555 <-gsub("_.*","",dt$Endpoint555)
  dt$Endpoint647 <-gsub("_.*","",dt$Endpoint647)
  #Create ArrayRow and ArrayColumn variables
  dt$ArrayRow <- dt$iixii %>%
    str_remove("x.*") %>%
    str_remove("._") %>%
    as.integer()
  dt$ArrayColumn <- dt$iixii %>%
    str_remove(".*_") %>%
    as.integer()
  #Add a WellSpace spot index that recognizes the arrays are rotated 180 degrees
  dt$PrintSpot <- dt$Spot
  nrArrayRows <- max(dt$ArrayRow)
  nrArrayColumns <- max(dt$ArrayColumn)
  dt$PrintSpot[grepl("B", dt$Well)] <- (nrArrayRows*nrArrayColumns+1)-dt$PrintSpot[grepl("B", dt$Well)]
  #Add a drug concentration in uM/L
  ConcFactor <- str_extract(dt$Drug1ConcUnit,".*_per_.|.molar|volume_percent") %>%
    sapply(., function(x){
      switch(x,
             mol_per_L = 1e6,
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
  dt$Drug1Conc <- dt$Drug1Conc*ConcFactor
  dt$Drug1Conc[is.na(dt$Drug1Conc)] <- 0
  return(dt)
}

#' Read in and merge the Omero URLs
#' 
#' Adds Omero image IDs based on the WellName values
#' 
#' @param path The path to the file named barcode_imageIDs.tsv
#' @return a datatable with WellIndex, ArrayRow, ArrayColumn and ImageID columns
#' @export
getOmeroIDs <- function(path){
  dt <- fread(path)[,list(WellName,Row,Column,ImageID)]
  #Extract well index and convert to alphanumeric label
  if(grepl("nd2",dt$WellName[1])){
    dt[,WellIndex := gsub(".*_Well","",WellName)]
    dt[,WellIndex := gsub("_.*","",WellIndex)]
    dt[,WellIndex := (match(str_extract(WellIndex,"[[:alpha:]]"),LETTERS)-1)*12+as.integer(str_extract(WellIndex,"[[:digit:]]+"))]
  } else {
    dt[,WellIndex := as.integer(gsub(".*_Well","",WellName))]
  }
  setnames(dt,"Row","ArrayRow")
  setnames(dt,"Column","ArrayColumn")
  dt[,WellName := NULL]
  return(dt)
}
