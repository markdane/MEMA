
#'Read the well metadata from a multi-sheet Excel file
#'
#' @param fn full path and file name for the well metadata file
#' @return a datatable with the well level metadata 
#' @export 
getWellMetadata <- function(fn){
  if(!length(fn)==1) stop(paste("There must be 1 xlsx metadata file in the",barcode, "Analysis folder"))
  data.table(readMetadata(fn), key="Well")
}

#'Get deposition values from an Aushon log dta files
#'
#' @param fn full path and file name for the log metadata file
#' @return a datatable with an integer Depostion column
#' @export 
getLogMetadata <- function(fn){
  if(!length(fn)==1) stop(paste("There must be 1 xml file in the",barcode, "Analysis folder"))
  readLogData(fn)
}
#'Read the spot metadata from a gal file
#'
#' @param fn full path and file name for the spot metadata file
#' @return a datatable with the spot level metadata 
#' @export 
getSpotMetadata <- function(fn){
  if(!length(fn)==1) stop(paste("There must be 1 gal file in the",barcode, "Analysis folder"))
  smd <- readSpotMetadata(fn)
  setnames(smd, "Name", "ECMp")
  return(smd)
}
#'Create a single datatable for all wells in an 8 well plate using!An! metadata
#'
#'8 well MEMAs are printed with the B row rotated 180 degrees from the A row. This routine rotates
#'the B rows and creates. The Spot values are indices in the physical coordinate system.
#'PrintSpot values are indices in the rotated printing coordinate system.
#'@param spotMetadata A datatable of spot level metadata
#'@param wellMetadata A datatable of well level metadata
#'@return A spot level datatable that has the well level metadata properly merged
#'@export
mergeSpot8WellMetadata <- function(spotMetadata,wellMetadata){
  wmdL <- apply(wellMetadata,1,  function(x){
    if(grepl("A",x[["Well"]])){
      dt <- cbind(spotMetadata,data.frame(t(x), stringsAsFactors = FALSE))
      dt <- dt[,PrintSpot := Spot]
    } else {
      dt <- cbind(rotateMetadata(spotMetadata),data.frame(t(x), stringsAsFactors = FALSE))
      dt <- dt[,PrintSpot := max(Spot)+1-Spot]
    }
  })
  rbindlist(wmdL)
}
#' Merge the spot well level metadata for a 96 well plate
#' 
#' This function assumes a 12 pin printhead in 6 row and 2 column configuration.
#' The printhead is oriented 90 ccw from the well plate. The 12 well pattern is
#' duplicated 8 times in the well plate to print 8 identical sets of 12 arrays
#'@param spotMetadata A datatable of spot level metadata
#'@param wellMetadata A datatable of well level metadata
#'@return A spot level datatable that has the well level metadata properly merged
#'@export
mergeSpot96WellMetadata <- function(spotMetadata,wellMetadata){
  #The ArrayRow and ArrayColumn indices are oriented with A01 well in the upper left
  #These are coordinates within each well
  #The gal file coordinates are rotated 90 degrees ccw from the array coordinates
  nrArrayRow <- max(spotMetadata$Column)
  nrArrayColumn <- max(spotMetadata$Row)
  spotMetadata <- spotMetadata[,ArrayRow := (nrArrayColumn+1-Column)]
  spotMetadata <- spotMetadata[,ArrayColumn := Row]
  #Assign Spots as sequential indices within each array
  spotMetadata <- spotMetadata[,Spot := ArrayColumn+(ArrayRow-1)*nrArrayRow]
  #PrintSpots are the same as spot in these plates
  spotMetadata <- spotMetadata[,PrintSpot := Spot]
  #Copy the 12 pin spot metadata 8 times to fill the 96 well plate
  smdPlate <- rbindlist(lapply(1:(length(unique(wellMetadata$Well))/length(unique(spotMetadata$Block))),function(x){
    cbind(spotMetadata,PlatePrintIndex = x)
  }))
  ###Assign PrintHead rows and columns in order to define wells
  #PlatePrint coordinates are for the print head sectors in the 96 wellplate
  nrPlatePrintCol <- 2
  nrPlatePrintRow <- max(smdPlate$PlatePrintIndex)/nrPlatePrintCol
  smdPlate$PlatePrintRow <- ceiling(smdPlate$PlatePrintIndex/nrPlatePrintCol)
  smdPlate$PlatePrintCol <- ((smdPlate$PlatePrintIndex-1) %% nrPlatePrintCol)+1
  #PrintHead coords are the rows and columns of pins/blocks within the print head
  nrPrintHeadCol <- 2
  nrPrintHeadRow <- max(smdPlate$Block)/nrPrintHeadCol
  smdPlate$PrintHeadRow <- ceiling(smdPlate$Block/nrPrintHeadCol)
  smdPlate$PrintHeadCol <- ((smdPlate$Block-1) %% nrPrintHeadCol)+1
  #Use PlatePrint and PrintHead to determine WellIndex
  smdPlate$PlateRow <- (smdPlate$PlatePrintRow-1)*nrPrintHeadCol+(3-smdPlate$PrintHeadCol)
  smdPlate$PlateCol <- smdPlate$PlateCol <- (smdPlate$PlatePrintCol-1)*nrPrintHeadRow+smdPlate$PrintHeadRow
  nrPlateRow <- max(smdPlate$PlateRow)
  nrPlateCol <- max(smdPlate$PlateCol)
  smdPlate <- smdPlate[,WellIndex := (PlateRow-1)*nrPlateCol+PlateCol]
  smdPlate <-smdPlate[,Well:=wellAN(nrPlateRow, nrPlateCol)[smdPlate$WellIndex]]
  #Cleanup working columns used to create Well values
  smdPlate <- smdPlate[,PlatePrintRow :=NULL]
  smdPlate <- smdPlate[,PlatePrintCol :=NULL]
  smdPlate <- smdPlate[,PrintHeadRow :=NULL]
  smdPlate <- smdPlate[,PrintHeadCol :=NULL]
  smdPlate <- smdPlate[,WellIndex :=NULL]
  #Merge in the well metadata
  metadata <- merge(smdPlate,wellMetadata,by="Well")
}

#' Get metadata from either An! or !An! files
#' @param metadataFiles A list of of full file names for metadata. The items in 
#' the list must have names of annotMetadata, logMetadata, spotMetadata or wellMetadata
#' @param useAnnotMetadata Logical indicating whether to use metadata in Annot files
#' or excel files
#' @return A datatable containing the metadata
#' 
#' @export
getMetadata <- function(metadataFiles, useAnnotMetadata=TRUE){
  #Use metadata from an2omero files
  if(useAnnotMetadata){
    metadata <- processan2omero(metadataFiles[["annotMetadata"]])
    MEMA8Well <- length(unique(metadata$Well))==8
    MEMA96Well <- length(unique(metadata$Well))==96
  } else { #Process xml, gal and excel files to get all metadata
    #Read the log metadata from an Aushon log file
    ldf <- getLogMetadata(metadataFiles[["logMetadata"]])
    #Read in the spot metadata from the gal file
    galmd <- getSpotMetadata(metadataFiles[["spotMetadata"]])
    spotMetadata <- merge(galmd,ldf, by = c("Row","Column"), all=TRUE)
    #Read the well metadata from a multi-sheet Excel file
    wellMetadata <- getWellMetadata(metadataFiles[["wellMetadata"]])
    #Determine if this is metadata for an 8 or 96 well plate 
    MEMA8Well <- setequal(unique(wellMetadata$Well),c("A01","A02","A03","A04","B01","B02","B03","B04"))
    MEMA96Well <- length(unique(wellMetadata$Well))==96
    if(MEMA8Well){
      #Create a single datatable for all wells in an 8 well plate using!An! metadata
      metadata <- mergeSpot8WellMetadata(spotMetadata,wellMetadata)
    } else if (MEMA96Well) {
      # Merge the spot well level metadata for a 96 well plate
      metadata <- mergeSpot96WellMetadata(spotMetadata,wellMetadata)
    } else {
      stop("Only 8 well and 96 well plates are supported")
    }
    #Add MEP and convenience labels for wells and ligands
    metadata <- metadata[,MEP:=paste(ECMp,Ligand,sep = "_")]
    metadata <- metadata[,Well_Ligand:=paste(Well,Ligand,sep = "_")]
    metadata <- metadata[,MEP_Drug:=paste(MEP,Drug,sep = "_")]
    metadata <- metadata[,Barcode := barcode]
    # Eliminate Variations in the Endpoint metadata
    endpointNames <- grep("End",colnames(metadata), value=TRUE)
    endpointWL <- regmatches(endpointNames,regexpr("[[:digit:]]{3}|DAPI",endpointNames))
    setnames(metadata,endpointNames,paste0("Endpoint",endpointWL))
    #Create short display names, then replace where not unique
    #Use entire AnnotID for ligands with same uniprot IDs
    # metadata$Ligand[grepl("NRG1",metadata$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = metadata$Ligand[grepl("NRG1",metadata$Ligand)])
    # metadata$Ligand[grepl("TGFB1",metadata$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = metadata$Ligand[grepl("TGFB1",metadata$Ligand)])
    # metadata$Ligand[grepl("CXCL12",metadata$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = metadata$Ligand[grepl("CXCL12",metadata$Ligand)])
    # metadata$Ligand[grepl("IGF1",metadata$Ligand)] <- simplifyLigandAnnotID(ligand = "IGF1",annotIDs = metadata$Ligand[grepl("IGF1",metadata$Ligand)])
  }
  return(metadata)
}

#' Get the data from the Cell Profiler pipeline
#' @export
getCPData <- function(dataBWInfo, curatedOnly=TRUE, curatedCols= "ImageNumber|ObjectNumber|AreaShape|_MedianIntensity_|_IntegratedIntensity_|_Center_|_PA_|Texture", verbose=FALSE){
  dtL <- lapply(unique(dataBWInfo$Well), function(well){
    if(verbose) message(paste("Reading and annotating data for",barcode, well,"\n"))
    nuclei <- convertColumnNames(fread(dataBWInfo$Path[grepl("Nuclei",dataBWInfo$Location)&grepl(well,dataBWInfo$Well)], verbose=FALSE, showProgress=FALSE))
    if (curatedOnly) nuclei <- nuclei[,grep(curatedCols,colnames(nuclei)), with=FALSE]
    setnames(nuclei,paste0("Nuclei_",colnames(nuclei)))
    setnames(nuclei,"Nuclei_CP_ImageNumber","Spot")
    setnames(nuclei,"Nuclei_CP_ObjectNumber","ObjectNumber")
    setkey(nuclei,Spot,ObjectNumber) 
    
    if(any(grepl("Cells",dataBWInfo$Location)&grepl(well,dataBWInfo$Well))){
      cells <- convertColumnNames(fread(dataBWInfo$Path[grepl("Cells",dataBWInfo$Location)&grepl(well,dataBWInfo$Well)], verbose=FALSE, showProgress=FALSE))
      if (curatedOnly) cells <- cells[,grep(curatedCols,colnames(cells)), with=FALSE]
      setnames(cells,paste0("Cells_",colnames(cells)))
      setnames(cells,"Cells_CP_ImageNumber","Spot")
      setnames(cells,"Cells_CP_ObjectNumber","ObjectNumber")
      setkey(cells,Spot,ObjectNumber)
    } 
    
    if(any(grepl("Cytoplasm",dataBWInfo$Location)&grepl(well,dataBWInfo$Well))){
      cytoplasm <- convertColumnNames(fread(dataBWInfo$Path[grepl("Cytoplasm",dataBWInfo$Location)&grepl(well,dataBWInfo$Well)], verbose=FALSE, showProgress=FALSE))
      if (curatedOnly) cytoplasm <- cytoplasm[,grep(curatedCols,colnames(cytoplasm)), with=FALSE]
      setnames(cytoplasm,paste0("Cytoplasm_",colnames(cytoplasm)))
      setnames(cytoplasm,"Cytoplasm_CP_ImageNumber","Spot")
      setnames(cytoplasm,"Cytoplasm_CP_ObjectNumber","ObjectNumber")
      setkey(cytoplasm,Spot,ObjectNumber)
    } 
    
    #Merge the data from the different locations if it exists
    if(exists("cells")&exists("cytoplasm")) {
      dt <- cells[cytoplasm[nuclei]]
    } else {
      dt <- nuclei
    }
    #Add the well name and barcode as parameters
    dt <- dt[,Well := well]
    dt <- dt[,Barcode := barcode]
    #Scale the intensity values
    intensityNames <- grep("Intensity",colnames(dt), value=TRUE)
    scaledInts <- dt[,intensityNames, with=FALSE]*2^16
    dt <- cbind(dt[,!intensityNames, with=FALSE],scaledInts)
    return(dt)
  })
}
#' Read and convert INCell data to CP format
#' @param cellDataFilePaths A character vector of the file paths to the cell level data
#' @export
getICData <- function(cellDataFilePaths, endPoint488, endPoint555, endPoint647, verbose=FALSE){
  if(verbose) message(paste("Reading and annotating INCell data for",barcode,"\n"))
  #Read and combine the 2 header rows after the summary information
  hdrRows <- read.csv(cellDataFilePaths,skip = 18, nrows=2, header=FALSE,stringsAsFactors = FALSE)
  hdr <- sub("^_","",paste(hdrRows[1,],hdrRows[2,],sep="_"))
  #Read the cell level and spot summary data
  df <- read.csv(cellDataFilePaths,skip = 20, header=FALSE, stringsAsFactors = FALSE)
  #remove the spot summary data and convert to a data.table
  if(any(which(df$V1==""))) {
    df <- df[-(min(which(df$V1=="")):nrow(df)),]
  }
  dt <- data.table(df)
  #Name the columns
  setnames(dt,names(dt), hdr)
  if("NA_NA" %in% colnames(dt)) dt <- dt[,NA_NA := NULL]
  #Create a spot column from the field value
  dt$Spot <- gsub(".*fld ","",dt$Well)
  dt$Spot <- as.numeric(gsub(")","",dt$Spot))
  #Convert well names to alphanumeric with 2 digit columns
  wellRow <- str_match(dt$Well,"[:alpha:]")
  wells <- str_match(dt$Well,"[:digit:][:digit:]?") %>%
    as.numeric() %>%
    sprintf("%02d",.) %>%
    paste0(wellRow,.)
  
  #Convert all columns besides the Well to numeric values
  dt <- dt[,lapply(.SD, as.numeric),.SDcols = grep("Well",colnames(dt),value=TRUE,invert=TRUE)]
  dt$Well <- wells
  # dt$PlateRow <- wellRow
  # dt <- dt[,PlateCol := as.numeric(gsub("[[:alpha:]]","",dt$Well))]
  #Convert INCell names to CP versions
  dt <- convertColumnNames(dt)
  #Assume first nuclear channel is DAPI
  #Create a list of lists with IC names and corresponding CP names if available
  ICtoCPNames <- list(
    list("Nuclei_NucIntensity","Nuclei_CP_Intensity_MedianIntensity_Dapi"),
    list("Nuclei_NucArea","Nuclei_CP_AreaShape_Area"),
    list("Nuclei_NuccgX","Nuclei_CP_AreaShape_Center_X"),
    list("Nuclei_NuccgY","Nuclei_CP_AreaShape_Center_Y"),
    list("Nuclei_NucElongation","Nuclei_CP_AreaShape_Eccentricity"),
    list("Nuclei_NucCellIntensity","Cell_CP_Intensity_MedianIntensity_Dapi"),
    list("Nuclei_IxANuc","Nuclei_CP_Intensity_IntegratedIntensity_Dapi"),
    list("Cells_CellIntensity",paste0("Cytoplasm_CP_Intensity_MedianIntensity_",endPoint488)),
    list("Reference1_CellIntensity",paste0("Cytoplasm_CP_Intensity_MedianIntensity_",endPoint555)),
    list("Reference2_NucIntensity",paste0("Nuclei_CP_Intensity_MedianIntensity_",endPoint647)),
    list("Cell","ObjectNumber")
  )
  #Change names within dt to match downstream CP pipeline
  foo <- sapply(ICtoCPNames,function(x) {
    if(x[[1]] %in% colnames(dt)) setnames(dt,x[[1]],x[[2]])
  })
  rm(foo)
  
  dt$Barcode <- barcode
  dtL <-list(dt)
}
#'Clean up issues that have been problems in previous pipelines
#'@param dt A datatable of MEMA data and metadata
#'@return The same datatable with legacy issues corrected
#'@export
cleanLegacyIssues <- function(dt){
  #Remove problematic features
  dt <- dt[,grep("Euler",colnames(dt),invert=TRUE), with=FALSE]
  #Change Edu back to EdU
  if(any(grepl("Edu",colnames(dt)))){
    edUNames <- grep("Edu",colnames(dt),value=TRUE)
    setnames(dt,edUNames,gsub("Edu","EdU",edUNames))
  }
  #Add the pin diameter metadata in microns
  if(any(grepl("MCF7|PC3|YAPC",unique(dt$CellLine)))){
    dt$PinDiameter <- 180
  } else {
    dt$PinDiameter <- 350
  }
  return(dt)
}

#` Filter our debris and cell clusters
#'@export
filterObjects <- function(dt,nuclearAreaThresh=50,nuclearAreaHiThresh=4000){
  if(any(grepl("Nuclei_CP_AreaShape_Area",colnames(dt)))){
    dt <- dt[dt$Nuclei_CP_AreaShape_Area > nuclearAreaThresh,]
    dt <- dt[dt$Nuclei_CP_AreaShape_Area < nuclearAreaHiThresh,]
  }
  return(dt)
}
#' Add local XY and polar coordinates
#' @export
addPolarCoords <- function(dt){
  if(any(grepl("Nuclei_CP_AreaShape_Center",colnames(dt)))){
    #Add the local polar coordinates
    dt <- dt[,Nuclei_PA_Centered_X :=  Nuclei_CP_AreaShape_Center_X-median(Nuclei_CP_AreaShape_Center_X), by=c("Well","Spot")]
    dt <- dt[,Nuclei_PA_Centered_Y :=  Nuclei_CP_AreaShape_Center_Y-median(Nuclei_CP_AreaShape_Center_Y), by=c("Well","Spot")]
    dt <- dt[, Nuclei_PA_AreaShape_Center_R := sqrt(Nuclei_PA_Centered_X^2 + Nuclei_PA_Centered_Y^2), by=c("Well","Spot")]
    dt <- dt[, Nuclei_PA_AreaShape_Center_Theta := calcTheta(Nuclei_PA_Centered_X, Nuclei_PA_Centered_Y), by=c("Well","Spot")]
  }
  return(dt)
}

#'Add spot level normalizations for median intensities
#'@export
spotNormIntensities <- function(dt){
  intensityNamesAll <- grep("_CP_Intensity_Median",colnames(dt), value=TRUE)
  intensityNames <- grep("Norm",intensityNamesAll,invert=TRUE,value=TRUE)
  for(intensityName in intensityNames){
    #Median normalize the median intensity at each spot
    setnames(dt,intensityName,"value")
    dt <- dt[,paste0(intensityName,"_SpotNorm") := medianNorm(value),by="Barcode,Well,Spot"]
    setnames(dt,"value",intensityName)
  }
  return(dt)
}
#' Add measures of adjacency based on cell positions in the population
#' 
#' @param neighborhoodNucleiRadii Defines the neighborhood annulus
#' @param neighborsThresh Gates sparse cells on a spot
#' @param wedgeAngs Size in degrees of spot wedges used in perimeter gating
#' @param outerThresh Defines outer cells used in perimeter gating
#'@export
calcAdjacency <- function(dt, neighborhoodNucleiRadii=7,neighborsThresh=0.4,wedgeAngs=20, outerThresh=0.5){
  densityRadius <- sqrt(median(dt$Nuclei_CP_AreaShape_Area, na.rm = TRUE)/pi)
  #Count the number of neighboring cells
  dt <- dt[,Nuclei_PA_AreaShape_Neighbors := cellNeighbors(.SD, radius = densityRadius*neighborhoodNucleiRadii), by = "Barcode,Well,Spot"]
  #Rules for classifying perimeter cells
  dt <- dt[,Spot_PA_Sparse := Nuclei_PA_AreaShape_Neighbors < neighborsThresh]
  #Add a local wedge ID to each cell based on conversations with Michel Nederlof
  dt <- dt[,Spot_PA_Wedge:=ceiling(Nuclei_PA_AreaShape_Center_Theta/wedgeAngs)]
  #Define the perimeter cell if it exists in each wedge
  #Classify cells as outer if they have a radial position greater than a thresh
  dt <- dt[,Spot_PA_OuterCell := labelOuterCells(Nuclei_PA_AreaShape_Center_R, thresh=outerThresh),by="Barcode,Well,Spot"]
  #Require a perimeter cell not be in a sparse region
  denseOuterDT <- dt[!dt$Spot_PA_Sparse  & dt$Spot_PA_OuterCell]
  denseOuterDT <- denseOuterDT[,Spot_PA_Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Spot_PA_Wedge"]
  setkey(dt,Barcode,Well,Spot,ObjectNumber)
  setkey(denseOuterDT,Barcode,Well,Spot,ObjectNumber)
  dt <- denseOuterDT[,list(Barcode,Well,Spot,ObjectNumber,Spot_PA_Perimeter)][dt]
  dt$Spot_PA_Perimeter[is.na(dt$Spot_PA_Perimeter)] <- FALSE
  return(dt)
}

#' Gate cells based on specific stains if they are in the dataset
#' 
#' @param dt A datatable that contains cell level data and metadata
#' @export
gateCells <- function(dt, synId="syn10138929"){
  #Read parameter file from Synapse
  library(synapseClient)
  synapseLogin()
  parms <- synGet(synId) %>%
    synapseClient::getFileLocation() %>%
    data.table::fread(check.names = TRUE) %>%
    dplyr::filter(Barcode==unique(dt$Barcode))
  overrideParms <- nrow(parms)==1
  #Set 2N and 4N DNA status
  #Use override parameter if it exists
  if(overrideParms) {
    if(!parms$CycleStateThresh==0){
      #if CycleStateThresh == 0 then do not generate a Cycle State signal
      dt$Nuclei_PA_Cycle_State <- 1
      dt$Nuclei_PA_Cycle_State[dt$Nuclei_CP_Intensity_IntegratedIntensity_Dapi > parms$CycleStateThresh] <- 2
    }
  } else {
    dt <- dt[,Nuclei_PA_Cycle_State := gateOnlocalMinima(Nuclei_CP_Intensity_IntegratedIntensity_Dapi)]
  }
  
  #Create staining set specific derived parameters
  if("Nuclei_CP_Intensity_MedianIntensity_EdU" %in% colnames(dt)){
    #Use override parameters if they exist
    if(overrideParms) {
      if(!parms$EdUPositiveThresh==0){
        #if EdUPositiveThresh == 0 then do not generate a Gated EdU Positive signal
        dt$Nuclei_PA_Gated_EdUPositive <- 0
        dt$Nuclei_PA_Gated_EdUPositive[dt$Nuclei_CP_Intensity_MedianIntensity_EdU > parms$EdUPositiveThresh] <- 1
      }
    } else {
      #Use the entire plate to set the autogate threshold if there's no control well
      if(any(grepl("FBS",unique(dt$Ligand)))){
        dt <- dt[,Nuclei_PA_Gated_EdUPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_EdU",ctrlLigand = "FBS"), by="Barcode"]
      } else {
        dt <- dt[,Nuclei_PA_Gated_EdUPositive := kmeansCluster(.SD,value =  "Nuclei_CP_Intensity_MedianIntensity_EdU",ctrlLigand = "."), by="Barcode"]
      }
      #Modify the auto gate threshold to be 3 sigma from the EdU- median
      EdUDT <- dt[Nuclei_PA_Gated_EdUPositive==0,.(EdUMedian=median(Nuclei_CP_Intensity_MedianIntensity_EdU, na.rm=TRUE),EdUSD=sd(Nuclei_CP_Intensity_MedianIntensity_EdU, na.rm=TRUE)),by="Barcode,Well"]
      EdUDT <- EdUDT[,Median3SD := EdUMedian+3*EdUSD]
      dt <- merge(dt,EdUDT,by=c("Barcode","Well"))
      dt$Nuclei_PA_Gated_EdUPositive[dt$Nuclei_CP_Intensity_MedianIntensity_EdU>dt$Median3SD]<-1
    }
  }
  
  #Gate KRT5 using override parameters if they exist
  if("Cytoplasm_CP_Intensity_MedianIntensity_KRT5" %in% colnames(dt)){  
    if(overrideParms) {
      if(!parms$KRT5PositiveThresh==0){
        #if KRT5PositiveThresh == 0 do not generate a KRT5 gated Positive signal
        dt$Cytoplasm_PA_Gated_KRT5Positive <- 0
        dt$Cytoplasm_PA_Gated_KRT5Positive[dt$Cytoplasm_CP_Intensity_MedianIntensity_KRT5 > parms$KRT5PositiveThresh] <- 1
      }
    } else if(grepl("HMEC",unique(dt$CellLine))) {
      dt <- dt[,Cytoplasm_PA_Gated_KRT5Positive := gateOnQuantile(Cytoplasm_CP_Intensity_MedianIntensity_KRT5,probs=.02),by="Barcode"]
    } else {
      dt <- dt[,Cytoplasm_PA_Gated_KRT5Positive := kmeansCluster(.SD,value =  "Cytoplasm_CP_Intensity_MedianIntensity_KRT5",ctrlLigand = "."), by="Barcode"]
    }
  }
  
  #Gate KRT14 using override parameters if they exist
  if("Cytoplasm_CP_Intensity_MedianIntensity_KRT14" %in% colnames(dt)){  
    if(overrideParms) {
      if(!parms$KRT14PositiveThresh==0){
        #if KRT14PositiveThresh == 0 do not generate a KRT14 gated Positive signal
        dt$Cytoplasm_PA_Gated_KRT14Positive <- 0
        dt$Cytoplasm_PA_Gated_KRT14Positive[dt$Cytoplasm_CP_Intensity_MedianIntensity_KRT14 > parms$KRT14PositiveThresh] <- 1
      }
    } else {
      dt <- dt[,Cytoplasm_PA_Gated_KRT14Positive := gateOnQuantile(Cytoplasm_CP_Intensity_MedianIntensity_KRT14,probs=.02),by="Barcode"]
    }
  }
  
  #Gate KRT19 using override parameters if they exist
  if("Cytoplasm_CP_Intensity_MedianIntensity_KRT19" %in% colnames(dt)){  
    if(overrideParms) {
      if(!parms$KRT19PositiveThresh==0){
        #if KRT19PositiveThresh == 0 do not generate a KRT19 gated Positive signal
        dt$Cytoplasm_PA_Gated_KRT19Positive <- 0
        dt$Cytoplasm_PA_Gated_KRT19Positive[dt$Cytoplasm_CP_Intensity_MedianIntensity_KRT19 > parms$KRT19PositiveThresh] <- 1
      }
    } else {
      dt <- dt[,Cytoplasm_PA_Gated_KRT19Positive := kmeansCluster(.SD, "Cytoplasm_CP_Intensity_MedianIntensity_KRT19",ctrlLigand = "."), by="Barcode"]
    }
  }
  
  #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
  if ("Cytoplasm_CP_Intensity_MedianIntensity_KRT19" %in% colnames(dt)&
      "Cytoplasm_CP_Intensity_MedianIntensity_KRT5" %in% colnames(dt)){
    dt <- dt[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
    
    #Determine the class of each cell based on KRT5 and KRT19 class
    #0 double negative
    #1 KRT5+, KRT19-
    #2 KRT5-, KRT19+
    #3 KRT5+, KRT19+
    dt <- dt[,Cytoplasm_PA_Gated_KRTClass := Cytoplasm_PA_Gated_KRT5Positive+2*Cytoplasm_PA_Gated_KRT19Positive]
  }
  
  if ("Cytoplasm_CP_Intensity_MedianIntensity_KRT14" %in% colnames(dt)&
      "Cytoplasm_CP_Intensity_MedianIntensity_KRT19" %in% colnames(dt)){
    #Calculate a lineage ratio of luminal/basal or KRT19/KRT14
    dt <- dt[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT14]
    #Determine the class of each cell based on KRT5 and KRT19 class
    #0 double negative
    #1 KRT14+, KRT19-
    #2 KRT14-, KRT19+
    #3 KRT14+, KRT19+
    dt <- dt[,Cytoplasm_PA_Gated_KRTClass := Cytoplasm_PA_Gated_KRT14Positive+2*Cytoplasm_PA_Gated_KRT19Positive]
  }
  
  if ("Cytoplasm_PA_Gated_KRT5Positive" %in% colnames(dt)&
      "Nuclei_PA_Gated_EdUPositive" %in% colnames(dt)){
    #Determine the class of each cell based on KRT5 and EdU class
    #0 double negative
    #1 KRT5+, EdU-
    #2 KRT5-, EdU+
    #3 KRT5+, EdU+
    #Use override parameters if they exist
    dt <- dt[,Cells_PA_Gated_EdUKRT5Class := Cytoplasm_PA_Gated_KRT5Positive+2*Nuclei_PA_Gated_EdUPositive]
  }
  
  if ("Cytoplasm_PA_Gated_KRT14Positive" %in% colnames(dt)&
      "Nuclei_PA_Gated_EdUPositive" %in% colnames(dt)){
    #Determine the class of each cell based on KRT5 and EdU class
    #0 double negative
    #1 KRT14+, EdU-
    #2 KRT14-, EdU+
    #3 KRT14+, EdU+
    #Use override parameters if they exist
    dt <- dt[,Cells_PA_Gated_EdUKRT14Class := Cytoplasm_PA_Gated_KRT14Positive+2*Nuclei_PA_Gated_EdUPositive]
  }
  
  if("Nuclei_PA_Gated_EdUPositive" %in% colnames(dt)&
     "Cytoplasm_PA_Gated_KRT19Positive" %in% colnames(dt)){
    #Determine the class of each cell based on KRT19 and EdU class
    #0 double negative
    #1 KRT19+, EdU-
    #2 KRT19-, EdU+
    #3 KRT19+, EdU+
    dt <- dt[,Cells_PA_Gated_EdUKRT19Class := Cytoplasm_PA_Gated_KRT19Positive+2*Nuclei_PA_Gated_EdUPositive]
  }
  
  return(dt)
}

#' Write the cell level raw data and metadata Level 1 data. 
#' @param dt A datatable of cell level data to be written to disk
#' @param path The path to write the file
#' @param barcode The barcode that will be encoded into the file name
#'@return  None called for the side effect of writing to disk
#'@export
writeCellLevel <- function(dt,path,barcode, verbose=FALSE){
  if(verbose) message(paste("Writing",barcode,"level 1 full file to disk\n"))
  writeTime<-Sys.time()
  fwrite(dt, paste0(path,barcode, "/Analysis/", barcode,"_","Level1.tsv"), sep = "\t", quote=FALSE)
  if(verbose) message(paste("Write time:", Sys.time()-writeTime,"\n"))
}  

#' Format 96 well plate raw data from the Cell Profiler Pipeline
#' 
#' Reformats the metadata for one row of a 96 well plate to match the pipeline. 
#' The Spot column is converted to an index within each well and a label for the Well is added
#' to the output data.table
#' 
#' @param dt data.table with a Spot column
#' @param nrArrayRows The number of rows in the printed array
#' @param nrArrayColumns The number of columns in the printed array
#' @return A data.table with an alphanumeric Well column and reindexed Spot column
#' @export
formatCP96Well <- function(dt, nrArrayRows, nrArrayColumns){
  nrArraySpots <-nrArrayRows*nrArrayColumns
  dt$WellIndex <- (dt$Spot-1) %/% nrArraySpots +1
  dt$Spot <- ((dt$Spot-1) %% nrArraySpots)+1
  dt$Well <-  dt$WellIndex %>%
    sprintf("%02d",.) %>%
    paste0(gsub("Row","",dt$Well),.)
  dt[,WellIndex := NULL]
  return(dt)
}
