# Helper function to create level 3 and 4 data
#` Calculate standard error of the mean
#' 
#' Omit na values if present
#' @param x a numeric vector
#' @return the standard error of the mean for x as a numeric value
#'@export
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

#' Summarize cell level data to the spot level
#' 
#' Median summarize the cell level normalized values and the most biologically 
#' interpretable raw data at each spot, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param dt The datatable of cell level data to be summarized
#'  @return A datatable of spot-level, median summarized data with standard error values and 
#'  metadata
#' @export
summarizeToSpot <- function(dt, seNames=NULL){
  #Summarize cell data to medians of the spot parameters
  parameterNames<-grep(pattern="(Children|_CP_|_PA_)",x=names(dt),value=TRUE)
  
  #Remove spot-normalized or summarized parameters
  parameterNames <- grep("SpotNorm|Wedge|Sparse|OuterCell|Center|^Nuclei_PA_Gated_EduPositive$",parameterNames,value=TRUE,invert=TRUE)
  
  slDT<-dt[,lapply(.SD, numericMedian), by="Barcode,Well,Spot", .SDcols=parameterNames]
  #Use seNames to select the parameters that get SE values
  if(!is.null(seNames)){
    seNamesPattern<-paste(seNames,collapse="|")
    seNames <- grep(seNamesPattern,parameterNames,value=TRUE)
    slDTse <- dt[,lapply(.SD,se), by="Barcode,Well,Spot", .SDcols=seNames]
  } else{
    slDTse <- dt[,lapply(.SD,se), by="Barcode,Well,Spot"]
  }
  
  #Add _SE to the standard error column names
  setnames(slDTse, grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the spot and well metadata
  metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Drug|Endpoint|ECMp|MEP|Well_Ligand|ECM|ImageID|Barcode|^Well$|^PrintSpot$|^Spot$|Pin|Lx|[[:digit:]]{3}nm$|StainingSet)", x=colnames(dt), value=TRUE)
  setkey(dt,Barcode, Well,Spot)
  mDT <- dt[,metadataNames,keyby="Barcode,Well,Spot", with=FALSE]
  slDT <- mDT[slDT, mult="first"]
  #Merge in the standard err values
  setkey(slDT, Barcode, Well, Spot)
  setkey(slDTse, Barcode, Well, Spot)
  slDT <- slDT[slDTse]
  return(slDT)
}

#' QA the spot level data
#' 
#' Median summarize the cell level normalized values and the most biologically 
#' interpretable raw data at each spot, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param dt Datatable of spot level data
#' @param lthresh The threshold used in the loess model to define low cell count
#'  regions
#'  @return A datatable of spot-level, median summarized data with standard error values and 
#'  metadata
#' @export
QASpotData <- function(dt, lthresh = lthresh){
  #Add a count of replicates
  dt <- dt[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp"]
  #Add the loess model of the SpotCellCount on a per well basis
  dt <- dt[,Spot_PA_LoessSCC := loessModel(.SD, value="Spot_PA_SpotCellCount", span=.5), by="Barcode,Well"]
  #Add well level QA Scores to spot level data
  dt <- dt[,QAScore := calcQAScore(.SD, threshold=lthresh, maxNrSpot = max(dt$ArrayRow)*max(dt$ArrayColumn),value="Spot_PA_LoessSCC"),by="Barcode,Well"]
  return(dt)
}

#' Add the porportions of gated cells in the spot populations
#' @param dt A datatable 
#' @return nothing is returned, the proportions are added to the input datatable
#' @export
addSpotProportions <- function(dt){
  #Calculate DNA proportions based on cell cycle state
  dt <- dt[,Nuclei_PA_Cycle_DNA2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
  dt <- dt[,Nuclei_PA_Cycle_DNA4NProportion := 1-Nuclei_PA_Cycle_DNA2NProportion]
  
  #Create proportions and logits of mutlivariate gates
  if ("Cytoplasm_PA_Gated_KRTClass" %in% colnames(dt)){
    #Determine the class of each cell based on KRT5 and KRT19 class
    #0 double negative
    #1 KRT5+, KRT19-
    #2 KRT5-, KRT19+
    #3 KRT5+, KRT19+
    #Calculate gating proportions for EdU and KRT19
    dt <- dt[,Cytoplasm_PA_Gated_KRT5RT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
    dt <- dt[,Cytoplasm_PA_Gated_KRT5NegativeKRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
    dt <- dt[,Cells_PA_Gated_KRT5PositiveKRT19NegativeProportion := sum(Cytoplasm_PA_Gated_KRTClass==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
    dt <- dt[,Cells_PA_Gated_KRT5PositivedKRT19PositiveProportion := sum(Cytoplasm_PA_Gated_KRTClass==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
  }
  
  if ("Cells_PA_Gated_EdUKRT5Class" %in% colnames(dt)){
    #Determine the class of each cell based on KRT5 and EdU class
    #0 double negative
    #1 KRT5+, EdU-
    #2 KRT5-, EdU+
    #3 KRT5+, EdU+
    #Calculate gating proportions for EdU and KRT5
    dt <- dt[,Cells_PA_Gated_EdUKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==0)/length(ObjectNumber),by="Barcode,Well,Spot"]
    dt <- dt[,Cells_PA_Gated_EdUPositiveKRT5NegativeProportion := sum(Cells_PA_Gated_EdUKRT5Class==2)/length(ObjectNumber),by="Barcode,Well,Spot"]
    dt <- dt[,Cells_PA_Gated_EdUNegativeKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==1)/length(ObjectNumber),by="Barcode,Well,Spot"]
    dt <- dt[,Cells_PA_Gated_EdUKRT5PositiveProportion := sum(Cells_PA_Gated_EdUKRT5Class==3)/length(ObjectNumber),by="Barcode,Well,Spot"]
  }
}



# #' @export
# summarizeFBS <- function(dt){
#   #Summarize all by ligand and ECMp
#   #This will only change the FBS data
#   dtFBS <- dt[grepl("FBS",dt$Ligand),]
#   #Remove anything after FBS
#   dtFBS$Ligand <- gsub("FBS.*","FBS", dtFBS$Ligand)
#   #Create the new MEP name
#   dtFBS$MEP <- paste(dtFBS$ECMp,dtFBS$Ligand,sep="_")
#   if("StainingSet" %in% colnames(dtFBS)){
#     #Find the medians or the unique metadata values
#     dtFBSMedians <- dtFBS[, lapply(.SD, numericMedianUniqueMetadata), keyby = "Ligand,ECMp,StainingSet"]
#   } else {
#     #Find the medians or the unique metadata values
#     dtFBSMedians <- dtFBS[, lapply(.SD, numericMedianUniqueMetadata), keyby = "Ligand,ECMp"]
#   }
#   #If the FBS is from multiple plates its barcode has been reset to NA.
#   #replace this with the word Multiple
#   dtFBSMedians$Barcode[is.na(dtFBSMedians$Barcode)] <- "Multiple"
#   #Delete the FBS wells from the original dt
#   dt <- dt[!grepl("FBS",dt$Ligand),]
#   #Bind in the summarized FBS values
#   dt <- rbind(dt,dtFBSMedians)
#   return(dt)
# }

# #' Summarize spot level data to the MEP level
# #' 
# #' Median summarize the spot level normalized values the most biologically 
# #' interpretable raw data at each spot, calculate the standard errors and add
# #' SE columns for all median summarized data
# #' @param l3 The datatable of spot level data to be summarized
# #' @return A datatable of MEP level, median summarized data with standard error values and 
# #'  metadata
# #' @export
# createl4 <- function(l3, seNames=NULL){
#   #Add a count of replicates
#   l3 <- l3[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp"]
#   l4Names<-grep("Loess$|RUV|Norm|^Ligand$|^ECMp|Barcode|Spot_PA_SpotCellCount$|Spot_PA_ReplicateCount$", x=names(l3),value=TRUE)
#   #remove the _SE values
#   l4Names <- grep("_SE|NormMethod|AnnotID",l4Names, value = TRUE, invert = TRUE)
#   l4Keep<-l3[,l4Names,with=FALSE]
#   l4DT<-l4Keep[,lapply(.SD,numericMedian),keyby="Ligand,ECMp,Barcode"]
#   #Use seNames to select the parameters that get SE values
#   if(!is.null(seNames)){
#     seNamesPattern<-paste(seNames,collapse="|")
#     seNames <- grep(seNamesPattern,l4Names,value=TRUE)
#     l4DTse <- l4Keep[,lapply(.SD,se),keyby="Ligand,ECMp,Barcode", .SDcols=seNames]
#   } else{
#     l4DTse <- l4Keep[,lapply(.SD,se),keyby="Ligand,ECMp,Barcode"]
#   }
#   
#   #Add _SE to the standard error column names
#   setnames(l4DTse, grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(l4DTse), value = TRUE, invert = TRUE),"_SE"))
#   
#   l3Names <- grep("Barcode|Well|CellLine|Ligand|ECM|Endpoint488|Endpoint555|Endpoint647|EndpointDAPI|ECMp|MEP|Lx|PinDiameter", colnames(l3), value=TRUE)
#   #Merge back in the replicate metadata
#   mDT <- l3[,l3Names,keyby="Ligand,ECMp,Barcode", with=FALSE]
#   setkey(mDT,Ligand,ECMp,Barcode)
#   l4DT <- mDT[l4DT, mult="first"]
#   l4DT <- l4DTse[l4DT]
#   l4DT <- summarizeFBS(l4DT)
#   return(l4DT)
# }#End of createl4

#' Summarize spot level data to the MEP level
#' 
#' Median summarize the spot level normalized values the most biologically 
#' interpretable raw data at each spot, calculate the standard errors and add
#' SE columns for all median summarized data
#' @param dt The datatable of spot level data to be summarized
#' @return A datatable of MEP level, median summarized data with standard error values and 
#'  metadata
#' @export
preprocessLevel4 <- function(dt, seNames=NULL){
  #Add a count of replicates
  dt <- dt[,Spot_PA_ReplicateCount := .N,by="Ligand,ECMp,Drug,CellLine"]
  rawSignalNames <- grep("_SE",grep("Log2|Logit|_PA_|Intensity|AreaShape",colnames(dt), value=TRUE), value=TRUE, invert=TRUE)
  l4Signals<- dt[,lapply(.SD, numericMedian), by="Ligand,ECMp,Drug,CellLine", .SDcols=rawSignalNames]
  
  #Use seNames to select the parameters that get SE values
  if(!is.null(seNames)){
    seNamesPattern<-paste(seNames,collapse="|")
    seSignalNames <- grep(seNamesPattern,rawSignalNames,value=TRUE)
    l4Ses <- dt[,lapply(.SD,MEMA:::se),keyby="Ligand,ECMp,Drug,CellLine", .SDcols=seSignalNames]
  } else{
    l4Ses <- dt[,lapply(.SD,MEMA:::se),keyby="Ligand,ECMp,Drug,CellLine"]
  }
  
  #Add _SE to the standard error column names
  setnames(l4Ses, grep("Ligand|ECMp|Drug|CellLine",colnames(l4Ses), value = TRUE, invert = TRUE), paste0(grep("Ligand|ECMp|Drug|CellLine",colnames(l4Ses), value = TRUE, invert = TRUE),"_SE"))
  
  #Merge back in the replicate metadata
  metadataNames <- grep("_SE|Barcode|^BW$|ArrayRow|ArrayColumn|^Well$|^Spot$|^PrintSpot$|^Well_Ligand$|ImageID|QA_|ECMSet|^Row$|^Column$|^Block$|PlateRow|^QAScore$|^ID$",colnames(dt), value=TRUE,invert=TRUE) %>%
    setdiff(rawSignalNames)
  mdDT <- unique(dt[,metadataNames, with=FALSE])
  setkey(l4Signals,Ligand,ECMp,Drug,CellLine)
  setkey(l4Ses,Ligand,ECMp,Drug,CellLine)
  setkey(mdDT,Ligand,ECMp,Drug,CellLine)
  l4DT <- merge(mdDT,merge(l4Signals,l4Ses))
  return(l4DT)
}#End of preprocesslevel4
