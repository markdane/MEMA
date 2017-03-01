
# #Parse names for ECM, Ligand and Concentration values
# #
# #\code{parseContents} returns a datatable with new columns for the ECM, ligand
# #and concentration at each spot.
# #
# #@param dt a datatable with a Name column to be parsed
# #@return The input datatable with new columns for the ECM, Ligand and
# #  Concentration values.
# #@section Usage: This function assumes a Name column exists and that its
# #  format is Concentration_ECM_Ligand or ligand_ECM or Col I_AF 488 or ECM. If Name only
# #  contains an ECM, a NULL string is assigned to the Ligand and a
# #  concentration value of 0 is assigned to Concentration. If there are 2 words
# #  separated by an underscore, the first is the Ligand and the second is the
# #  ECM. The concentration is assigned a value of 0. If there are 3 words,
# #  the first is assigned to the concentration.
# #
# #  In the special case of the AF 488 fiducial, the Ligand is assigned the value
# #  AF 488.
# #
# parseContents <- function(dt){
#   contentList <- lapply(strsplit(dt$Name,"_"),function(n){
#     if(length(n)==1){ECM<-n[[1]]
#     Ligand<-""
#     Concentration<-0
#     return(list(Ligand = Ligand, ECM = ECM,Concentration=Concentration))
# 
#     } else if(length(n)==2) {
#       Ligand<-n[[1]]
#       ECM<-n[[2]]
#       if(Ligand=="Col I"&ECM=="AF 488") {
#         Ligand<-"AF 488"
#         ECM<-"Col I"
#       }
#       Concentration <- 0
#       return(list(Ligand = Ligand, ECM = ECM,Concentration = Concentration))
# 
#     } else{
#       Concentration <- as.integer(sub("X","",x = n[[1]]))
#       Ligand        <- n[[2]]
#       ECM           <- n[[3]]
#       return(list(Ligand = Ligand, ECM = ECM,Concentration = Concentration))
#     }
#   }
#   )
#   contents <- data.table::rbindlist(contentList)
#   contents$ConcentrationRank <- match(contents$Concentration,sort(unique(contents$Concentration)))
#   dtA <- cbind(dt,contents)
# }
# 
# #Parse BufferWinners names for ECM, Ligand and Concentration values
# #
# #\code{parseBufferWinnerContents} returns a datatable with new columns based on the Name column.
# #
# #@param dt a datatable with a Name column to be parsed
# #@return The input datatable with new columns of Condition1, Condition2, Condition3 and Condition4.
# #@section Usage: This function assumes a Name column exists and uses an underscore as field separators.
# #
# parseBufferWinnerContents <- function(dt) {
#   #parse the content names in the gal file for Ligand, ECM and concentrations
#   contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
#     return(list(Condition1 = n[[1]], Condition2 = n[[2]],Condition3 = n[[3]],Condition4 = n[[4]]))
#   }
#   )
#   contents <- data.table::rbindlist(contentList)
#   dtA <- cbind(dt,contents)
#   return(dtA)
# }
# 
# #Parse 4wellvalidationContents for ECM, Glycerol, Triton, EDTA and Buffer values
# #
# #\code{parse4wellvalidationContents} returns a datatable with new columns based on the Name column.
# #
# #@param dt a datatable with a Name column to be parsed
# #@return The input datatable with new columns of ECM, Glycerol, Triton, EDTA and Buffer.
# #@section Usage: This functions assumes a Name column exists and uses an underscore as field separators.
# #
# parse4wellvalidationContents <- function(dt){
#   #parse the content names in the gal file for Ligand, ECM and concentrations
#   contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
#     return(list(ECM = n[[1]], Gycerol = n[[2]],Triton = n[[3]],EDTA = n[[4]],Buffer = n[[5]]))
#   }
#   )
#   contents <- data.table::rbindlist(contentList)
#   dtA <- cbind(dt,contents)
#   return(dtA)
# }

#Parse simulated data for ECM, Ligand and ConcentrationRank values
#
#\code{parseSimulatedContents} returns a datatable with new columns based on
#the Name column.
#
#@param dt a datatable with a Name column to be parsed
#@return The input datatable with new columns of ECM, Ligand and ConcentrationRank.
#@section Usage: This function assumes a Name column exists and uses an
#  underscore as the field separators. The format must be either
#  ECM or ECM_Ligand_ConcentrationRank. ConcentrationRank will be coereced into
#  an integer by removing any non-numeric values.
#
parseSimulatedContents <- function(dt){
  #parse the content names in the gal file for Ligand, ECM and concentrations
  contentList <- lapply(strsplit(dt$Name,"_"),function(n) {
    if (length(n)==2) stop("Invalid format in gal Name field. There should be 1 ECM or 1 ECM with 1 Ligand and a concentration value.")
    if(length(n)==1){ECM<-n[[1]]
    Ligand<-""
    ConcentrationRank<-0
    } else {
      ECM               <- n[[1]]
      Ligand            <- n[[2]]
      ConcentrationRank <- as.integer(gsub("[^[:digit:]]*","",n[[3]]))
    }

    return(list(Ligand = Ligand, ECM = ECM,ConcentrationRank = ConcentrationRank))
  }
  )
  contents <- data.table::rbindlist(contentList)
  dtA <- cbind(dt,contents)
}

# #' Coerce the All ECM contents to standard format
# #'
# #' \code{coerceAllECM}
# #'
# #'@param df A dataframe read from a gal file
# #'@return A dataframe with the Name column reformatted to match the Controlled Vocabulary standard.
# #'
# #'@export
# coerceAllECM <- function(df){
#   df$Name <- gsub("__","_",df$Name)
#   df$Name <- gsub("PBSCOL I","PBS_none_COL1_UNKNOWN",df$Name)
#   df$Name <- gsub("^PBS$","PBS_none",df$Name)
#   df$Name <- gsub("COL I$","COL1_UNKNOWN",df$Name)
#   df$Name <- gsub("^-$","blank_none",df$Name)
# 
#   contentMatrix <- limma::strsplit2(x = df$Name,split="_")
# 
#   df$Name <- apply(contentMatrix,1,function(x){
#     pUIDs <- paste0(x[seq(1,length(x),2)],"_",x[seq(2,length(x),2)])
#     pUIDs <- gsub("^_$","",pUIDs)
#     pUID <- paste0(pUIDs,collapse="-")
#     pUID <- gsub("[-]*$","",pUID)
#   }
#   )
#   return(df)
# }
# 
# #Parse controlled vocabulary data for ECM, Ligand and ConcentrationRank values
# #
# #\code{parseCVContents} returns a datatable with new columns based on the Name
# #column.
# #
# #@param dt A dataframe or datatable with a Name column to be parsed.
# #@return The input datatable with a new column of the ECM names pasted
# #  togethers with an underscore.
# #@section Usage: This functions assumes a Name column exists and uses a dash
# #  symbol - as a field separator between each protein and an underscore between
# #  a protein's name and its Uniprot ID.
# #
# parseCVContents <- function(dt){
#   dt$ECM <- lapply(dt$Name,function(x){
#     tmp<-strsplit(unlist(strsplit(x,split = "[-]")),split="_")
#     ECMNames<-lapply(tmp, function(x){
#       x[1]
#     })
#     paste(ECMNames,collapse="_")
#   })
# }
# 
# 
# # Extract Intensity Endpoints
# #
# extractIntsEndpts<-function(DT){
#   #Use the metadata in the columns Endpoint.xxx to extract the columns of
#   #mean intensity data and endpoint name
#   epColNames<-grep("(Endpoint)",colnames(DT),value = TRUE)
#   intColNames<-paste0("Mean.Intensity.Alexa.",sub("Endpoint.","",epColNames))
#   selectColNames<-c(epColNames,intColNames)
#   ieDT<-DT[,selectColNames,with=FALSE]
#   return(ieDT)
# }
# 
# # Assign endpoint names
# #
# assignEndPoints<-function(iedt,dt){
#   Endpoints<-unique(unlist(iedt[,grep("(Endpoint)",colnames(dt),value = TRUE),with=FALSE]))
#   #TODO Generalize this code to determine how many endpoints are in the staining set
#   #This code now assumes there are 3 stained endpoints and ignores DAPI
#   ei<-do.call("cbind",lapply(Endpoints,function(ep,ieDT){
#     #create a column of intensities for each endpoint based on the Endpoint column
#     endpointIntsPtrs<-which(ieDT[,1:3,with=FALSE]==ep,arr.ind = TRUE)
#     endpointIntsPtrs<-endpointIntsPtrs[order(endpointIntsPtrs[,1]),]
#     iValuesM<-as.matrix(ieDT[,4:6,with=FALSE])
#     epValues<-iValuesM[endpointIntsPtrs]
#   },ieDT=iedt))
#   colnames(ei)<-Endpoints
#   DT<-cbind(dt,ei)
# }

# # Summarize the cell level data to the well level
# #
# wellLevelData<-function(cd){
#   #browser()
#   #Create a datatable of the columns to be summariazed
#   intNames<-grep(pattern="(Intensity|Area|Well$|WellCellCount)",x=names(cd),value=TRUE)
#   #Get the metadata column names which are unconstrained
#   if("Spot" %in% names(cd)){
#     mdNames<-grep(pattern="(Intensity|Area|WellCellCount|Object.ID|Parent.Object.ID..MO|Elongation.Factor|Circularity.Factor|Perimeter|Max.Feret.Diameter|Center.X|X|Center.Y|Y|Position|WellIndex|Spot|Grid|Column|Row|ID|Name|ArrayRow|ArrayColumn)",x=names(cd),value=TRUE,invert=TRUE)
#   } else { mdNames<-grep(pattern="(Intensity|Area|WellCellCount|Object.ID|Parent.Object.ID..MO|Elongation.Factor|Circularity.Factor|Perimeter|Max.Feret.Diameter|Center.X|X|Center.Y|Y|Position|WellIndex)",x=names(cd),value=TRUE,invert=TRUE) }
#   keep<-cd[,intNames,with=FALSE]
#   keepMd<-cd[,mdNames,with=FALSE]
#   #Take the mean of each column stratified by well
#   wd<-keep[,lapply(.SD,mean),by=Well]
#   #Get the unique metadata value for each well
#   #The metadata of a spotted well is not unique
#   md<-keepMd[,lapply(.SD,unique),by=Well]
#   data.table::setkey(wd,Well)
#   data.table::setkey(md,Well)
#   all<-wd[md]
#   return(all)
# }

# Add the array row, column and index to a GAL file
#
addArrayPositionNoRotate<-function(df,gridsPerRow=4){
  #Handle dataframes that have the name block instead of grid
  blockToGrid<-FALSE
  if("Block" %in% colnames(df))
  {
    blockToGrid<-TRUE
    data.table::setnames(df,"Block","Grid")
  }
  #Return the dataframe in the correct space
  #These values are in the printer space
  rowsPerGrid<-max(df$Row)
  colsPerGrid<-max(df$Column)
  #Assign the grid row. These will run from 1:12 when using a 4x12 print head
  df$arrayGridRow<-ceiling(df$Grid/gridsPerRow)
  #Assign the ArrayRow which is in the imaging space
  df$ArrayRow<-(df$arrayGridRow-1)*rowsPerGrid+df$Row
  #Assign the ArrayColumn in imaging space (No rotation)
  df$ArrayColumn<-((df$Grid-1)%%(gridsPerRow))*colsPerGrid+df$Column
  #Order the array by ArrayRow then ArrayColumn
  df<-df[order(df$ArrayRow,df$ArrayColumn),]
  #Remove the arrayGridRow column used for calculations
  df<-subset(df,select=-c(arrayGridRow))
  #Assign a spot number in sequential order by row then column
  df$Spot<-1:(max(df$ArrayColumn)*max(df$ArrayRow))
  #Handle dataframes that have the name block instead of grid
  if(blockToGrid) data.table::setnames(df,"Grid","Block")
  return(df)
}


#' Wrapper function to log intensity values in a data.table
#' @export
logIntensities <- function(dt){
  intensityNames <- grep("Intensity",colnames(dt), value=TRUE, ignore.case = TRUE)
  dtLog <- dt[,lapply(.SD,log2),.SDcols=intensityNames]
  setnames(dtLog,colnames(dtLog),paste0(colnames(dtLog),"Log2"))
  return(cbind(dt,dtLog))
}

#' Utility function that returns log2 with lower non-infinite bound  
#' @export
boundedLog2 <- function(x){
  #use the smallest non zero value as the minimum
  xMin <- min(x, na.rm=TRUE)
  if (xMin==0) xMin <- unique(x[order(x)])[2]
  x[x==0]<- xMin
  xLog2 <- log2(x)
  return(xLog2)
}

#' Create a median normalized loess model of an array
#'
#'@param data A dataframe with ArrayRow, ArrayColumn and signal intensity columns
#'@param value The column name of the signal intensity column
#'@param span The span value passed to loess. Values between 0 and 1 determine the
#'proportion of the population to be included in the loess neighborhood.
#'@return a vector of median normalized loess values of the signal
#'@export
loessModel <- function(data, value, span){
  dataModel <- loess(as.formula(paste0(value," ~ ArrayRow+ArrayColumn")), data,span=span)
  dataPredicted <- predict(dataModel)
  predictedMedian <- median(dataPredicted, na.rm = TRUE)
  dataNormed <- dataPredicted/predictedMedian
}

#' Returns a character vector of alphanumeric well names
#'
#' \code{wellAN} is a convenience function that provides alphanumeric names
#' for well plates or arrays.
#' @param nrRows The number of rows in the plate or array.
#' @param nrCols The number of columns in the plate or array.
#' @return A character vector of the well names
#' @export
wellAN<-function(nrRows,nrCols){
  if(nrRows>702)stop("Too many rows to convert. Well alphanumerics are limited to 2 letters")
  Well=unlist(c(lapply(1:nrRows, function(x){
    paste0(paste0(LETTERS[(x-1)%/%26],LETTERS[(1+(x-1)%%26)]), lapply(1:nrCols,function(x){
      if(x<10) x<-paste0("0",x)
      return(as.character(x))
    }))
  })
  ))
  return(Well)
}

#' Change the periods in the column names to spaces
#'
#'\code{cleanColumnNames} substitues a space for any single period in the colmun names of a datatable. This is useful at the end of an anlaysis to have human readable display names.
#' @param dt a datatable
#' @return the same datable with spaces in place of periods.
#' @export
cleanColumnNames<-function(dt){
  data.table::setnames(dt,gsub("[.]"," ",make.names(colnames(dt))))
  return(dt)
}

#' Get barcodes from the Synapse Plate Tracker files
#' @param studyName Character string of the study name
#' @return Barcodes of the plates in the study
#' @export
getBarcodes <- function(studyName){
  library(synapseClient)
  
  synapseLogin()
  barcodes <- synGet("syn8313413") %>%
    synapseClient::getFileLocation() %>%
    data.table::fread(check.names = TRUE) %>%
    dplyr::filter(Study.Name == studyName) %>%
    dplyr::select(Plate.IDs) %>%
    stringr::str_split(",") %>%
    unlist()
  return(barcodes)
}