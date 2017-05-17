## ----setup---------------------------------------------------------------
library(MEMA)
library(synapseClient)
library(parallel)
library(stringr)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)


## ----login---------------------------------------------------------------
suppressMessages(synapseLogin())

## ----getMetadataFromSynapse----------------------------------------------

barcode <- "LI8X00651"
metadataq <- sprintf("select id from syn7494072 WHERE DataType='Metadata' AND Barcode='%s'",barcode)
metadataTable <- synTableQuery(metadataq, verbose=FALSE)@values
metadataFiles <- lapply(metadataTable$id, synGet)
metadataFiles <- list(annotMetadata=getFileLocation(metadataFiles[[1]]))

## ----readAndProcessMetadata----------------------------------------------

metadata <- getMetadata(metadataFiles, useAnnotMetadata = TRUE)


## ----showMetadata--------------------------------------------------------
head(metadata[,.(Barcode, Well, Spot, Ligand, ECMp, CellLine, StainingSet)])

## ----getRawDataFromSynapse-----------------------------------------------

q <- sprintf("select id,name,Barcode,Level,Well,StainingSet,Location from syn9838977 WHERE Level='0' AND Barcode='%s'", barcode)
rawFiles <- synTableQuery(q)
dataBWInfo <- rawFiles@values

# Download raw files, or get from cache if already downloaded
res <- lapply(dataBWInfo$id, synGet)

# Get the path on disk to each file
cellDataFilePaths <- unlist(lapply(res, getFileLocation))
dataBWInfo$Path <- cellDataFilePaths

dtL <- getCPData(dataBWInfo = dataBWInfo, verbose=FALSE)


## ----headRawData---------------------------------------------------------
head(dtL[[1]][,.(Nuclei_CP_Intensity_MedianIntensity_Dapi,Nuclei_CP_AreaShape_Area)])

## ----cleanAndProcess-----------------------------------------------------
dtL <- lapply(dtL, function(dt){
  filterObjects(dt,nuclearAreaThresh = 50, nuclearAreaHiThresh = 4000)})

## ----addDerivedFeatures--------------------------------------------------
# Add local XY and polar coordinates
dtL <- lapply(dtL, addPolarCoords)

# Add spot level normalizations for median intensities
dtL <- lapply(dtL,spotNormIntensities)

# Add adjacency values
dtL <- lapply(dtL, calcAdjacency)

#Clean up legacy issues in column names and some values
dtL <- lapply(dtL, cleanLegacyIssues)

## ----mergeDataAndMetadata------------------------------------------------

cDT <- merge(rbindlist(dtL),metadata,by=c("Barcode","Well","Spot"))

## ----gate----------------------------------------------------------------
# Gate cells where possible
cDT <- gateCells(cDT)

## ----DNAHistogram, fig.width=6, fig.height=4-----------------------------
dt <- cDT[grepl("COL1",cDT$ECMp),]
dt <- dt[,TotalDNANormed := Nuclei_CP_Intensity_IntegratedIntensity_Dapi/median(Nuclei_CP_Intensity_IntegratedIntensity_Dapi[Nuclei_PA_Cycle_State==1]), by="Barcode,Well"]
dt <- dt[,DNAThresh := min(TotalDNANormed[Nuclei_PA_Cycle_State==2]), by="Barcode,Well"]
dt <- dt[,Well_Ligand :=  paste(Well,Ligand, sep="\n")]
binwidth=.04

  p <- ggplot(dt, aes(x=TotalDNANormed))+
    geom_histogram(binwidth = binwidth)+
    coord_cartesian(x=c(0,4))+
    geom_vline(data = dt, aes(xintercept = DNAThresh), colour = "blue")+
    facet_wrap(~Well_Ligand, nrow=2) +
    ylab("Count")+xlab("Total DAPI Intensity Normalized to 2n")+ggtitle(paste("Total DAPI Intensity:",barcode))+
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=rel(.5)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(.8)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.3)),strip.text.x = element_text(size = rel(.4)))
  suppressWarnings(print(p))



## ----showLevel1----------------------------------------------------------
head(cDT[,.(ObjectNumber, Barcode, Well, Spot, MEP, Nuclei_CP_Intensity_IntegratedIntensity_Dapi)])


## ----preprocessLevel2----------------------------------------------------
cDT <- cDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]


## ----addProportions------------------------------------------------------
#Add proportions for signals with multivariate gating and non-conforming gate values
addSpotProportions(cDT)

#Calculate proportions for binary gated signals
gatedSignals <- grep("Proportion", grep("Positive|High",colnames(cDT), value=TRUE), value=TRUE, invert=TRUE)
if(length(gatedSignals)>0){
  proportions <- cDT[,lapply(.SD, calcProportion),by="Barcode,Well,Spot", .SDcols=gatedSignals]
  setnames(proportions,
           grep("Gated",colnames(proportions),value=TRUE),
           paste0(grep("Gated",colnames(proportions),value=TRUE),"Proportion"))
}

## ----summarizeToSpot-----------------------------------------------------
#Define the regex patterns for the signals that get standard error values.
seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin")

#median summarize the rest of the signals to the spot level
signals <- summarizeToSpot(cDT, seNames)

if(exists("proportions")) {
  spotDT <- merge(signals,proportions)
} else {
  spotDT <- signals
}

#Set a threshold for the loess well level QA Scores
if(sum(c("ArrayRow","ArrayColumn") %in% colnames(spotDT))==2) spotDT <- QASpotData(spotDT, lthresh = .6)


## ----showLevel2----------------------------------------------------------
head(spotDT[,.(Barcode, Well, Spot, MEP, Spot_PA_SpotCellCount)])


## ----pseudoimaesQAScores, fig.width=3.7,fig.height=4---------------------
#Set a threshold for the loess well level QA Scores
spotDT <- QASpotData(spotDT, lthresh = .6)
dt <- spotDT[,Well_Ligand :=  paste(Well,Ligand, sep="\n")]

plotSCCHeatmapsQAHistograms(dt, barcode,lthresh = 0.6)


## ------------------------------------------------------------------------
sessionInfo()

