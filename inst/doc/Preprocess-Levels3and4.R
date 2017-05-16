## ----setup---------------------------------------------------------------
library(MEMA)
library(RUVnormalize)
library(ruv)
suppressPackageStartupMessages(library(synapseClient))
library(parallel)
library(stringr)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
suppressPackageStartupMessages(library(plotly))
library(d3heatmap)
library(RColorBrewer)


## ----login---------------------------------------------------------------
suppressMessages(synapseLogin())

## ----getLevel2DataFromSynapse--------------------------------------------
studyName <- "HMEC240L_SS4"
q <- sprintf("select id from syn7494072 WHERE Level='2' AND Study='%s'", studyName)
level2File <- synTableQuery(q)
dataSynID <- level2File@values$id

# Download level 2 file or get from cache if already downloaded
res <- lapply(dataSynID, synGet)

# Get the path on disk to the file
dataFilePath <- unlist(lapply(res, getFileLocation))

#Read level 2 data
slDT <- getSpotLevelData(dataFilePath)

## ----normalizeLevel2Data-------------------------------------------------

signalsMinimalMetadata <- grep("_SE",grep("_SpotCell|EdU|QA|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^Drug1Conc$|^ArrayRow$|^ArrayColumn$|^CellLine$",colnames(slDT), value=TRUE), value=TRUE, invert=TRUE)

#Define the maximum number of factors to be removed by the RUV normalization
k <- 256
#RUVLoess normalize all signals
nDT <- normRUVLoessResiduals(slDT[,signalsMinimalMetadata, with = FALSE], k)
nDT$NormMethod <- "RUVLoessResiduals"
slDT$k <- k
slDT <- merge(slDT, nDT, by = c("BW","PrintSpot"))

## ----QASpotlevelData-----------------------------------------------------
#Add QA flags to the data
slDT <- QASpotLevelData(slDT, lowSpotCellCountThreshold=5,
                        lowRegionCellCountThreshold = 0.4,
                        lowWellQAThreshold = .7)

## ----showLevel3----------------------------------------------------------
head(slDT[,.(Barcode, Well, Spot, MEP, Spot_PA_SpotCellCount, Spot_PA_SpotCellCountNorm)])


## ----createLevel4--------------------------------------------------------
#Summarize to the MEP_Drug level
mepDT <- preprocessLevel4(slDT,seNames=c("DNA2N","SpotCellCount","EdU","MitoTracker","KRT","Lineage","Fibrillarin"))

## ----QALevel4------------------------------------------------------------
#Add in the barcodes for each MEP_Drug
mepDT <- addBarcodes(dt3 = slDT, dt4 = mepDT)
# Add a QA flag for spots with few replicates
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < 3

## ----showLevel4----------------------------------------------------------
head(mepDT[,.(Barcode, MEP, Spot_PA_SpotCellCount, Spot_PA_SpotCellCountNorm)])

## ----SCCBoxplots, fig.width=8, fig.height=5------------------------------

dt <- mepDT[,list(Barcode,Ligand,Spot_PA_SpotCellCount, Spot_PA_SpotCellCountNorm)]
dt <- melt(dt,id.vars=c("Barcode","Ligand"),measure.vars=c("Spot_PA_SpotCellCount","Spot_PA_SpotCellCountNorm"), variable.name = "Processed",value.name = "SpotCellCount", variable.factor = FALSE)
dt$Processed[!grepl("Norm",dt$Processed)] <- "Raw"
dt$Processed[grepl("Norm",dt$Processed)] <- "Normalized"
dt$Processed <- factor(dt$Processed,levels = c("Raw","Normalized"))

p <- ggplot(dt, aes(x = reorder(Ligand, SpotCellCount, FUN=median), y = SpotCellCount, colour = Barcode))+geom_boxplot()+
  ggtitle(paste("Raw and Normalized Spot Cell Count by Ligand"))+
  xlab("")+ylab("")+
  geom_hline(yintercept = median(dt$SpotCellCount), colour="blue", alpha=.5)+
  facet_wrap(~Processed, ncol=1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(.8)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.5)),legend.title=element_text(size = rel(.8)))
p


## ----MEPSbySpotCellCount, fig.width=8, fig.height=5----------------------

dt <- mepDT

p <- ggplot(dt, aes(x =reorder(MEP, Spot_PA_SpotCellCountNorm), y = Spot_PA_SpotCellCountNorm))+
  geom_errorbar(aes(ymin=Spot_PA_SpotCellCountNorm-Spot_PA_SpotCellCountNorm_SE, ymax=Spot_PA_SpotCellCountNorm+Spot_PA_SpotCellCountNorm_SE), width=.01, colour="black") +
  xlab("MEP")+ylab("Normalized Spot Cell Count")+
  theme( axis.text.x=element_blank(), axis.ticks=element_blank(),  panel.grid.major = element_blank())+
  ggtitle("MEPs Ordered by Normalized Spot Cell Count")

p <- p + geom_point(aes(y=Spot_PA_SpotCellCountNorm),colour = "darkblue", alpha = .5)

ggplotly(p)

## ----SCCHeatmapFull, fig.width=8, fig.height=5---------------------------
#Cast to get ligands into columns
df <- dcast(data.frame(mepDT[,list(ECMp,Ligand,Spot_PA_SpotCellCountNorm,Barcode)]),ECMp~Ligand, value.var = "Spot_PA_SpotCellCountNorm",fun=median)

rownames(df) <- df$ECMp
df <- df[,!grepl("ECMp",names(df))]

hmcols<-colorRampPalette(c("blue","white","red"))(16)
try(d3heatmap(dfZoom(df, .05, .95), colors=hmcols, xaxis_font_size="6pt", yaxis_font_size="5pt"), TRUE)


## ----edUNormLigandBoxplots, fig.width=8, fig.height=5--------------------
p <- ggplot(mepDT, aes(x=Ligand, y=Nuclei_PA_Gated_EdUPositiveProportionNorm))+
  geom_boxplot(outlier.colour = NA, fill=NA)+
  geom_jitter(size=rel(.4))+
  #coord_cartesian(ylim = c(0,.5))+
  guides(colour=FALSE)+
  xlab("Ligand")+ylab("Normalized EdU+")+
  ggtitle("MEP EdU+ Response by Ligand")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(.7)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(1)),legend.title=element_text(size = rel(1)))

print(p)

## ------------------------------------------------------------------------
sessionInfo()

