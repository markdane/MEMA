#Plotting functions in MEMA reports

#' Plot a heatmap with ligand rows and ECM protein columns
#' 
#' @param dt Datatable with Ligand, ECMp and a column named by the fill value
#' @param hclustfun Optional clustering function to overide heatmap.2 default
#' @param fill Name of column with numeric values for coloring the heatmap
#' @param titleExpression Title for heatmap
#' @param limits Length 2 integer or numeric vector used to squish the colors of extreme values
#' @param xAxisSize Optional scale factor for x axis text
#' @param yAxisSize Optional scale factor for y axis text
#' @return none called for its side effect of plotting the heatmap
#' @export
plotLEHmap <- function(dt, hclustfun = NULL, fill, titleExpression, limits, xAxisSize=.6, yAxisSize=.5){
  if(is.null(hclustfun)){
    p <- ggplot(dt, aes_string(x="Ligand", y="ECMp", fill = fill))+
      geom_tile()+
      scale_fill_gradient(low="white",high="red",oob = scales::squish,
                          limits=limits)+
      ggtitle(titleExpression)+
      xlab("")+ylab("")+
      guides(fill=FALSE)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(xAxisSize)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(yAxisSize)), plot.title = element_text(size = rel(1)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.7)),panel.grid.major = element_line(linetype = 'blank'),panel.background = element_rect(fill = "dimgray"))
    suppressWarnings(print(p))
  } else {
    hmcols<-colorRampPalette(c("blue","white","red"))(16)
    #Cast to get ligands into columns
    df <- data.table::dcast(data.frame(dt[,c("ECMp","Ligand",fill), with=FALSE]),ECMp~Ligand, value.var = fill, fun.aggregate=median)
    
    rownames(df) <- df$ECMp
    df <- df[,!grepl("ECMp",names(df))]
    par(cex.main=0.7)
    try(heatmap.2(as.matrix(dfZoom(df, limits[1], limits[2])), hclustfun=get(hclustfun), col=hmcols, trace="none", cexRow=.6, cexCol=.6, main=titleExpression, xlab="Ligands", ylab="ECM Proteins"), TRUE)
    
  }
  
}