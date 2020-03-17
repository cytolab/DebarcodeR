#' Writes FCS files for each well/level and returns debarcoded flowSet
#'
#' @param flowFrameFCB flowFrameFCB object with barcoded flowframe and uptake flowframe post deskewing,
#' clustering (at least one barcodes slot filled), and assigning
#' @param simplify if TRUE, return only the debarcoded flowSet; if FALSE, return
#' debarcoded flowSet as slot in flowFrameFCB (default is TRUE)
#' @return a matrix of probabilities, with ncol = levels, and nrow = legnth(vec).
#' @export
#' @import plyr tidyverse

as.flowSet <- function(flowFrameFCB) {

  barcoded.data = as.data.frame(flowFrameFCB@barcoded.ff@exprs)
  for (i in 1:length(flowFrameFCB@barcodes)){
    barcoded.data[,ncol(barcoded.data)+1]<- flowFrameFCB@barcodes[[i]]$assignment$values
    colnames(barcoded.data)[ncol(barcoded.data)]<-paste0("levels_",names(flowFrameFCB@barcodes)[i])}

  for (s in 1:length(flowFrameFCB@barcodes)){
    barcoded.data = subset(barcoded.data, barcoded.data[,c(ncol(barcoded.data)-s+1)] != 0)}

  barcode.test <- barcoded.data %>%
    select(contains("levels")) %>%
    unite("well", remove = F)

  barcoded.data$well = (barcode.test$well)

  debarcoded.data <- barcoded.data %>%
    dplyr::select(-contains("level"), -contains("well")) %>%
    split(barcoded.data$well)

  ff_list <- vector("list",length(debarcoded.data))
  for(f in 1:length(ff_list)) {
    ff.i <- flowFrameFCB@barcoded.ff
    exprs(ff.i) <- as.matrix(as.data.frame(debarcoded.data[[f]]))
    ff_list[[f]]<- ff.i
    names(ff_list)[[f]] <- names(debarcoded.data)[f]
  }

  dir.create(paste0("debarcoded_", Sys.Date()))
  orig.name <-sub('.fcs', "", sub("data/", "", ff_list[[1]]@description$FILENAME))


  debarcoded.fs = flowSet(ff_list)

  for(r in 1:length(ff_list)){
    print(r)
    well.i <- names(ff_list)[[r]]
    ff<- flowSet(ff_list[[r]])
    flowframe.to.export = ff[[1]]
    write.FCS(flowframe.to.export,
              filename =  paste0("debarcoded_", Sys.Date(), "/", orig.name,"_",well.i, ".fcs"))
  }
    return(debarcoded.fs)

}
