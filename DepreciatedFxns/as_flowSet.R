#' Generic function for coercing flowSet like objects to flowSets
#'
#' @param x flowSet like object
#' @param ... other parameters passed to specific methods
#' @export
#' @import methods

setGeneric("as.flowSet")
#' Generic function for coercing flowSet like objects to flowSets
#'
#' @param x ncdfFlowSet
#' @param ... other parameters passed to specific methods
#' @export
#' @importFrom ncdfFlow as.flowSet
setMethod("as.flowSet", "ncdfFlowSet", function(x, ...)  {
  ncdfFlow::as.flowSet(x, ...)
})

#' Writes FCS files for each well/level and returns debarcoded flowSet
#'
#' @param fcbFlowFrame fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing,
#' clustering (at least one barcodes slot filled), and assigning
#' @param simplify if TRUE, return only the debarcoded flowSet; if FALSE, return
#' debarcoded flowSet as slot in fcbFlowFrame (default is TRUE)
#' @return a matrix of probabilities, with ncol = levels, and nrow = legnth(vec).
#' @export
#' @import plyr tidyverse


setMethod("as.flowSet", "fcbFlowFrame", function(x, ...) {

  if (!any(class(x) == "fcbFlowFrame")) {
    stop("Input must be an object of class fcbFlowFrame")
  }

  if (length(x@barcodes) == 0) {
    stop(
      "Input must have channels in the barcodes slot that have been run through deskew_fcbFlowFrame"
    )
  }

  if (length(x@barcodes[[1]]) == 1) {
    stop(
      "Input must have channels in the barcodes slot that have been run through cluster_fcbFlowFrame"
    )
  }

  if (length(x@barcodes[[1]]) == 2) {
    stop(
      "Input must have channels in the barcodes slot that have been run through assign_fcbFlowFrame"
    )
  }

  barcoded.data = as.data.frame(x@barcoded.ff@exprs)
  for (i in 1:length(x@barcodes)) {
    barcoded.data[,ncol(barcoded.data) + 1] <- x@barcodes[[i]]$assignment$values
    colnames(barcoded.data)[ncol(barcoded.data)] <- paste0("levels_",names(x@barcodes)[i])}

  for (s in 1:length(x@barcodes)) {
    barcoded.data = subset(barcoded.data, barcoded.data[,c(ncol(barcoded.data) - s + 1)] != 0)}

  barcode.test <- barcoded.data %>%
    select(contains("levels")) %>%
    unite("well", remove = F)

  barcoded.data$well = (barcode.test$well)

  debarcoded.data <- barcoded.data %>%
    dplyr::select(-contains("level"), -contains("well")) %>%
    split(barcoded.data$well)

  ff_list <- vector("list",length(debarcoded.data))
  for (f in 1:length(ff_list)) {
    ff.i <- x@barcoded.ff
    exprs(ff.i) <- as.matrix(as.data.frame(debarcoded.data[[f]]))
    ff_list[[f]] <- ff.i
    names(ff_list)[[f]] <- names(debarcoded.data)[f]
  }

  dir.create(paste0("debarcoded_", Sys.Date()))
  orig.name <- sub('.fcs', "", sub("data/", "", ff_list[[1]]@description$FILENAME))


  debarcoded.fs = flowSet(ff_list)

  for (r in 1:length(ff_list)) {
    print(r)
    well.i <- names(ff_list)[[r]]
    ff <- flowSet(ff_list[[r]])
    flowframe.to.export = ff[[1]]
    write.FCS(flowframe.to.export,
              filename =  paste0("debarcoded_", Sys.Date(), "/", orig.name,"_",well.i, ".fcs"))
  }
    return(debarcoded.fs)

})
