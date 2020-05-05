#' Generic function for coercing flowFrame like objects to flowFrame
#'
#' @param x flowFrame like object
#' @param ... other parameters passed to specific methods
#' @export
#' @import methods
setGeneric("as.flowFrame", function(x, ...) {
  standardGeneric("as.flowFrame")
})

#' Coerces fcbFlowFrames to flowFrames by extracting the barcoded.ff
#'
#' @param fcbFlowFrame fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing,
#' clustering (at least one barcodes slot filled), and assigning
#' @return the barcoded flowframe from slot 1
#' @export
setMethod("as.flowFrame", "fcbFlowFrame", function(x, ...) {

  if (!any(class(x) == "fcbFlowFrame")) {
    stop("Input must be an object of class fcbFlowFrame")
  }
  x@barcoded.ff
})
