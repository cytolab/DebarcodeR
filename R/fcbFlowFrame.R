#' fcbFlowFrame
#' -----------------------------------------------------------------------------
#' A containiner for barcoded flow cytometry data, with slots for the barcoded
#' flowFrame, a single level 'uptake control', and a slot to contain the results
#' of the debarcoding functions contained within debarcoder.
#'
#' @name fcbFlowFrame-class
#' @slot barcoded.ff {Object of class\code{flowFrame} containing the barcoded
#' data, approriately compensated, transformed, and gated}
#' @slot uptake.ff {object of class \code{flowFrame} containing cells which
#' have been barcoded with a single level of each barcoding channel as well
#' as stained with the approriate barcoding controls}
#' @slot barcodes {a list of barcodes, each one named for the channel from which
#' it was derived}
#' @slot platemap {a platemap for conditions per barcode level}
#' @import methods
#' @importClassesFrom flowCore flowFrame
#' @export
fcbFlowFrame <- setClass("fcbFlowFrame",
          contains = "flowFrame",
          slots = c(barcodes = "list")
         )

#' @export
#' @importClassesFrom flowCore flowFrame
fcbFlowFrame <- function(x,  barcodes = list()) {
  if (!any(class(x) == "flowFrame")) {
    stop("x must be an object of class FlowFrame")
    #should maybe leave the possibility of attempting to coerce to flowFrame
  }
#  print(31)
  as(x, "fcbFlowFrame")
}

