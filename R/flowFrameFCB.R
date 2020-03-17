#' flowFrameFCB
#' -----------------------------------------------------------------------------
#' A containiner for barcoded flow cytometry data, with slots for the barcoded
#' flowFrame, a single level 'uptake control', and a slot to contain the results
#' of the debarcoding functions contained within debarcoder.
#'
#' @name flowFrameFCB-class
#' @slot barcoded.ff {Object of class\code{flowFrame} containing the barcoded
#' data, approriately compensated, transformed, and gated}
#' @slot uptake.ff {object of class \code{flowFrame} containing cells which
#' have been barcoded with a single level of each barcoding channel as well
#' as stained with the approriate barcoding controls}
#' @slot barcodes {a list of barcodes, each one named for the channel from which
#' it was derived}
#' @slot platemap {a platemap for conditions per barcode level}
#' @importClassesFrom flowCore flowFrame
#' @export
setClass("flowFrameFCB",
         representation = representation(barcoded.ff = "flowFrame",
                                         uptake.ff = "flowFrame",
                                         barcodes = "list",
                                         platemap = "list")
         )

## constructor
#' @export
flowFrameFCB <- function(barcoded.ff, uptake.ff = NULL) {
  if(!any(class(barcoded.ff) == "flowFrame")){
    stop("barcoded.ff must be an object of class FlowFrame")
  }

  if(!any(c(class(uptake.ff) == "flowFrame", is.null(uptake.ff)))){
    stop("uptake.ff must be an object of class FlowFrame")
  }

  if(is.null(uptake.ff)) {
    warning("No uptake flowFrame supplied, will use barcoded sample as uptake control")
    uptake.ff <- barcoded.ff
  }

  return(new("flowFrameFCB",barcoded.ff = barcoded.ff,
                            uptake.ff = uptake.ff,
                            barcodes = list(),
                            platemap = list()))
}
