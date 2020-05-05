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
#' @importClassesFrom flowCore flowFrame
#' @export
setClass("fcbFlowSet",
         representation = representation(fcbFrames = "list",
                                         uptake.ff = "flowFrame")
)

## constructor
#' @export

fcbFlowSet <- function(fcbFrames, uptake = NULL) {
  if (!any(class(fcbFrames) == "list")) {
    #needs to actually check that all the items are actually fcbFlowFrames
    stop("fcbFlowFrames must be a list of fcbFlowFrames")
  }

  if (!any(c(class(uptake.ff) == "flowFrame", is.null(uptake.ff)))) {
    stop("uptake.ff must be an object of class FlowFrame")
  }

  if (is.null(uptake.ff)) {
    warning("No uptake flowFrame supplied, will check each fcbFlowFrame for uptake control")
    #do something here
  }

  return(new("fcbFlowSet", fcbFrames = fcbFrames, uptake = uptake))
}
