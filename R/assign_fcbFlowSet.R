#' Defines populations on barcoded datasets
#'
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing
#' and clustering (at least one barcodes slot filled)
#' @param channel The name (string) of the channel that has been corrected and clustered
#' @param likelihoodcut numeric, a likelihood cutoff for discarding unlikely cells, less than 1/k as likely
#' as the most likely cell from that population
#' @param ambiguitycut numeric from 0 to 1, threshhold below which to discard ambigious cells, eg: 0.02,
#'  discards cells with more than 2\% chance of originating from another population

#' @return a fcbFlowFrame object with a barcode slot filled with deskewing, clustering, cell assignment as
#' a vector of integers from 0:ncol(probs), cells assigned a classification of 0 remained unassigned,
#' otherwise number corresponds to the barcoding level assignment of that cell
#' @export
assign_fcbFlowSet <- function(fcbFlowSet,
                                channel,
                                likelihoodcut = 8 ,
                                ambiguitycut = 0.02){

  if (class(fcbFlowSet) != "fcbFlowSet") {
    stop("Input must be a fcbFlowSet")
  }

  fcbFlowSet.assigned <- fsApply(fcbFlowSet, assign_fcbFlowFrame,
                                 channel = channel,
                                 likelihoodcut = likelihoodcut,
                                 ambiguitycut = ambiguitycut)

  return(fcbFlowSet(fcbFlowSet.assigned))
}
