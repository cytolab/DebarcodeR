#' Defines populations on barcoded datasets
#'
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing
#' and clustering (at least one barcodes slot filled)
#' @param channel The name (string) of the channel that has been corrected and clustered
#' @param likelihoodcut numeric, a likelihood cutoff for discarding unlikely cells, less than 1/k as likely
#' as the most likely cell from that population
#' @param ambiguitycut numeric from 0 to 1, threshhold below which to discard ambigious cells, eg: 0.02,
#'  discards cells with more than 2% chance of originating from another population

#' @return a fcbFlowFrame object with a barcode slot filled with deskewing, clustering, cell assignment as
#' a vector of integers from 0:ncol(probs), cells assigned a classification of 0 remained unassigned,
#' otherwise number corresponds to the barcoding level assignment of that cell
#' @export

assign_fcbFlowFrameOld <- function(fcbFlowFrame,
                                channel,
                                likelihoodcut = 8 ,
                                ambiguitycut = 0.02){

  if (!any(class(fcbFlowFrame) == "fcbFlowFrame")) {
    stop("Input must be an object of class fcbFlowFrame")
  }

  if (length(fcbFlowFrame@barcodes) == 0) {
    stop(
      "Input must have channels in the barcodes slot that have been run through deskew_fcbFlowFrame"
    )
  }

  if (length(fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]]) == 1) {
    stop(
      "Input must have channels in the barcodes slot that have been run through cluster_fcbFlowFrame"
    )
  }

  probs =  fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][["clustering"]][["probabilities"]]


  row.max <-  apply(probs, 1, sum)
  probs.norm.row <- sweep(probs, 1, row.max, FUN = "/")

  col.max <-  apply(probs, 2, max)
  probs.norm.col <- sweep(probs, 2, col.max, FUN = "/")

  classif <- rep(0, nrow(fcbFlowFrame@barcoded.ff))

  if (ncol(probs) > 1) { # if assigning more than one level
    classif <- as.numeric(apply(probs.norm.row, 1, which.max))
  } else {# if assigning only one level
    classif <- rep(1, nrow(probs))
  }

  classif[which(is.na(classif))] <- 0 #catches the few cells with 0 probability of belonging to any pop
  likely <- probs.norm.col > 1/likelihoodcut

  if (ncol(probs) > 1) { # if assigning more than one level
    likely.sum <- apply(likely, 1, sum) #converts logical to numeric
    print(paste0(round(sum(apply(likely, 1, any))/nrow(likely)*100, 3), "% above likelihood cutoff"))
    classif[which(likely.sum < 1)] <- 0
    non.ambigious <- apply(probs.norm.row, 1, max) > (1 - ambiguitycut)
    classif[which(!non.ambigious)] <- 0
  } else {
    likely.sum <- as.numeric(likely)
    classif[which(likely.sum != 1)] <- 0
  }

  fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][[3]] <- list(values = classif,
                                                                                       ambiguity = ambiguitycut,
                                                                                       likelihood = likelihoodcut)
  names(fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]])[3] <-
    "assignment"

  return(fcbFlowFrame)
}

#' Ambiguity cutoff
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing
#' and clustering (at least one barcodes slot filled)
#' @return the prbability matrix normalized by row
#' @export

calculate.ambiguity <- function(fcbFlowFrame,channel)
{
  if (!any(class(fcbFlowFrame) == "fcbFlowFrame")) {
    stop("Input must be an object of class fcbFlowFrame")
  }

  if (length(fcbFlowFrame@barcodes) == 0) {
    stop(
      "Input must have channels in the barcodes slot that have been run through deskew_fcbFlowFrame"
    )
  }

  if (length(fcbFlowFrame@barcodes[[1]]) == 1) {
    stop(
      "Input must have channels in the barcodes slot that have been run through cluster_fcbFlowFrame"
    )
  }

  probs =  fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][["clustering"]][["probabilities"]]

  row.max <-  apply(probs, 1, sum)
probs.norm.row <- sweep(probs, 1, row.max, FUN = "/")

return(probs.norm.row)}

#' Likelihood cutoff
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing
#' and clustering (at least one barcodes slot filled)
#' @return the prbability matrix normalized by column
#' @export

calculate.likelihood <- function(fcbFlowFrame,channel){
  if (!any(class(fcbFlowFrame) == "fcbFlowFrame"))
    {
    stop("Input must be an object of class fcbFlowFrame")
  }

if (length(fcbFlowFrame@barcodes) == 0) {
  stop(
    "Input must have channels in the barcodes slot that have been run through deskew_fcbFlowFrame"
  )
}

if (length(fcbFlowFrame@barcodes[[1]]) == 1) {
  stop(
    "Input must have channels in the barcodes slot that have been run through cluster_fcbFlowFrame"
  )
}

  probs =  fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][["clustering"]][["probabilities"]]

  col.max <-  apply(probs, 2, max)
probs.norm.col <- sweep(probs, 2, col.max, FUN = "/")
return(probs.norm.col)}

