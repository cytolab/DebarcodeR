#' Defines populations on barcoded datasets
#'
#'This function allows you to calculate the probability of a cell originating from a given population using
#'either gaussian mixture modeling or jenks natural breaks classification
#'
#' @param fcbFlowSet a fcbFlowSet object with barcoded flowframe and uptake flowframe post deskewing (at least one barcodes slot filled)
#' @param channel The name (string) of the channel to be clustered
#' @param ret.model Option to retain the model for deskewing
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI
#' @param levels integer, the number of barcoding intensities present in the vector
#' @param opt string, either "mixture" (default) for gaussian mixture modeling, or "fisher" for fisher-jenks natural breaks optimization
#' @param dist string in c("Normal, Skew.normal, Tdist"), passed to mixsmsn
#' @param subsample Integer, number of cells to subsample, defaults to 10,000
#' @param trim numberic between 0, 1; used to trim the upper and lower extremes to exlcude outliers (eg. trim = 0.01 exludes most extreme 1\% of data)
#'
#' @return a fcbFlowFrame with deskewed barcodes slot and clustering slot with a matrix of probabilities, with ncol = levels, and nrow = legnth(vec).
#' If gaussian mixture modeling is used the probailities correspond to the probability
#' of the cell originaiting that level under the distrubtion specified by the mixture model
#' If jenks natural breaks optimization is used, the probability is estimated empirically based on a histogram
#'
#' @seealso \code{\link{deskew_fcbFlowFrame}}
#' @export
#' @import classInt mixsmsn sn

# aspirational --> v2?
# find sd and mean of uptake - use for probabilities
# bounds as 5th and 95th, separation of each level --> find probability
# uptake for two channels - use as model (Prior)
# READ smsn.mix paper and function

cluster_fcbFlowSet <- function(fcbFlowSet, #flowFrame FCB, output of deskwe_fcbFlowFrame
                       channel, #channel name (char)
                       levels, #number of levels
                       opt = "mixture", #mixture (guassian mixture models) or fisher (univariate k-means)
                       dist = NULL, #for gaussian mixture models, Skew.normal, normal, T.dist
                       subsample = 10e3,
                       trim = 0,
                       ret.model = TRUE,
                       updateProgress = NULL){



  #validation of inputs -------------------------
  if (class(fcbFlowSet) != "fcbFlowSet") {
    stop("Input must be a fcbFlowSet")
  }

  fcbFlowSet.clustered <- fsApply(fcbFlowSet, cluster_fcbFlowFrame,
                                 channel = channel,
                                 levels = levels,
                                 opt = opt,
                                 dist = dist,
                                 subsample = subsample,
                                 trim = trim,
                                 ret.model = ret.model,
                                 updateProgress = updateProgress)

  return(fcbFlowSet(fcbFlowSet.clustered))
}

# plot as histogram and show fit overlay (Ben has code?)
# overlay gaussian fit on histogram
# look at colored count vs PO plot or count vs PB plot
