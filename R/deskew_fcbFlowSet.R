#' Corrects morphology based on scatter or uptake control for fcbFlowFrame
#'
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe, post compensation and preprocessing
#' @param channel The name (string) of the channel to be corrected, ie. the column name in fcbFlowFrame barcoded.ff exprs
#' @param method The name of the morphology correction model to use. Choose between earth, lm (linear model), or knijnenburg
#' @param predictors The vector of channel names to be used to build the regression model
#' @param subsample Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.
#' @param ret.model Option to retain the model for deskewing
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI
#'
#' @return a fcbFlowFrame with barcode slots added for selected channel corrected for predictors chosen
#' @import earth janitor
#' @export

deskew_fcbFlowSet <- function(fcbFlowSet,
                                uptake = NULL,
                                channel,
                                #channel name (char)
                                method = "earth",
                                #default to earth
                                predictors = c('fsc_a', 'ssc_a'),
                                #defaults to fsc/ssc
                                subsample = 20e3,
                                ret.model = TRUE,
                                verbose = FALSE,
                                updateProgress = NULL,
                                ...)
{


  #validation of inputs -------------------------
  if (class(fcbFlowSet) == "flowSet") {
    fcbFlowSet <- fcbFlowSet(fcbFlowSet)
  } else if (class(fcbFlowSet) != "fcbFlowSet") {
    stop("Input must be a flowSet or fcbFlowSet")
  }

  fcbFlowSet.deskewed <- fsApply(fcbFlowSet, deskew_fcbFlowFrame,
          uptake = uptake, channel = channel,
          method = method,
          predictors = predictors,
          subsample = subsample,
          ret.model = ret.model,
          verbose = verbose,
          updateProgress = updateProgress,
          ...)

return(fcbFlowSet(fcbFlowSet.deskewed))
}
