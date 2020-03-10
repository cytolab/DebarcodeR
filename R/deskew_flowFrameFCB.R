#' Corrects morphology based on scatter or uptake control for flowFrameFCB
#'
#' @param flowFrameFCB a flowFrameFCB object with barcoded flowframe and uptake flowframe, post compensation and preprocessing
#' @param channel The name (string) of the channel to be corrected, ie. the column name in flowFrameFCB barcoded.ff exprs
#' @param method The name of the morphology correction model to use. Choose between earth, lm (linear model), or knijnenburg
#' @param predictors The vector of channel names to be used to build the regression model
#' @param subsample Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.
#' @param ret.model Option to retain the model for deskewing
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI
#'
#' @return a flowFrameFCB with barcode slots added for selected channel corrected for predictors chosen
#' @import earth janitor
#' @export

deskew_flowFrameFCB <- function(flowFrameFCB,
                            channel, #channel name (char)
                            method = "earth", #default to earth
                            predictors = c('fsc_a', 'ssc_a'), #defaults to fsc/ssc
                            subsample = 20e3,
                            ret.model = TRUE,
                            verbose = FALSE,
                            updateProgress = NULL,
                            ...){

  #validation of inputs -------------------------
  if(!any(class(flowFrameFCB) == "flowFrameFCB")){
    stop("Input must be an object of class flowFrameFCB")
  }

  methods <- c("earth", "knijnenburg", "lm")
  method_selected <- match.arg1(method, methods)

  if (identical(flowFrameFCB@barcoded.ff@exprs,flowFrameFCB@uptake.ff@exprs)==TRUE) {
    warning("Barcoded sample being used as uptake control,
            `knijnenburg` method may provide best results'")
  }

  # fcb sample extracted
  fcb = janitor::clean_names(as.data.frame(flowFrameFCB@barcoded.ff@exprs))
  # uptake sample extracted
  uptake = janitor::clean_names(as.data.frame(flowFrameFCB@uptake.ff@exprs))

  # earth model
  if(method_selected == "earth") {
    fcb2 <- morphology_corr.earth(
      fcb = fcb,
      uptake = uptake,
      channel = channel,
      predictors = predictors,
      subsample = subsample,
      ret.model = ret.model,
      updateProgress = updateProgress,
      ...
    )

    # knijnenburg model
  } else if(method_selected == "knijnenburg") {
    fcb2 <- morphology_corr.knijnenburg(
      fcb = fcb,
      uptake = uptake,
      channel = channel,
      fsc_ssc = predictors,
      subsample = subsample,
      ret.model = ret.model,
      updateProgress = updateProgress
    )

    # linear model
  }  else if(method_selected == "lm") {
    fcb2 <- morphology_corr.lm(
      fcb = fcb,
      uptake = uptake,
      channel = channel,
      predictors = predictors,
      ret.model = ret.model,
      slope = 1,
      updateProgress = updateProgress
    )
  }

  # create barcode slots per channel flowFrameFCB with deskewed and model data
  if (length(flowFrameFCB@barcodes) == 0){
    slot(flowFrameFCB, "barcodes") <- list(list(list()))
    flowFrameFCB@barcodes[[1]][[1]] <- fcb2
    names(flowFrameFCB@barcodes)[[1]] <- channel
    names(flowFrameFCB@barcodes[[1]]) <- "deskewing"
  }else{
    if (length(which(names(flowFrameFCB@barcodes) == channel)) == 1) {
      flowFrameFCB@barcodes[[which(names(flowFrameFCB@barcodes) == channel)]][[1]] <- fcb2
      names(flowFrameFCB@barcodes)[[which(names(flowFrameFCB@barcodes) == channel)]] <-
        channel
      names(flowFrameFCB@barcodes[[which(names(flowFrameFCB@barcodes) == channel)]]) <-
        "deskewing"
    }
    else {
      flowFrameFCB@barcodes[[(length(flowFrameFCB@barcodes) + 1)]] <-
        list(fcb2)
      names(flowFrameFCB@barcodes)[[length(flowFrameFCB@barcodes)]] <-
        channel
      names(flowFrameFCB@barcodes[[length(flowFrameFCB@barcodes)]]) <-
        "deskewing"

    }
  }

return(flowFrameFCB)
}

#' Option selector.
#'
#' @param arg a choice or a vector of choices
#' @param choices a character vector of choices
#'
#' @return returns the first choice is no choice is made, otherwise returns
#'  a the choice, or an error if the choice was invaldi

match.arg1 <- function (arg, choices)
{
  if (missing(choices)) {
    formal.arg <-
      formals(sys.function(sys.parent()))[[deparse(substitute(arg))]]
    if (length(formal.arg)==3 && formal.arg[[1]]=="[" &&
        formal.arg[[3]]==1)
      formal.arg <- formal.arg[[2]]
    choices <- eval(formal.arg)
  }
  if (all(arg == choices))
    return(choices[1])
  i <- pmatch(arg, choices)
  if (is.na(i))
    stop(paste("Choice should be one of", paste(choices, collapse = ", "),
               sep = " "))
  if (length(i) > 1)
    stop("there is more than one match in match.arg")
  choices[i]
}
