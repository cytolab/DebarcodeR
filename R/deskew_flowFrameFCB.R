#' Corrects morphology based on scatter or uptake control for flowFrameFCB
#'
#' @param fcb the barcoded dataframe, post compensation and preprocessing
#' @param uptake Optional: a dataframe consisting of all cells barcoded with a single level of the barcoding dye
#' @param channel The name (string) of the channel to be corrected, ie. the column name in 'fcb_df'
#' @param predictors The vector of channel names to be used to build the regression model
#' @param subsample Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI#'
#'
#' @return a tibble/data.frame with the selected channel corrected for fsc and ssc
#' @import earth janitor
#' @export
#' @export

deskew_flowFrameFCB <- function(flowFrameFCB,
                            channel, #channel name (char)
                            method = c("earth", "knignenburg", "lm"), #default to earth
                            predictors = c('fsc_a', 'ssc_a'), #defaults to FSC/SSC
                            subsample = 30e3,
                            ret.model = TRUE,
                            verbose = FALSE,
                            updateProgress = NULL,
                            ...) {

  #validation of inputs -------------------------
  if(!any(class(flowFrameFCB) == "flowFrameFCB")){
    stop("Input must be an object of class flowFrameFCB")
  }

  methods <- c("earth", "knignenburg", "lm")
  method_selected <- match.arg1(method, methods)

  if (identical(flowFrameFCB@barcoded.ff@exprs,flowFrameFCB@uptake.ff@exprs)==TRUE) {
    warning("Barcoded sample being used as uptake control,
            `knignenburg` method may provide best results'")
  }

  # fcb sample extracted
  fcb = janitor::clean_names(as.data.frame(flowFrameFCB@barcoded.ff@exprs))
  # uptake sample extracted
  uptake = janitor::clean_names(as.data.frame(flowFrameFCB@uptake.ff@exprs))

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

  } else if(method_selected == "knignenburg") {
    fcb2 <- morphology_corr.knignenburg(
      fcb = fcb,
      uptake = uptake,
      channel = channel,
      fsc_ssc = predictors,
      subsample = subsample,
      updateProgress = updateProgress
    )

  }  else if(method_selected == "lm") {
    fcb2 <- morphology_corr.lm(
      fcb = fcb,
      uptake = uptake,
      channel = channel,
      predictors = predictors,
      slopes = 1,
      updateProgress = updateProgress
    )
  }

slot(flowFrameFCB, "barcodes") <- list(list())
if(ret.model == TRUE){
  flowFrameFCB@barcodes[[1]]<-fcb2
  names(flowFrameFCB@barcodes)[[1]]<-channel
}

#else flowFrameFCB@barcodes[[i]]<-list(fcb = fcb[,channels[i]])
#names(flowFrameFCB@barcodes)[[i]]<-channels[i]}
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
    stop(paste("Method should be one of", paste(choices, collapse = ", "),
               sep = " "))
  if (length(i) > 1)
    stop("there is more than one match in match.arg")
  choices[i]
}



