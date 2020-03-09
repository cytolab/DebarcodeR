#' Corrects morphology based on scatter or uptake control.
#'
#' @param fcb the barcoded dataframe, post compensation and preprocessing
#' @param uptake Optional: a dataframe consisting of all cells barcoded with a single level of the barcoding dye
#' @param channel The name (string) of the channel to be corrected, ie. the column name in 'fcb_df'
#' @param predictors A single predictor
#' @param subsample Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI#'
#'
#' @return a tibble/data.frame with the selected channel corrected for fsc and ssc
#' @export

morphology_corr.lm <- function(fcb,
                               uptake,
                               channel,
                               predictors = NULL,
                               subsample = 10e3,
                               slope = 1,
                               updateProgress = NULL,
                               ret.model = FALSE) {

  if(length(predictors) != 1){
    stop('Please select a single predictor')}

  if(is.null(slope)){
    lm.formula <- as.formula(paste(channel, "~", predictors))
  } else {
    lm.formula <- as.formula(paste0(channel,
                                    " ~ 1 + offset(",
                                    slope, "*",
                                    predictors, ")"))
  }
  lm.model <- lm(lm.formula, data = uptake)
  fcb[,channel]<- fcb[,channel] - predict(lm.model, newdata = fcb) +   median(unlist(fcb[,channel]))

  if(ret.model == FALSE){
    return(list(values = fcb[,channel]))
  } else{
    return(list(values = fcb[,channel],
                model = lm.model))
  }
}
