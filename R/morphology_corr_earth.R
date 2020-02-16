#' Corrects morphology using multivariate adaptive regression splines as implemented in the earth pacakge
#'
#' @param fcb the barcoded dataframe, post compensation and preprocessing
#' @param uptake Optional: a dataframe consisting of all cells barcoded with a single level of the barcoding dye
#' @param channel The name (string) of the channel to be corrected, ie. the column name in 'fcb_df'
#' @param predictors The vector of channel names to be used to build the regression model
#' @param subsample Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI#'
#'
#' @return a tibble/data.frame with the selected channel corrected for fsc and ssc
#' @import earth
#' @export
#'
#'
morphology_corr.earth <- function(fcb,
                                  uptake,
                                  ret.model = FALSE,
                                  what = c('x', 'x + se'),
                                  nfold = 1,
                                  ncross = 0,
                                  channel,
                                  predictors = c('FSC-A', 'SSC-A'),
                                  subsample = 30e3,
                                  updateProgress = NULL,
                                  ...) {
  what.options <- c("x", "x + se")
  what <- match.arg1(what, what.options)
  print(what)
  #print(channel)
  if (is.function(updateProgress)) {
    updateProgress(detail = "Training adaptive splines...")}

  lhs <- paste0("`", predictors, "`", collapse = " + ")
  earth.formula <- paste(channel, '~', lhs)
  if(what == "x"){
    earth.model <- earth(as.formula(earth.formula),
                         degree = 2,
                         nprune = 21,
                         nfold = nfold,
                         ncross = ncross,
                         keepxy = TRUE,
                         data = uptake,
                         ...)

    if (is.function(updateProgress)) {
      updateProgress(detail = "Fitting fcb data...")}
    fcb[,channel] <- fcb[,channel] - predict(earth.model, fcb) + median(unlist(uptake[,channel]))

  } else if(what == "x + se") {
    earth.model <- earth(as.formula(earth.formula),
                         degree = 2,
                         nprune = 21,
                         nfold = nfold,
                         ncross = ncross,
                         keepxy = TRUE,
                         varmod.method = "x.earth",
                         data = uptake,
                         trace = 0.3)

    fcb[,channel] <- fcb[,channel] - predict(earth.model, fcb) + median(unlist(uptake[,channel]))
    fcb[,paste0(channel,"se")]<- predict(earth.model, newdata = fcb, interval = "se")

  }


  if(ret.model == FALSE){
    return(fcb)
  } else{
    return(list(fcb = fcb[,channel],
                model = earth.model))
  }
}


