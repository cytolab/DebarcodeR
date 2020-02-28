#' Corrects morphology using multivariate adaptive regression splines as implemented in the earth pacakge for flowFrameFCB
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
deskew_earth <- function(flowFrameFCB,ret.model = FALSE,
                         what = c('x', 'x + se'),
                         nfold = 1,
                         ncross = 0,
                         channels,
                         predictors = c('fsc_a', 'ssc_a'),
                         subsample = 30e3,
                         updateProgress = NULL,
                         ...){

  # fcb sample extracted
  fcb = janitor::clean_names(as.data.frame(flowFrameFCB@barcoded.ff@exprs))

  # uptake sample extracted
  uptake = janitor::clean_names(as.data.frame(flowFrameFCB@uptake.ff@exprs))

  slot(flowFrameFCB, "barcodes") <- list(list())

for (i in 1:length(channels)){
  # making model
  what.options <- c("x", "x + se")
  what <- match.arg1(what, what.options)
  print(what)
  print(channels[i])
  if (is.function(updateProgress)) {
    updateProgress(detail = "Training adaptive splines...")}
  lhs <- paste0("`", predictors, "`", collapse = " + ")
  earth.formula <- paste(channels[i], '~', lhs)
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
    fcb[,channels[i]] <- fcb[,channels[i]] - predict(earth.model, fcb) + median(unlist(uptake[,channels[i]]))

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

    fcb[,channels[i]] <- fcb[,channels[i]] - predict(earth.model, fcb) + median(unlist(uptake[,channels[i]]))
    fcb[,paste0(channels[i],"se")]<- predict(earth.model, newdata = fcb, interval = "se")

  }
  if(ret.model == TRUE){
flowFrameFCB@barcodes[[i]]<-list(fcb = fcb[,channels[i]],
     model = earth.model)
names(flowFrameFCB@barcodes)[[i]]<-channels[i]
}else flowFrameFCB@barcodes[[i]]<-list(fcb = fcb[,channels[i]])
  names(flowFrameFCB@barcodes)[[i]]<-channels[i]}

    return(flowFrameFCB)
}
