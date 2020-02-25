deskew_earth <- function(flowFrameFCB,ret.model = FALSE,
                         what = c('x', 'x + se'),
                         nfold = 1,
                         ncross = 0,
                         channel,
                         predictors = c('fsc_a', 'ssc_a'),
                         subsample = 30e3,
                         updateProgress = NULL,
                         ...){

  # fcb sample extracted
  fcb = clean_names(as.data.frame(flowFrameFCB@barcoded.ff@exprs))

  # uptake sample extracted
  uptake = clean_names(as.data.frame(flowFrameFCB@uptake.ff@exprs))

  # making model
  what.options <- c("x", "x + se")
  what <- match.arg(what, what.options)
  print(what)
  print(channel)
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

slot(flowFrameFCB, "barcodes") <- list(list())
flowFrameFCB@barcodes[[1]]<-list(fcb = fcb[,channel],
     model = earth.model)
names(flowFrameFCB@barcodes)[[1]]<-channel


  if(ret.model == FALSE){
    return(fcb)
  } else{
    return(flowFrameFCB)
  }
}
