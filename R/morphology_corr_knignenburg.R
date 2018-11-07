#' Corrects morphology based on scatter or uptake control.
#'
#' @param fcb the barcoded dataframe, post compensation and preprocessing
#' @param uptake Optional: a dataframe consisting of all cells barcoded with a single level of the barcoding dye
#' @param channel The name (string) of the channel to be corrected, ie. the column name in 'fcb_df'
#' @param predictors The vector of channel names to be used to build the regression model, must be FSC, SSC.
#' @param subsample Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI#'
#'
#' @return a tibble/data.frame with the selected channel corrected for fsc and ssc
#' @export
morphology_corr.knignenburg <- function(fcb,
                                        uptake,
                                        channel,
                                        fsc_ssc = c(fsc = 'FSC-A', ssc = 'SSC-A'),
                                        subsample = 10e3,
                                        updateProgress = NULL) {

  #print(channel)
  if (is.function(updateProgress)) {
    updateProgress(detail = "Mapping cellular Density")}
  #print("26")
  area_density <- selectDenseScatterArea(uptake,
                                         subsample = subsample,
                                         fsc_ssc = fsc_ssc)
  #print("30")

  if (is.function(updateProgress)) {
    updateProgress(detail = "Performing Morphology Correction")}
  regression.output <- doRegressContrained(uptake,
                                           fcb,
                                           fsc_ssc = fsc_ssc,
                                           Loc = area_density$loc,
                                           weight = area_density$c,
                                           trans = "none",
                                           columns = c(channel),
                                           monodir = c(1,1))
  cor.data <- regression.output[[1]]
  fcb[, channel]<- cor.data[, channel] - median(unlist(cor.data[, channel])) + median(unlist(fcb[, channel]))
  return(fcb)
}


