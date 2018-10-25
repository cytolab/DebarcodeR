#' Corrects morphology based on scatter or uptake control.
#'
#' @param fcb_df the barcoded dataframe, post compensation and preprocessing
#' @param bc_single_level Optional: a dataframe consisting of all cells barcoded with a single level of the barcoding dye
#' @param channel The name (string) of the channel to be corrected, ie. the column name in 'fcb_df'
#' @param opt String, Morphology correction apporach, either "regression" (default) or "control dye"
#' @param subsample Integer, number of cells to sample (with replacement) for the morphology correction, defaults to 10,000.
#' @param trans String, transformation to apply to the barcoding channels, defaults to arcsinh, also log10
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI
#' @param cofactor_bc1 Numeirc, Cofactor used for the arcsinh transformation on the barcoding channels
#' @param cofactor_update Numeric, Cofactor used for the arcsinh transformation of the uptake channel
#' @param uptake_channel String, uptake channel if "control dye" option is used as morphology correction method
#'
#'
#' @return numeric vector of length nrow(fcb_df) representing the morphology corrected channel
#'
#' @seealso \code{\link{selectDenseScatterArea}} \code{\link{doRegressConstrained}}
#' @export
#' @examples

morphology_corr <- function(fcb_df, bc_single_level = NULL, channel,
                            opt = "regression", subsample = 10e3, trans = 'arcsinh',
                            updateProgress = NULL,
                            cofactor_bc1 = NULL, cofactor_uptake = NULL,
                            uptake_channel = NULL) {
  # fcb_df
  # bc_single_level <- fcb_df
  # channel <- "Pacific-Orange-A"
  # levels <- 6
  # cofactor_bc1 <- 150
  # opt <- "regression"
  # updateProgress <- NULL
  # subsample <- 10e3
  # trans  <- "arcsinh"
  if(opt == "regression") {
    if (is.function(updateProgress)) {
      updateProgress(detail = "Mapping cellular Density")
    }

    area_density <- selectDenseScatterArea(bc_single_level, subsample = subsample)

    if (is.function(updateProgress)) {
      updateProgress(detail = "Performing Morphology Correction")
    }

    regression.output <- doRegressContrained(bc_single_level, fcb_df, Loc = area_density$loc, weight = area_density$c,
                                             trans = trans, columns = c(channel), monodir = c(1,1), cofactor = cofactor_bc1)
    cor.data <- regression.output[[1]]
    fcb_df2 <- fcb_df
    fcb_df2[, channel]<- cor.data[, channel]
    ####################################
  } else if(opt == "controldye") {
    if(is.null(uptake_channel)){
      warning("No uptake control channel specified")
      return()
    }
    if(is.null(cofactor_uptake) & trans == "arcsinh"){
      warning("No uptake control cofactor specified")
      return()
    }
    if (trans == "log10"){
      Y <- log10(fcb_df[,channel])
      X <- log10(fcb_df[,uptake_channel])
    } else if (trans == "arcsinh"){
      Y <- asinh(fcb_df[,channel]/cofactor_bc1)
      X <- asinh(fcb_df[,uptake_channel]/cofactor_uptake)
    }

    vec <- Y - X
    vec <- vec + mean(Y)-mean(vec)

    if (trans == "log10"){
      vec <- vec^10
    } else if (trans == "arcsinh"){
      vec <- sinh(vec)*cofactor_bc1
    }
    fcb_df2 <- fcb_df
    fcb_df2[, channel]<- vec


  }
  return(fcb_df2[, channel])
}
