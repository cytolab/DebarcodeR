#' find the area which contains 95\% of the data
#'
#' @param data text
#' @param fsc text
#' @param ssc text
#' @param subsample text
#' @return some text
#'
#' @seealso text
#' @examples text
# fsc_ssc <- c(fsc = 'FSC-A', ssc = 'SSC-A')
# fsc_ssc[['fsc']]
#fsc_ssc[['fsc']]

selectDenseScatterArea <- function (data,
                                    fsc_ssc = c(fsc = 'FSC-A', ssc = 'SSC-A'),
                                    subsample = 10e3) {
  data <- as.data.frame(data)
  fsc <- fsc_ssc[1]
  ssc <- fsc_ssc[2]
  if(!is.null(subsample) & (nrow(data) > subsample)){
    data <- data[sample(1:nrow(data), subsample),]
  }

  fsc_limits <- c(min(data[,fsc]), max(data[,fsc]))
  ssc_limits <- c(min(data[,ssc]), max(data[,ssc]))

  S <- 8
  area <- 0.95
  FGA <- matrix(0, 2^S, 2^S)
  AP <- matrix(1, 2^S, 2^S)
  #print("31")
  N <- nrow(data)
  #as.data.frame(uptake[,'fsc_a'])

  P <- MASS::kde2d(data[,fsc],
                   data[,ssc],
                   n = 2^S, lims = c(fsc_limits, ssc_limits))
  x_d <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
  y_d <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
  z_d <- P[[3]]
  z_d <- z_d/sum(z_d)
  ft <- function(x) abs(sum(z_d[z_d>x])- area)
  tmin <- optimize(ft, c(0,1e-4), tol = 1e-9)$objective
  AP <- AP*(z_d<tmin)
  FGA <- FGA + z_d
  Y1 <- x_d[!AP]
  Y2 <- y_d[!AP]
  Loc <- cbind(Y1, Y2)
  c <- FGA[!AP]

  return(list("loc" = Loc , "c" = c))

}
