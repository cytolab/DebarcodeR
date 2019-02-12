#' Defines populations on barcoded datasets
#'
#'This function allows you to calculate the probability of a cell originating from a given population using
#'either gaussian mixture modeling or jenks natural breaks classificatio.n
#'
#' @param vec a vector of barcode intensities, usually post morphology correction
#' @param levels integer, the number of barcoding intensities present in the vector
#' @param opt string, either "mixture" (default) for gaussian mixture modeling, or "fisher" for fisher-jenks natural breaks optimization
#' @param dist string in c("normal, skew.Normal, Tdist"), passed to mixsmsn
#' @param trans string, transformation to apply to the barcoding channels, defaults to arcsinh, also log10
#' @param cofactor Numeirc, Cofactor used for the arcsinh transformation on vec
#' @param subsample Integer, number of cells to subsample, defaults to 10,000
#' @param trim numberic between 0, 1; used to trim the upper and lower extremes to exlcude outliers (eg. trim = 0.01 exludes most extreme 1% of data)
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI
#'
#'
#' @return a matrix of probabilities, with ncol = levels, and nrow = legnth(vec).
#' If gaussian mixture modeling is used the probailities correspond to the probability
#' of the cell originaiting that level under the distrubtion specified by the mixture model
#' If jenks natural breaks optimization is used, the probability is estimated empirically based on a histogram
#'
#' @seealso \code{\link{morphology_corr}}
#' @export
#' @examples
#'

fit_models <- function(vec, #vector of barcoding intensities, output of morphology.corr
                       levels, #number of levels
                       opt = "mixture", #mixture (guassian mixture models) or fisher (univariate k-means)
                       dist = NULL, #for gaussian mixture models, Skew.normal, normal, T.dist
                       trans = "arcsinh", #asinh or log10
                       cofactor = NULL,
                       subsample = 10e3,
                       trim = 0,
                       updateProgress = NULL){#cofactor for asinh transofrmation

  if (trans == "log10"){
    vec <- log10(vec)
  } else if (trans == "arcsinh"){
    vec <- asinh(vec/cofactor)
  }

  quantiles <- quantile(vec, c(trim/2, 1- trim/2))
  vec.trim <- vec[quantiles[1] < vec & quantiles[2] > vec]
  vecss <- sample(vec.trim, subsample, replace = TRUE)

  if(opt == "mixture") {
    #cofactor <- 150


    if (levels > 1) {
      mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")
      classif <- sapply(vecss, function(x) pracma::findintervals(x, mod.int$brks))
      classif <- levels + 1 - classif
      mu.i <- as.numeric(unlist(lapply(split(vecss, classif), median))[-1])
    }
    if (is.function(updateProgress)) {
      updateProgress(detail = "Optimizing Mixture Model")
    }
    # dist <- "Skew.normal"
    # levels <- 4
    # ?mixsmsn::smsn.mix
    # ptm <- proc.time()
    Snorm.analysis <- mixsmsn::smsn.mix(vecss, nu = 3, g = levels, criteria = TRUE,
                                        get.init = TRUE, group = TRUE, family = dist, calc.im = FALSE, obs.prob = TRUE,
                                        kmeans.param = list(iter.max = 20, n.start = 10, algorithm = "Hartigan-Wong"))
    # proc.time() - ptm
    # mix.hist(vecss, Snorm.analysis, breaks = 50)
    # (Snorm.analysis)

    loc <- Snorm.analysis$mu
    scale <- sqrt(Snorm.analysis$sigma2)
    shape <- Snorm.analysis$shape

    Snorm.df <- data.frame(loc, scale, shape)
    probs <- data.frame(x = vec)

    for (i in (1:nrow(Snorm.df))){
      probs[,as.character(i)] <- sn::dsn(vec, dp = as.numeric(Snorm.df[i,]))
    }

    #print(str(probs))

    if(levels > 1) {
      probs.scaled <- apply(probs[,-1], 1, function(vec) { vec * Snorm.analysis$pii})
      probs.scaled.df <- as.data.frame(t(probs.scaled))[,rev(order(loc))] #order(loc) reorders the columns in descending order
    } else {
      probs.scaled <- probs[,-1]
      probs.scaled.df <- as.data.frame(probs.scaled)
    }


  } else if(opt == "fisher") {
    mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")

    classif <- lapply(vec, FUN = function(x) {findInterval(x, mod.int$brks)})

    classif <- unlist(classif)

    classif <- levels + 1 - classif
    classif[classif > levels] <- 0


    # ggplot(max(vec.split[[2]]))
    # preplot <-ggplot(data.frame(x = vec.split[[3]]), aes(x = x)) +
    #   geom_histogram(bins = 100, col = "black", fill = "grey99") +
    #   geom_vline(data = data.frame(x = mod.int$brks), aes(xintercept = x),
    #              linetype = 2, size = 1)
    # preplot <- ggplot_build(preplot)
    # preplot$data[[1]][["count"]]
    # ggplot(data.frame(x = vec.split[[3]]), aes(x = x)) +
    #   geom_histogram(bins = 100, col = "black", fill = "grey99") +
    #   geom_vline(data = data.frame(x = mod.int$brks), aes(xintercept = x),
    #              linetype = 2, size = 1) +
    #   geom_hline(yintercept = max(preplot$data[[1]][["count"]])/8,
    #              linetype = 4, size = 1, col = "blue") +
    #   theme_classic()

    #calculate the empircal probability for each cell belonging to each
    #population based on the histogram

    vec.split <- split(vec, classif)

    hist.probs <- list()
    for (i in as.character(1:max(as.numeric(names(vec.split))))){
      myhist <- hist(vec.split[[i]],100, plot = FALSE)
      binprobs <- myhist$counts/sum(myhist$counts)
      hist.probs.i<- rep(0, times = length(vec))
      bin.assingments <- findInterval(vec, myhist$breaks)
      hist.probs.i[which(bin.assingments != 0)] <- binprobs[bin.assingments]
      hist.probs.i[which(is.na(hist.probs.i))] <- 0
      hist.probs[[i]] <- hist.probs.i
    }
    hist.probs.m <- do.call(cbind, hist.probs)
    probs.scaled.df <- as.data.frame(hist.probs.m)

  }

  if (levels > 1) {
   # colMax <- apply(probs.scaled.df, 2, max)
    #print(str(colMax))
   # rowMax <- apply(probs.scaled.df, 1, max)
    #print(str(rowMax))
    #probs.rescale.col <- sweep(probs.scaled.df, 2, colMax, FUN="/")
    #probs.rescale.row <- sweep(probs.scaled.df, 1, rowMax, FUN="/")

  } else {
   # colMax <- max(probs.scaled.df)
   # probs.rescale.col <- probs.scaled.df/colMax
  }
  return(probs.scaled.df)
}
