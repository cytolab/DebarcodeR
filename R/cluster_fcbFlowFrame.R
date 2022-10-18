#' Defines populations on barcoded datasets
#'
#'This function allows you to calculate the probability of a cell originating from a given population using
#'either gaussian mixture modeling or jenks natural breaks classification
#'
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing (at least one barcodes slot filled)
#' @param channel The name (string) of the channel to be clustered
#' @param ret.model Option to retain the model for deskewing
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI
#' @param levels integer, the number of barcoding intensities present in the vector
#' @param opt string, either "mixture" (default) for gaussian mixture modeling, or "fisher" for fisher-jenks natural breaks optimization, or "manual.breaks" for manual specification of breakpoints
#' @param dist string in c("Normal, Skew.normal, Tdist"), passed to mixsmsn
#' @param subsample Integer, number of cells to subsample, defaults to 10,000
#' @param manbreaks Vector, length levels + 1, used to specifiy manual breakpoints for levels
#' @param trim numberic between 0, 1; used to trim the upper and lower extremes to exlcude outliers (eg. trim = 0.01 exludes most extreme 1\% of data)
#'
#' @return a fcbFlowFrame with deskewed barcodes slot and clustering slot with a matrix of probabilities, with ncol = levels, and nrow = legnth(vec).
#' If gaussian mixture modeling is used the probailities correspond to the probability
#' of the cell originaiting that level under the distrubtion specified by the mixture model
#' If jenks natural breaks optimization is used, the probability is estimated empirically based on a histogram
#'
#' @seealso \code{\link{deskew_fcbFlowFrame}}
#' @export
#' @import classInt mixsmsn sn

# aspirational --> v2?
# find sd and mean of uptake - use for probabilities
# bounds as 5th and 95th, separation of each level --> find probability
# uptake for two channels - use as model (Prior)
# READ smsn.mix paper and function
cluster_fcbFlowFrame <- function(fcbFlowFrame, #flowFrame FCB, output of deskwe_fcbFlowFrame
                       channel, #channel name (char)
                       levels, #number of levels
                       opt = "mixture", #mixture (guassian mixture models) or fisher (univariate k-means)
                       dist = NULL, #for gaussian mixture models, Skew.normal, normal, T.dist
                       subsample = 3e3,
                       trim = 0,
                       ret.model = TRUE,
                       manbreaks = NULL,
                       updateProgress = NULL){


  # match.arg1 here for options and distributions (normal, skew.normal) - dist and opt
  if (!any(class(fcbFlowFrame) == "fcbFlowFrame")) {
    stop("Input must be an object of class fcbFlowFrame")
  }

  if (length(fcbFlowFrame@barcodes) == 0) {
    stop(
      "Input must have channels in the barcodes slot that have been run through deskew_fcbFlowFrame"
    )
  }

  options <- c("mixture","fisher", "manual.breaks")
  opt_selected <- match.arg1(opt, options)

  distributions <- c("Normal","Skew.normal","Tdist")
  dist_selected <- match.arg1(dist, distributions)


  vec <-  fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][["deskewing"]][["values"]]

  quantiles <- quantile(vec, c(trim/2, 1 - trim/2))
  vec.trim <- vec[quantiles[1] < vec & quantiles[2] > vec]
  vecss <- sample(vec.trim, subsample, replace = TRUE)

  if (opt_selected == "mixture") {

    if (levels > 1) {
      #?classInt::classIntervals
      mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")   #fisher-jenks (breaks in data)
      classif <- sapply(vecss, function(x) pracma::findintervals(x, mod.int$brks))
      classif <- unlist(classif)
      classif <- levels + 1 - classif
      mu.i <- as.numeric(unlist(lapply(split(vecss, classif), median))[-1])
    }
    if (is.function(updateProgress)) {
      updateProgress(detail = "Optimizing Mixture Model")
    }
    Snorm.analysis <- mixsmsn::smsn.mix(vecss, nu = 3,
                                        g = levels,
                                        shape = rep(0, levels),
                                        mu = mu.i,
                                        get.init = TRUE, group = TRUE, family = dist_selected, calc.im = FALSE,
                                        kmeans.param = list(iter.max = 20, n.start = 10, algorithm = "Hartigan-Wong"))   #sigma 2 parameter = mu.i
    loc <- Snorm.analysis$mu
    scale <- sqrt(Snorm.analysis$sigma2)
    shape <- Snorm.analysis$shape
    Snorm.df <- data.frame(loc, scale, shape)
    probs.x <- data.frame(x = vec)
    probs.y <- apply(Snorm.df, 1, function(i) sn::dsn(vec, dp = as.numeric(i)))
    colnames(probs.y) <- as.character(1:nrow(Snorm.df))
    probs <- cbind(probs.x, probs.y)
    #for (i in (1:nrow(Snorm.df))) {
     # probs[,as.character(i)] <- sn::dsn(vec, dp = as.numeric(Snorm.df[i,]))
    #}

    if (levels > 1) {
      probs.scaled <- t(t(as.matrix(probs[,-1])) * Snorm.analysis$pii)
      probs.scaled.df <- as.data.frame(t(probs.scaled))[,rev(order(loc))]
    } else {
      probs.scaled <- probs[,-1]
      probs.scaled.df <- as.data.frame(probs.scaled)
    }

  } else if (opt_selected == "fisher") {
    mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")

    classif <- lapply(vec, FUN = function(x) {findInterval(x, mod.int$brks)}) ##

    classif <- unlist(classif)

    classif <- levels + 1 - classif
    classif[classif > levels] <- 0

    vec.split <- split(vec, classif)

    hist.probs <- list()
    for (i in as.character(1:max(as.numeric(names(vec.split))))) {
      myhist <- hist(vec.split[[i]],100, plot = FALSE)
      binprobs <- myhist$counts/sum(myhist$counts)
      hist.probs.i <- rep(0, times = length(vec))
      bin.assingments <- findInterval(vec, myhist$breaks)
      hist.probs.i[which(bin.assingments != 0)] <- binprobs[bin.assingments]
      hist.probs.i[which(is.na(hist.probs.i))] <- 0
      hist.probs[[i]] <- hist.probs.i
    }
    probs <- do.call(cbind, hist.probs)
    probs.scaled <- probs
  #  probs.scaled.df <- as.data.frame(hist.probs.m)

  } else if (opt_selected == "manual.breaks") {

    if(is.null(manbreaks)) {stop("Please specifiy manual breakpoints")}
    if(length(manbreaks) != (levels + 1)) {stop("Please ensure breakpoints n+1 breakpoints are provided as a vector")}

    classif <- lapply(vec, FUN = function(x) {findInterval(x, manbreaks)}) ##
    classif <- unlist(classif)

    classif <- levels + 1 - classif
    classif[classif > levels] <- 0

    vec.split <- split(vec, classif)

    hist.probs <- list()
    for (i in as.character(1:max(as.numeric(names(vec.split))))) {
      myhist <- hist(vec.split[[i]],100, plot = FALSE)
      binprobs <- myhist$counts/sum(myhist$counts)
      hist.probs.i <- rep(0, times = length(vec))
      bin.assingments <- findInterval(vec, myhist$breaks)
      hist.probs.i[which(bin.assingments != 0)] <- binprobs[bin.assingments]
      hist.probs.i[which(is.na(hist.probs.i))] <- 0
      hist.probs[[i]] <- hist.probs.i
    }
    probs <- do.call(cbind, hist.probs)
    probs.scaled <- probs
  }

  if (ret.model == TRUE & opt == "mixture") {
    fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][[2]] <- list(probabilities = probs.scaled, model = Snorm.analysis)
    names(fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]])[2] <-
      "clustering"
  } else if (ret.model == TRUE & opt == "fisher") {
    fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][[2]] <- list(probabilities = probs.scaled, model = mod.int)
    names(fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]])[2] <-
      "clustering"
  } else{
    fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]][[2]] <- list(probabilities = probs.scaled)
    names(fcbFlowFrame@barcodes[[which(names(fcbFlowFrame@barcodes) == channel)]])[2] <-
      "clustering"
  }
  return(fcbFlowFrame)
}

# plot as histogram and show fit overlay (Ben has code?)
# overlay gaussian fit on histogram
# look at colored count vs PO plot or count vs PB plot
