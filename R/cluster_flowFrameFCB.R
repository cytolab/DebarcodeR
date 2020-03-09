#' Defines populations on barcoded datasets
#'
#'This function allows you to calculate the probability of a cell originating from a given population using
#'either gaussian mixture modeling or jenks natural breaks classification
#'
#' @param flowFrameFCB a flowFrameFCB object with barcoded flowframe and uptake flowframe post deskewing (at least one barcodes slot filled)
#' @param channel The name (string) of the channel to be clustered
#' @param ret.model Option to retain the model for deskewing
#' @param updateProgress used in reactive context (shiny) to return progress information to GUI
#' @param levels integer, the number of barcoding intensities present in the vector
#' @param opt string, either "mixture" (default) for gaussian mixture modeling, or "fisher" for fisher-jenks natural breaks optimization
#' @param dist string in c("normal, skew.Normal, Tdist"), passed to mixsmsn
#' @param subsample Integer, number of cells to subsample, defaults to 10,000
#' @param trim numberic between 0, 1; used to trim the upper and lower extremes to exlcude outliers (eg. trim = 0.01 exludes most extreme 1% of data)
#'
#' @return a flowFrameFCB with deskewed barcodes slot and clustering slot with a matrix of probabilities, with ncol = levels, and nrow = legnth(vec).
#' If gaussian mixture modeling is used the probailities correspond to the probability
#' of the cell originaiting that level under the distrubtion specified by the mixture model
#' If jenks natural breaks optimization is used, the probability is estimated empirically based on a histogram
#'
#' @seealso \code{\link{deskew_flowFrameFCB}}
#' @export
#' @import classInt mixsmsn sn

# aspirational --> v2?
# find sd and mean of uptake - use for probabilities
# bounds as 5th and 95th, separation of each level --> find probability
# uptake for two channels - use as model (Prior)
# READ smsn.mix paper and function

cluster_flowFrameFCB <- function(flowFrameFCB, #flowFrame FCB, output of deskwe_flowFrameFCB
                       channel, #channel name (char)
                       levels, #number of levels
                       opt = "mixture", #mixture (guassian mixture models) or fisher (univariate k-means)
                       dist = NULL, #for gaussian mixture models, Skew.normal, normal, T.dist
                       subsample = 10e3,
                       trim = 0,
                       ret.model = TRUE,
                       updateProgress = NULL){


  # match.arg1 here for options and distributions (normal, skew.normal) - dist and opt
  if(!any(class(flowFrameFCB) == "flowFrameFCB")){
    stop("Input must be an object of class flowFrameFCB")
  }

  if (length(flowFrameFCB@barcodes) == 0) {
    stop(
      "Input must have channels in the barcodes slot that have been run through deskew_flowFrameFCB"
    )
  }


 vec =  flowFrameFCB@barcodes[[which(names(flowFrameFCB@barcodes) == channel)]][["deskewing"]][["values"]]

  quantiles <- quantile(vec, c(trim/2, 1- trim/2))
  vec.trim <- vec[quantiles[1] < vec & quantiles[2] > vec]
  vecss <- sample(vec.trim, subsample, replace = TRUE)

  if(opt == "mixture") {

    if (levels > 1) {
      mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")   #fisher-jenks (breaks in data)
      classif <- sapply(vecss, function(x) pracma::findintervals(x, mod.int$brks))
      classif <- unlist(classif)
      classif <- levels + 1 - classif
      mu.i <- as.numeric(unlist(lapply(split(vecss, classif), median))[-1])
    }
    if (is.function(updateProgress)) {
      updateProgress(detail = "Optimizing Mixture Model")
    }

    Snorm.analysis <- mixsmsn::smsn.mix(vecss, nu = 3, g = levels, criteria = TRUE,
                                        get.init = TRUE, group = TRUE, family = dist, calc.im = FALSE, obs.prob = TRUE,
                                        kmeans.param = list(iter.max = 20, n.start = 10, algorithm = "Hartigan-Wong"))   #sigma 2 parameter = mu.i

    loc <- Snorm.analysis$mu
    scale <- sqrt(Snorm.analysis$sigma2)
    shape <- Snorm.analysis$shape

    Snorm.df <- data.frame(loc, scale, shape)
    probs <- data.frame(x = vec)

    for (i in (1:nrow(Snorm.df))){
      probs[,as.character(i)] <- sn::dsn(vec, dp = as.numeric(Snorm.df[i,]))
    }

    if(levels > 1) {
      probs.scaled <- apply(probs[,-1], 1, function(vec) { vec * Snorm.analysis$pii})
      probs.scaled.df <- as.data.frame(t(probs.scaled))[,rev(order(loc))]
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

  if(ret.model == FALSE){
    flowFrameFCB@barcodes[[which(names(flowFrameFCB@barcodes) == channel)]][[2]] <- list(probabilities = probs.scaled.df)
    names(flowFrameFCB@barcodes[[which(names(flowFrameFCB@barcodes) == channel)]])[2] <-
      "clustering"
  } else{
    flowFrameFCB@barcodes[[which(names(flowFrameFCB@barcodes) == channel)]][[2]] <- list(probabilities = probs.scaled.df, model = Snorm.analysis)
    names(flowFrameFCB@barcodes[[which(names(flowFrameFCB@barcodes) == channel)]])[2] <-
      "clustering"
  }
  return(flowFrameFCB)
}

# plot as histogram and show fit overlay (Ben has code?)
# overlay gaussian fit on histogram
# look at colored count vs PO plot or count vs PB plot
