#' Defines populations on barcoded datasets
#'
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing, clustering, and assignment
#' @param ambiguitycut numeric from 0 to 1, threshhold below which to discard ambigious cells, eg: 0.02,
#'  discards cells with more than 2% chance of originating from another population
#' @param subsample number, number of cells to subsample during m step
#' @return a fcbFlowFrame object with a barcode slot filled with deskewing, clustering, cell assignment as
#' a vector of integers from 0:ncol(probs), cells assigned a classification of 0 remained unassigned,
#' otherwise number corresponds to the barcoding level assignment of that cell
#' @export
#' @import mclust matrixStats
#' @importFrom data.table tstrsplit
em_optimize <- function(fcbFlowFrame,
                        modelName = "VEE",
                        subsample  = 50e3,
                        verbose = TRUE,
                        niter = 1,
                        shrinkage = 0.05) {

  if (class(fcbFlowFrame) != "fcbFlowFrame") {
    stop("Input must be a fcbFlowFrame")
  }
  if (verbose) cat("Initializing...\n")

  barcodes <- fcbFlowFrame@barcodes[!(names(fcbFlowFrame@barcodes) == "wells")]
#  bc_initial.df <- as.data.frame(lapply(lapply(barcodes,`[[`, "assignment"), `[[`, "values"))
  #assigned.ind <- which(!apply(bc_initial.df == 0, 1 ,any))
  deskewed_cols <- as.data.frame(lapply(lapply(barcodes,`[[`, "deskewing"), `[[`, "values"))
  deskewed_models <- lapply(lapply(barcodes,`[[`, "clustering"), `[[`, "model")
  deskewed_probs  <- lapply(lapply(lapply(barcodes,`[[`, "clustering"), `[[`, "probabilities"), as.matrix)
  # as_tibble(x = barcodes$pacific_blue$deskewing$values) %>%
  #   bind_cols(as_tibble(barcodes$pacific_blue_a$clustering$probabilities)) %>%
  #   sample_n(1e3) %>%
  #   gather(cat, prob, -value) %>%
  #   ggplot(aes(x=value, y = prob, col = cat))  +
  #   geom_point()
  deskewed.dims <- lapply(deskewed_probs, ncol)
  deskewed.dims.vec <- as.numeric(deskewed.dims)
  deskewed.rep.mat <- kronecker(deskewed.dims.vec, t(rep(1, length(deskewed.dims.vec))))
  deskewed.rep.mat.each <- deskewed.rep.mat
  deskewed.rep.mat.times  <- deskewed.rep.mat
  deskewed.rep.mat.times[upper.tri(deskewed.rep.mat.times, diag = T)] <- 1
  deskewed.rep.mat.each[lower.tri(deskewed.rep.mat.times, diag = T)] <- 1
  deskewed.reps <- split(cbind(apply(deskewed.rep.mat.each, 2, prod),
                               apply(deskewed.rep.mat.times, 2, prod)),
                         seq(length(deskewed.dims.vec)))
  # deskewed.reps <- list()
  # for (i in seq(length(deskewed.dims.vec))) {
  #   deskewed.reps[[names(deskewed.dims)[i]]][['times']] <- prod(deskewed.dims.vec[seq(length(deskewed.dims.vec)) > i])
  #   deskewed.reps[[names(deskewed.dims)[i]]][['each']] <- prod(deskewed.dims.vec[seq(length(deskewed.dims.vec)) < i])
  # }
  # })
  #

  deskewed_names_expanded <-
    mapply(function(mat, reps) {
      t(rep(1, reps[[1]])) %x% mat %x% t(rep(1, reps[[2]]))
    },
    lapply(lapply(deskewed.dims, seq), `t`),
    deskewed.reps, SIMPLIFY = FALSE)
  deskewed_names_expanded$sep <- "."
  cl.names <- do.call(paste, deskewed_names_expanded)
  deskewed_probs_expanded <-
    mapply(function(mat, reps) {
      t(rep(1, reps[[1]])) %x% mat %x% t(rep(1, reps[[2]]))
    },
    deskewed_probs,
    deskewed.reps, SIMPLIFY = FALSE)
  cl.mat <- do.call('*', deskewed_probs_expanded)
  cl.mat <- cl.mat/rowSums(cl.mat)
  n <- round(nrow(deskewed_cols)*shrinkage)
  data.i <- deskewed_cols#[assigned.ind, ]
  # combn()
  #
  # #str((as.matrix(deskewed_probs$pacific_blue_a)) %o% (as.matrix(deskewed_probs$pacific_orange_a)))
  # #x <- 1:9; names(x) <- x
  if (shrinkage > 0) {
    shrinkage.list <- lapply(deskewed_models, function(mod) {
      mixsmsn::rmix(
        n,
        mod$pii,
        family = class(mod),
        arg = split(data.frame(mod[c("mu", "sigma2", "shape", "nu")]), seq(length(mod$pii))),
        cluster = TRUE
      )
    })
    shrinkage.y <- as.data.frame(lapply(shrinkage.list, `[[`, "y"))
    shrinkage.cl <- as.data.frame(lapply(shrinkage.list, `[[`, "cluster"))
    shrinkage.cl.mat <- mclust::unmap(as.factor(apply(shrinkage.cl, 1, paste0, collapse = ".")))
    data.ii <- rbind(data.i, shrinkage.y)
    cl.mat <-
      rbind(cl.mat, shrinkage.cl.mat)
  } else {
    cl.mat <- cl.mat
    data.ii <- data.i
  }
#  nrow(bc_initial.df_assigned)

#  bc_initial.f_assigned <- as.factor(apply(bc_initial.df_assigned, 1, paste0, collapse = "."))

#  cl <- bc_initial.f_assigned
#  cl.mat <- mclust::unmap(bc_initial.f_assigned)

  for (i in seq(niter)) {
    if (verbose) cat(paste("EM Round", i, "of", niter, "..."))
    ## EM fitting -------------------------------------------------------------
    msEst <- mclust::mstep(modelName = modelName,
                   data = data.ii,
                   z = cl.mat)
    esEst <- mclust::estep(modelName = msEst$modelName,
                   data = deskewed_cols,
                   parameters = msEst$parameters)
    if (verbose) cat(paste0(" loglik: ", round(esEst$loglik), "\n"))
  #  print(esEst$loglik)
    cl.mat <- esEst$z

    ## shrinkage---------------------------------------------------------------
    if (shrinkage > 0) {
      simdata <- mclust::sim(esEst$modelName, esEst$parameters, n)
      colnames(simdata)[-1] <- colnames(deskewed_cols)
      data.ii <- rbind(deskewed_cols, simdata[,-1])
      cl.mat <- rbind(cl.mat, unmap(simdata[,1]))
    } else {
      data.ii <- deskewed_cols
    }
  }
  if (verbose) cat(paste("Calculating probabilities..."))

  mvgmm.cdens <- esEst$parameters$pro * cdens(esEst$modelName, deskewed_cols, parameters = esEst$parameters) # component densities
  mvgmm.tdens <- dens(esEst$modelName, deskewed_cols, parameters = esEst$parameters) # total density
  probs <- mvgmm.cdens/sum(mvgmm.tdens) # density normalized to 1
  colnames(probs) <- cl.names
  # # Assignment------------------------------------------------------------------
  # cl <- apply(probs, 1, which.max)
  # cl <- cl.names[cl]
  #
  # probs.norm.row <- calculate.ambiguity(probs)
  # probs.norm.col <- calculate.likelihood(probs)
  #
  # non.ambigious <- apply(probs.norm.row, 1, max) > (1 - ambiguitycut)
  # cl[!non.ambigious] <- paste0(rep(0, ncol(deskewed_cols)), collapse = ".")
  # likely <- apply(probs.norm.col > 1/likelihoodcut, 1, any)
  # cl[!likely] <- paste0(rep(0, ncol(deskewed_cols)), collapse = ".")
 # new.assignments.l <- data.table::tstrsplit(cl, ".", fixed = TRUE, names = colnames(deskewed_cols), type.convert = T)

#
#   fcbFlowFrame@barcodes <- mapply(function(bc, assignments) {
#     bc[['assignment']][['values']] <- assignments
#     return(bc)
#   },
#   fcbFlowFrame@barcodes,
#   new.assignments.l,
#   SIMPLIFY = FALSE)

  # fcbFlowFrame@barcodes[['wells']][['assignment']][['values']] <- cl
  # fcbFlowFrame@barcodes[['wells']][['assignment']][['ambiguity']] <- ambiguitycut
  # fcbFlowFrame@barcodes[['wells']][['assignment']][['likelihoodcut']] <- likelihoodcut
  #
  #
  # Clustering------------------------------------------------------------------
  fcbFlowFrame@barcodes[['wells']][['clustering']][['probabilities']] <- probs
  fcbFlowFrame@barcodes[['wells']][['clustering']][['channels']] <- deskewed.dims
  fcbFlowFrame@barcodes[['wells']][['clustering']][['model']] <- esEst
  # fcbFlowFrame@barcodes[['wells']][['assignment']][['dims']] <- colnames(bc_initial.df)
  if (verbose) cat(paste("Done!"))

  return(fcbFlowFrame)
}
