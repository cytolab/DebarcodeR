#' Defines populations on barcoded datasets
#'
#' @param fcb_df dataframe, representing the barcoded data, columns for each channel, rows for each cell
#' @param probs matrix of probabilites under population models, returned by 'fit_models'
#' @param likelihoodcut numeric, a likelihood cutoff for discarding unlikely cells, less than 1/k as likely as the most likely cell from that population
#' @param ambiguitycut numeric from 0 to 1, threshhold below which to discard ambigious cells, eg: 0.02,
#'  discards cells with more than 2% chance of originating from another population
#' @param output string, default "classif," also "plot;" what to return?
#' @param channel string, channel which is barcoded, only needed if output == "plot"
#'
#' @return a vector of integers from 0:ncol(probs), cells assigned a classification of 0 remained unassigned,
#'  otherwise number corresponds to the barcoding level assignment of that cell
#' @export
#' @examples
#'

assign_cells <- function(fcb_df, probs, likelihoodcut = 8 , ambiguitycut = 0.02, 
                         output = "classif", channel = NULL){


  row.max <-  apply(probs, 1, sum)
  probs.norm.row <- sweep(probs, 1, row.max, FUN="/")

  col.max <-  apply(probs, 2, max)
  probs.norm.col <- sweep(probs, 2, col.max, FUN="/")


  classif <- rep(0, nrow(fcb_df))

  if(ncol(probs) > 1) { # if assigning more than one level
    classif <- as.numeric(apply(probs.norm.row, 1, which.max))
  } else { # if assigning only one level
    classif <- rep(1, nrow(probs))
  }

  classif[which(is.na(classif))] <- 0 #catches the few cells with 0 probability of belonging to any pop

  # length(as.numeric(apply(probs.norm.row, 1, function(vec){which(vec == 1)})))
  #
  # as.numeric(lapply(apply(probs.norm.row, 1, which.max), unname))
  #
  # length(unlist(apply(probs.norm.row, 1, which.max)))
  # str(t(apply(probs, 1, function(vec) {vec/sum(vec)})))'
  likely <- probs.norm.col > 1/likelihoodcut

  if(ncol(probs) > 1) { # if assigning more than one level
    likely.sum <- apply(likely, 1, sum) #converts logical to numeric
    print(paste0(round(sum(apply(likely, 1, any))/nrow(likely)*100, 3), "% above likelihood cutoff"))
    classif[which(likely.sum < 1)] <- 0
    #print(head(probs.norm.row))
    non.ambigious <- apply(probs.norm.row, 1, max) > (1 - ambiguitycut)
    classif[which(!non.ambigious)] <- 0
  } else {
    likely.sum <- as.numeric(likely)
    #print(str(likely.sum))
    classif[which(likely.sum != 1)] <- 0
  }

  if(output == "classif") {
    return(classif)
  } else {
   # print("reached 61")
    vec <- fcb_df[,channel]
    plot_df <- data.frame(vec = fcb_df[,channel], classif = factor(classif))
    preplot<- ggplot2::ggplot(plot_df, ggplot2::aes(x = vec, fill = classif)) +
      ggplot2::geom_histogram(bins = 400) +
      ggplot2::scale_x_continuous(trans = asinh_trans(150)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::theme_classic()
    preplot
    preplot.build<- ggplot2::ggplot_build(preplot)
    preplot.data <- preplot.build$data[[1]]
    preplot.data.list <- split(preplot.data, preplot.data$group)
    #print("reached 73")
    bounds <- list()
    for (i in 1:ncol(probs)) {
      bounds[[i]]<- c(classif = i, min = min(vec[likely[,i]]), max = max(vec[likely[,i]]))
    }
    #print(str(bounds))
    if(ncol(probs) > 1) { # if assigning more than one level
      thresholds<- data.frame(do.call(rbind, bounds))
      thresholds$classif <- as.factor(thresholds$classif)
     # print("reached 82")
      if(all(classif)) {
        thresholds$y<- do.call(rbind, lapply(preplot.data.list, function(mydf) max(mydf[,"count"])/likelihoodcut))
      } else {
        thresholds$y<- do.call(rbind, lapply(preplot.data.list, function(mydf) max(mydf[,"count"])/likelihoodcut))[-1]

      }

    } else { # if assigning more than one level
      thresholds<- data.frame(do.call(rbind, bounds))
      #print(str(thresholds))
      thresholds$classif <- as.factor(thresholds$classif)
      if(all(classif)) { #if no uncertain cells
        thresholds$y<- unlist(lapply(preplot.data.list, function(mydf) max(mydf[,"count"])/likelihoodcut))
      } else {
        thresholds$y<- unlist(lapply(preplot.data.list, function(mydf) max(mydf[,"count"])/likelihoodcut))[-1]
      }
    }
    #print("reached 93")
    myplot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = vec, fill = classif)) +
      ggplot2::geom_histogram(bins = 400) +
      ggplot2::geom_segment(data = thresholds, ggplot2::aes(x = min, xend = max,
                                          y = y, yend = y),
                   col = "black", inherit.aes = FALSE) +
      ggplot2::scale_x_continuous(trans = asinh_trans(150)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::theme_classic()

  }

  #print(str(classif))
  if(output == "plot"){
    return(myplot)
  } else if (output == "both"){
    return(list(classif = classif,
                plot = myplot))
  }
  # hist(apply(probs.norm.row,1, max), n = 21)
}
