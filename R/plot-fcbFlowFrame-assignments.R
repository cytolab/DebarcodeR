#' Plots assignemnts
#'
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing, clustering, and assignment
#' @param plot character, which plot to generate so far, only assignemnts have been implemented
#' @return a fcbFlowFrame object with a barcode slot filled with deskewing, clustering, cell assignment as
#' a vector of integers from 0:ncol(probs), cells assigned a classification of 0 remained unassigned,
#' otherwise number corresponds to the barcoding level assignment of that cell
#' @export
#' @import ggplot2 ggnewscale scales ggforce
plot.fcbflowframe <- function(fcbFlowFrame, plot = "assignments", seed = 1) {

  if (class(fcbFlowFrame) != "fcbFlowFrame") {
    stop("Input must be a fcbFlowFrame")
  }
  barcodes <- fcbFlowFrame@barcodes[!(names(fcbFlowFrame@barcodes) == "wells")]
  deskewed <- as.data.frame(lapply(barcodes, function(bc) bc[['deskewing']][['values']]))
  assignment <- as.data.frame(lapply(barcodes, function(bc) bc[['assignment']][['values']]))
  assignment <- apply(assignment, 1, paste0, collapse = ".")
  assignment[grepl(0, assignment)] <- "0.0"
  assignment[grepl("0", assignment)] <- "0.0"
  set.seed(seed)
  mypal <- hue_pal()(length(table(assignment)) - 1)
  mypal <- c("grey50", sample(mypal))
  mydata <- cbind(deskewed, assignment = assignment)
  mydata.assigned <- mydata[mydata$assignment != "0.0",]

  if (plot == 'data') {
    myplot <- mydata.assigned
  }
  if (plot == "assignments") {

    myplot <- ggplot(mydata,
      aes_string(
        y = colnames(deskewed)[1],
        x = colnames(deskewed)[2],
        col  = 'assignment'
      )
    ) +
      geom_point(shape = ".") +
      scale_color_manual(values = mypal)
  } else if (plot == 'density') {
    myplot <- ggplot(
      mydata,
      aes_string(
        y = colnames(deskewed)[1],
        x = colnames(deskewed)[2]
      )
    ) +
      stat_pointdensity(geom = ggrastr:::GeomPointRast,
                        adjust = 0.1,
                        shape = ".",
                        n = 256,
                        method = "kde2d",
                        raster.width = 4, raster.height = 4,
      )  +
      scale_color_viridis_c(option = "A", guide = FALSE) +
      ggnewscale::new_scale_color() +
      stat_density_2d(aes(fill = assignment, alpha = after_stat(level)),
                      data = mydata.assigned,
                      geom = "polygon", color = NA, bins = 8) +
      scale_fill_manual(values = mypal[-1])


      myplot +
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.position = "none")
      ggsave('density-cl.png', width = 4, height = 4, units = 'in', dpi = 300)
  } else if (plot == 'chull') {
    mydata.chull <- as_tibble(mydata.assigned) %>%
      group_by(assignment) %>%
      slice(chull(pacific_orange_a, pacific_blue_a))

    myplot <- ggplot(
      mydata,
      aes_string(
        y = colnames(deskewed)[1],
        x = colnames(deskewed)[2]
      )
    ) +
    stat_pointdensity(geom = ggrastr:::GeomPointRast,
                      adjust = 0.1,
                      shape = ".",
                      n = 256,
                      method = "kde2d",
                      raster.width = 4, raster.height = 4,
    )  +
    scale_color_viridis_c(option = "A", guide = FALSE) +
    ggnewscale::new_scale_color() +
    # geom_shape(data = mydata.chull, aes(color = assignment),
    #      expand = unit(0.25, "mm"),
    #      radius = unit(1, 'mm'), fill = NA) +
      geom_shape(data = mydata.chull, aes(fill = assignment),
                 expand = unit(0.25, "mm"),
                 radius = unit(1, 'mm'), alpha = 0.25) +
      # geom_mark_hull(aes(fill = assignment), data = mydata.assigned, alpha = 0.3,
      #                concavity = 0,
      #                expand  = unit(1, "mm"), radius = unit(2, "mm")) +
      #geom_voronoi_tile(aes(fill = assignment), data = mydata.assigned) +
      scale_fill_manual(values = mypal[-1])

    # myplot.polished  <- myplot +
    #   theme_bw() +
    #   theme(panel.grid = element_blank(),
    #         legend.position = "none")

    #ggsave('density-rounded-chull-filled.png', myplot.polished, width = 4, height = 4, units = 'in', dpi = 300)
  }

  return(myplot)
  }
