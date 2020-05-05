#' Splits FCB flowFrame by platemap
#'
#' @param fcbFlowFrame a fcbFlowFrame object with barcoded flowframe and uptake flowframe post deskewing
#' and clustering (at least one barcodes slot filled)

#' @return a fcbFlowSet with each flowFrame named by the assignments defined in the platemap
#' @export
assign_fcbFlowFrame <- function(fcbFlowFrame,
                                simplify = TRUE){

  cell_lut <- lapply(fcbFlowFrame@barcodes, function(x) x$assignment$values) %>%
    as_tibble() %>%
    left_join(fcbFlowFrame@platemap %>% janitor::clean_names()) %>%
    mutate(well = if_else(is.na(well), "Unassigned", well))


  fcbFlowFrame.list <- split(fcbFlowFrame@barcoded.ff, cell_lut$well)

  orig.desc <- description(fcbFlowFrame@barcoded.ff)
  orig.filename <- basename(orig.desc$FILENAME)
  split.names <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", orig.filename),
         "_",
         names(fcbFlowFrame.list))
  names(fcbFlowFrame.list) <- split.names

  if (simplify == TRUE) {
    flowSet.x <- flowSet(fcbFlowFrame.list)
    return(flowSet.x)
  } else {

  }
}
