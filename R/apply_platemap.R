#' Updates pData and filenames with well assignments from platemap
#'
#' @param fcbFlowSet debarcoded flowframe(s)
#' @param platemap data.frame, with the well assignments
#' @param drop0 logical, drop unassigned events, default = FALSE
#' @param prefix character, file prefix for new FCS files, if not specifed, will pull from FCS file $FIL keyword
#' @return an fcbFlowSet object with updated names and pData
#' @import flowCore janitor
#' @importFrom dplyr %>% left_join tibble mutate mutate_all
#' @export
apply_platemap <- function(fcbFlowSet, platemap, drop0 = FALSE, prefix = NA) {
  if (class(fcbFlowSet) != "fcbFlowSet") {
    stop("Input must be a fcbFlowSet")
  }
  barcoded.ff@description$FILENAME
  if (!any(class(platemap) == "data.frame")) {
    stop("Input must be a data.frame")
  }
  pData.orig <- pData(fcbFlowSet)
  suppressMessages({
    pData.new <- pData.orig %>%
      left_join(platemap %>%
                  janitor::clean_names() %>%
                  dplyr::mutate_all(as.character))
  })
  #rownames(pData.new) <- row.names(pData.orig)
  assigned.fs <- fcbFlowSet[!is.na(pData.new$well)]
  if (drop0 == FALSE) {
    unassigned.fs <- fcbFlowSet[is.na(pData.new$well)]
    unassigned.list <- flowSet_to_list(unassigned.fs)
    unassigned.concat <- do.call(rbind, lapply(unassigned.list, exprs))
    unassigned.ff <- unassigned.list[[1]]
    exprs(unassigned.ff) <- unassigned.concat
    out.list <- flowSet_to_list(assigned.fs)
    out.list <- c(out.list, "Unassigned" = unassigned.ff)
    out.fs <- flowSet(out.list)
    suppressMessages({suppressWarnings({
      pData.new <- left_join(pData(out.fs), pData(fcbFlowSet)) %>%
        left_join(platemap %>%
                    janitor::clean_names() %>%
                    mutate_all(as.character)) %>%
        mutate(well = if_else(is.na(well), "Unassigned", well))
    })})
    rownames(pData.new) <- rownames(pData(out.fs))
    pData(out.fs) <- pData.new
  } else {
    out.fs <- assigned.fs
  }
  if (is.na(prefix)) {
    pData(out.fs)$Prefix <- tools::file_path_sans_ext(unlist(fsApply(out.fs ,keyword, "FILENAME")))
  } else ({
    pData(out.fs)$Prefix <- prefix
  })
    pData(out.fs)$Filename <- apply(pData(out.fs)[, c("Prefix", "well")], 1, paste0, collapse = "_")
    sampleNames(out.fs) <- pData(out.fs)$Filename
    return(out.fs)
}
