#' Extracts assignments from fcbFlowFrame object
#'
#' @param fcbFlowFrame a fcbFlowFrame object that has been assigned using assign_fcbFlowFrame
#' @param platemap data.frame, a lookup table from barcoding levels to well assignments, (optional)
#' @param simplify logical, return factor instead of list
#' @return a list of factors, one for each barcoding channel that has been assigned
#' @export
getAssignments <- function(x, platemap = NULL, simplify = FALSE) {
  getAssignments.ff <- function(x) {
    assignments <- lapply(x@barcodes, `[[`, "assignment")
    x@barcodes$pacific_orange_a$deskewing
    assignments <- lapply(assignments, `[[`, "values")
    assignments <- lapply(assignments, as.factor)
  }

  if (class(x) == "fcbFlowFrame") {
    assignments <- getAssignments.ff(x)
  } else if (class(x) == "fcbFlowSet") {
    assignments <- fsApply(x, getAssignments.ff)
  }
  if (!is.null(platemap)) {
    #return wells instead of levels
  }

  if (simplify) {
    #some sort of unlist operation
  }
  return(assignments)
}
