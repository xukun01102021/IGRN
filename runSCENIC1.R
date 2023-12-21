#' Title
#'
#' @param loomFile
#'
#' @return
#' @export
#'
#' @examples
runSCENIC <- function(loomFile) {
  loom <- open_loom(loomFile)
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
  regulonAucThresholds <- get_regulon_thresholds(loom)
  regulons <- regulonsToGeneLists(regulons_incidMat)
  scenic_TF <- as.data.frame(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
  colnames(scenic_TF) <- "TF"
  write.csv(scenic_TF, "scenic_TF.csv", row.names = FALSE)
}
