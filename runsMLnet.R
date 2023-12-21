#' Title
#'
#' @param file
#' @param pval
#' @param logfc
#' @param LigClu
#' @param RecClu
#' @param workdir
#'
#' @return
#' @export
#'
#' @examples
runsMLnet <- function(file, pval, logfc, LigClu, RecClu, workdir) {
  afcell <- as.data.frame(file@active.ident)
  afcell <- afcell[rownames(file@meta.data),,drop=F]
  file$cell <- afcell[,1]

  BarCluFile <- as.data.frame(file@meta.data[,"cell",drop=F])
  colnames(BarCluFile) <- "Cluster"
  out <- cbind(Barcode=rownames(BarCluFile), BarCluFile)
  write.table(out, file="cellType.txt", sep="\t", quote=F, row.names=F)

  netList <- RunMLnet(file@assays$RNA@counts,
                      "cellType.txt", RecClu, LigClu, pval, logfc)
  DrawMLnet(netList, LigClu, RecClu, workdir)
}
