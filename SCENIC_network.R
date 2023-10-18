#' Title
#'
#' @param scenicLoomPath
#'
#' @return
#' @export
#'
#' @examples
SCENIC_network <- function(scenicLoomPath) {
  loom <- open_loom(scenicLoomPath, mode = "r")
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  regulonsAucThresholds <- get_regulon_thresholds(loom)

  TFs <- character()
  for (i in unique(TT$to)) {
    for (j in names(regulons)) {
      if (i %in% regulons[[j]]) {
        TFs <- c(TFs, j)
      } else {
        print(1)
      }
    }
  }

  RT$sign <- "(+)"
  RT <- unite(RT, "to_sign", to, sign, sep = "", remove = FALSE)
  TT$sign <- "(+)"
  TT <- unite(TT, "from_sign", from, sign, sep = "", remove = FALSE)
  TT <- unite(TT, "from_sign_to", from_sign, to, sep = "_", remove = FALSE)

  TF <- names(table(TFs))
  TF <- TF[TF %in% TT$from_sign]
  TFs_regulons <- regulons[names(regulons) %in% TF]
  xx1 <- list()
  for (i in names(TFs_regulons)) {
    xx1[[i]] <- data.frame(i, TFs_regulons[[i]])
  }
  names(xx1)
  yyy1 <- combine(xx1)
  colnames(yyy1) <- c("x", "y")
  yyy1 <- unite(yyy1, "x_y", x, y, sep = "_", remove = FALSE)

  TT_scenic <- TT[TT$from_sign_to %in% yyy1$x_y,]
  write.csv(TT_scenic, "TT_scenic.csv")

  RT_scenic <- RT[RT$to_sign %in% TT_scenic$from_sign,]
  write.csv(RT_scenic, "RT_scenic.csv")

  LR_scenic <- LR[LR$to %in% RT_scenic$from,]
  write.csv(LR_scenic, "LR_scenic.csv")
}
