#' Title
#'
#' @param ligands
#' @param lr
#' @param rl
#'
#' @return
#' @export
#'
#' @examples
cellphone_data <- function(ligands,lr,rl) {
  meansdf11 <- mymeans %>%
    dplyr::filter(interacting_pair %in% ligands) %>%
    dplyr::select("interacting_pair", starts_with(lr)) %>%
    reshape2::melt()
  meansdf1 <- mymeans %>%
    dplyr::filter(interacting_pair %in% ligands) %>%
    dplyr::select("interacting_pair", ends_with(rl)) %>%
    reshape2::melt()
  meansdf <- rbind(meansdf1, meansdf11)
  rm(meansdf1, meansdf11)
  colnames(meansdf) <- c("interacting_pair", "CC", "means")

  pvalsdf11 <- mypvals %>%
    dplyr::filter(interacting_pair %in% ligands) %>%
    dplyr::select("interacting_pair", starts_with(lr)) %>%
    reshape2::melt()
  pvalsdf1 <- mypvals %>%
    dplyr::filter(interacting_pair %in% ligands) %>%
    dplyr::select("interacting_pair", ends_with(rl)) %>%
    reshape2::melt()
  pvalsdf <- rbind(pvalsdf1, pvalsdf11)
  rm(pvalsdf1, pvalsdf11)
  colnames(pvalsdf) <- c("interacting_pair", "CC", "pvals")
  pvalsdf$joinlab <- paste0(pvalsdf$interacting_pair, "_", pvalsdf$CC)
  meansdf$joinlab <- paste0(meansdf$interacting_pair, "_", meansdf$CC)

  pldf <- merge(pvalsdf, meansdf, by = "joinlab")
  pldf <- pldf[pldf$means > 0,]
  write.csv(pldf, file = "LR.csv")

  tlung_significant_means_x <- significant_means[significant_means$interacting_pair %in% pldf$interacting_pair.x,]
  tlung_significant_means_y <- significant_means[significant_means$interacting_pair %in% pldf$interacting_pair.y,]
  tlung_significant_means_xy <- rbind(tlung_significant_means_x, tlung_significant_means_y)
  write.csv(data.frame(unique(tlung_significant_means_xy$interacting_pair)), "significant_means.csv", row.names = FALSE)
  write.csv(tlung_significant_means_xy, "significant_means.csv")
}
