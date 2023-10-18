#' Title
#'
#' @param LR_cellphoneDB
#' @param Lr
#' @param Rl
#'
#' @return
#' @export
#'
#' @examples
cellphone_network <- function(LR_cellphoneDB,Lr,Rl) {

  LR <- tidyr::unite(LR, "from_to", from, to, sep = "_", remove = FALSE)
  LR <- tidyr::unite(LR, "to_from", to, from, sep = "_", remove = FALSE)
  LR_verified_from_to <- LR[LR$from_to %in% LR_cellphoneDB$interacting_pair.x, ]
  LR_verified_to_from <- LR[LR$to_from %in% LR_cellphoneDB$interacting_pair.x, ]
  LR_verified <- rbind(LR_verified_from_to, LR_verified_to_from)
  write.csv(LR_verified, "LR_verified.csv")

  significant_means <- significant_means[, colnames(significant_means) %in% c("id_cp_interaction", "interacting_pair", "partner_a",
                                                                              "partner_b", "gene_a", "gene_b",
                                                                              "secreted", "receptor_a", "receptor_b",
                                                                              "annotation_strategy", "is_integrin", "rank",
                                                                              Lr, Rl)]
  lr1 <- subset(significant_means, interacting_pair %in% LR_verified$from_to)
  lr2 <- subset(significant_means, interacting_pair %in% LR_verified$to_from)
  lr <- rbind(lr1, lr2)
  write.csv(lr, 'LR_cellphone.csv')

  RT <- RT[RT$from %in% LR_verified$to, ]
  write.csv(RT, "RT_verified.csv")

  TT <- TFtarget[TFtarget$from %in% RT$to, ]
  write.csv(TT, "TT_verified.csv")
}
