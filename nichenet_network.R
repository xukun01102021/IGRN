#' Title
#'
#' @param ligand_target_matrix
#'
#' @return
#' @export
#'
#' @examples
nichenet_network <- function(ligand_target_matrix) {
  ligands_targets <- read.csv(ligand_target_matrix)

  names <- setdiff(colnames(ligands_targets), c("X"))
  ligands_targets <- melt(ligands_targets, id.vars = "X", measure.vars = names,
                          variable.names = "target", value.name = "expr")
  ligands_targets <- ligands_targets[ligands_targets$expr > 0,]
  write.csv(ligands_targets, "ligands_targets_score.csv")

  ligands_all <- unique(ligands_targets$X)
  targets_all <- as.character(unique(ligands_targets$variable))

  ligands_all <- ligands_all[which(ligands_all %in% colnames(ligand_tf_matrix))]
  targets_all <- targets_all[which(targets_all %in% gr_network$to)]

  active_signaling_network <- get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix,
                                                        ligands_all = ligands_all,
                                                        targets_all = targets_all,
                                                        weighted_networks = weighted_networks)

  active_signaling_network_min_max <- active_signaling_network

  graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max,
                                                     ligands_all = ligands_all,
                                                     targets_all = targets_all,
                                                     sig_color = "indianred",
                                                     gr_color = "steelblue")

  data_source_network <- infer_supporting_datasources(signaling_graph_list = active_signaling_network,
                                                      lr_network = lr_network,
                                                      sig_network = sig_network,
                                                      gr_network = gr_network)

  x1 <- bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"),
                  active_signaling_network$gr %>% mutate(layer = "regulatory"))
  write.csv(x1, "weighted_signaling_network.csv")

  x2 <- data_source_network
  write.csv(x2, "data_source_network.csv")

  specific_annotation_tbl = bind_rows(
    tibble(gene = ligands_all, annotation = "ligand"),
    tibble(gene = targets_all, annotation = "target"),
    tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
    tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
  )
  non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

  x3 <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
  write.csv(x3,"annotation_table.csv")

  ligands <- x3[which(x3$annotation %in% "ligand"),]
  receptor <- x3[which(x3$annotation %in% "receptor"),]
  TFs <- x3[which(x3$annotation %in% "transcriptional regulator"),]
  targets <- x3[which(x3$annotation %in% "target"),]
  write.csv(ligands,"ligands.csv")
  write.csv(receptor,"receptor.csv")
  write.csv(TFs,"TFs.csv")
  write.csv(targets,"targets.csv")

  LR <- x1[x1$from %in% ligands$gene,]
  LR <- LR[LR$to %in% receptor$gene,]
  write.csv(LR,"LR_nichenet.csv")
  RT <- x1[x1$from %in% LR$to,]
  RT <- RT[RT$to %in% TFs$gene,]
  write.csv(RT,"RT_nichenet.csv")
  TFtarget <- x1[x1$from %in% RT$to,]
  TFtarget <- TFtarget[TFtarget$to %in% targets$gene,]
  write.csv(TFtarget,"TFtarget_nichenet.csv")
}









