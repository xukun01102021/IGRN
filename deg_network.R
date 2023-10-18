#' Title
#'
#' @param cancer_count_file
#' @param cancer_data_file
#' @param mac_count_file
#' @param mac_data_file
#'
#' @return
#' @export
#'
#' @examples
deg_network <- function(cancer_count_file, cancer_data_file, mac_count_file, mac_data_file) {
  cancer_count_deg <- read.csv(cancer_count_file)
  cancer_data_deg <- read.csv(cancer_data_file)
  mac_count_deg <- read.csv(mac_count_file)
  mac_data_deg <- read.csv(mac_data_file)


  ligand_count_deg <- cancer_count_deg[cancer_count_deg$X %in% LR$from,]
  ligand_count_deg <- ligand_count_deg[ligand_count_deg$p_val_adj < 0.05,]
  ligand_count_deg <- ligand_count_deg[abs(ligand_count_deg$avg_log2FC) > 0.25,]

  ligand_data_deg <- cancer_data_deg[cancer_data_deg$X %in% LR$from,]
  ligand_data_deg <- ligand_data_deg[ligand_data_deg$p_val_adj < 0.05,]
  ligand_data_deg <- ligand_data_deg[abs(ligand_data_deg$avg_log2FC) > 0.25,]

  ligand_deg <- rbind(ligand_count_deg, ligand_data_deg)

  receptor_count_deg <- mac_count_deg[mac_count_deg$X %in% LR$to,]
  receptor_count_deg <- receptor_count_deg[receptor_count_deg$p_val_adj < 0.05,]
  receptor_count_deg <- receptor_count_deg[abs(receptor_count_deg$avg_log2FC) > 0.25,]

  receptor_data_deg <- mac_data_deg[mac_data_deg$X %in% LR$to,]
  receptor_data_deg <- receptor_data_deg[receptor_data_deg$p_val_adj < 0.05,]
  receptor_data_deg <- receptor_data_deg[abs(receptor_data_deg$avg_log2FC) > 0.25,]

  receptor_deg <- rbind(receptor_count_deg, receptor_data_deg)

  TF_count_deg <- mac_count_deg[mac_count_deg$X %in% RT$to,]
  TF_count_deg <- TF_count_deg[TF_count_deg$p_val_adj < 0.05,]
  TF_count_deg <- TF_count_deg[abs(TF_count_deg$avg_log2FC) > 0.25,]

  TF_data_deg <- mac_data_deg[mac_data_deg$X %in% RT$to,]
  TF_data_deg <- TF_data_deg[TF_data_deg$p_val_adj < 0.05,]
  TF_data_deg <- TF_data_deg[abs(TF_data_deg$avg_log2FC) > 0.25,]

  TF_deg <- rbind(TF_count_deg, TF_data_deg)

  target_count_deg <- mac_count_deg[mac_count_deg$X %in% TT$to,]
  target_count_deg <- target_count_deg[target_count_deg$p_val_adj < 0.05,]
  target_count_deg <- target_count_deg[abs(target_count_deg$avg_log2FC) > 0.25,]

  target_data_deg <- mac_data_deg[mac_data_deg$X %in% TT$to,]
  target_data_deg <- target_data_deg[target_data_deg$p_val_adj < 0.05,]
  target_data_deg <- target_data_deg[abs(target_data_deg$avg_log2FC) > 0.25,]

  target_deg <- rbind(target_count_deg,target_data_deg)

  LR_count_deg <- LR[LR$from %in% ligand_count_deg$X,]
  LR_data_deg <- LR[LR$from %in% ligand_data_deg$X,]
  L_Rdeg <- LR[LR$to %in% receptor_deg$X,]
  LR_degs <- rbind(LR_count_deg,LR_data_deg,L_Rdeg)
  LR_degs <- data.frame(unique(LR_degs$from_to))
  colnames(LR_degs) <- "L_R"
  LR_degs <- tidyr::separate(LR_degs,L_R,into = c("from","to"),sep = "_")
  write.csv(LR_degs,file = "LR_deg.csv")

  RT <- RT[RT$from %in% LR_degs$to,]
  Rsure_T  <- RT[RT$from %in% receptor_deg$X,]
  R_Tsure <- RT[RT$to %in% TF_deg$X,]
  Rsure_Tsure <- rbind(Rsure_T,R_Tsure)
  write.csv(Rsure_Tsure,"RT_deg.csv")

  TT <- TT[TT$from %in% Rsure_Tsure$to,]
  Tsure_Tsure <- TT[TT$from %in% R_Tsure$to,]
  T_Tsure <- TT[TT$from %in% Rsure_T$to,]
  write.csv(T_Tsure,"TT_deg.csv")}
