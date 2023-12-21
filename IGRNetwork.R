#' Title
#'
#' @param scenicTF
#' @param net_TfTar
#' @param net_RecTf
#' @param net_ligRec
#' @param DEGs
#' @param outputDir
#'
#' @return
#' @export
#'
#' @examples
GRNetwork <- function(scenicTF, net_TfTar, net_RecTf, net_ligRec, DEGs,outputDir) {
  scenic_TF <- read.csv(scenicTF)
  net_TF_Tar <- read.table(net_TfTar, sep = "\t")
  net_TF_Tar <- separate(net_TF_Tar, col = "V1", into = c("TF", "Target"), sep = "_", remove = FALSE)
  net_TF_Tar <- net_TF_Tar[, -1]
  net_TF_Tar$singn <- "(+)"
  net_TF_Tar <- tidyr::unite(net_TF_Tar, col = "TF", c("TF", "singn"), sep = "")
  net_TF_Tar <- subset(net_TF_Tar, TF %in% c(intersect(scenic_TF$TF, net_TF_Tar$TF)))

  DEGs <- read.csv(DEGs)
  net_TF_Tar <- subset(net_TF_Tar, Target %in% c(intersect(net_TF_Tar$Target, DEGs$X)))

  write.csv(net_TF_Tar, file.path(outputDir, "sub_TF_Target.csv"), row.names = FALSE)

  net_Rec_TF <- read.table(net_RecTf, sep = "\t")
  net_Rec_TF <- separate(net_Rec_TF, col = "V1", into = c("Rec", "TF"), sep = "_", remove = FALSE)
  net_Rec_TF <- net_Rec_TF[, -1]
  net_Rec_TF$singn <- "(+)"
  net_Rec_TF <- tidyr::unite(net_Rec_TF, col = "TF", c("TF", "singn"), sep = "")
  net_Rec_TF <- subset(net_Rec_TF, TF %in% c(net_TF_Tar$TF))

  write.csv(net_Rec_TF, file.path(outputDir, "sub_Rec_TF.csv"), row.names = FALSE)

  net_Lig_Rec <- read.table(net_ligRec, sep = "\t")
  net_Lig_Rec <- separate(net_Lig_Rec, col = "V1", into = c("Lig", "Rec"), sep = "_", remove = FALSE)
  net_Lig_Rec <- net_Lig_Rec[, -1]
  net_Lig_Rec <- subset(net_Lig_Rec, Rec %in% c(net_Rec_TF$Rec))

  write.csv(net_Lig_Rec, file.path(outputDir, "sub_Lig_Rec.csv"), row.names = FALSE)
}
