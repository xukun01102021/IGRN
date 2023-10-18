#Loading R packages
packages <- c("Seurat","IGRN","Matrix","parallel","scMLnet","SCENIC",
  "Seurat","SCopeLoomR","AUCell","dplyr","KernSmooth","RColorBrewer",
  "plotly","BiocParallel","grid","ComplexHeatmap","data.table","scRNAseq","tidyr"
)

for (package in packages) {
  library(package, character.only = TRUE)
}

#Constructing the initial IGRN
setwd("D:/IGRN/")
file <- readRDS('D:/IGRN/tumor.rds')

runsMLnet(file = file,pval = 0.05,logfc = 0.15,
          LigClu = 'Macrophage',RecClu = 'm_cluster5',
          workdir = './output')

#Analysis of transcription factors
loomFile <- 'D:/IGRN/scenic/m_cluster8/sample_SCENIC.loom'
runSCENIC(loomFile)

#Transcription factor and target gene screening
#The final file is saved under GRN Files
scenicTF <- 'D:/IGRN/scenic_TF.csv'
net_TfTar <- 'D:/IGRN/output/output/Macrophage_m_cluster5/TFTarGene.net.txt'
net_RecTf <- 'D:/IGRN/output/output/Macrophage_m_cluster5/RecTF.net.txt'
net_ligRec <- 'D:/IGRN/output/output/Macrophage_m_cluster5/LigRec.net.txt'
DEGs <- 'D:/Ryy/tumor/DEGs.markers.csv'
outputDir <- 'GRN'
GRNetwork(scenicTF, net_TfTar, net_RecTf, net_ligRec, DEGs,outputDir)

