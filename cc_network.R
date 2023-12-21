#Loading R packages
ps <- c('ccnetwork', 'nichenetr','tidyverse','foreach','Seurat','SCENIC','doParallel','pheatmap',
        'SCopeLoomR','AUCell','hdf5r','tidydr','dplyr','data.table','reshape','reshape2','tidyverse','tidyr','openxlsx','magrittr')
lapply(ps, function(x){library(x, character.only = T)})

#step1 Construct network (version=1)
#Documents required for loading nichenet, and 'ligand_target_matrix' analysed by the Institute.
setwd("./nichenet/")
weighted_networks <- readRDS("weighted_networks.rds")
ligand_tf_matrix <- readRDS("ligand_tf_matrix.rds")
lr_network <- readRDS("lr_network.rds")
signaling_network <- readRDS("signaling_network.rds")
gr_network <- readRDS("gr_network.rds")
ligands_targets <- "ligand_target_matrix.csv"

nichenet_network(ligands_targets)

#step2 Construct for the network (version=2) by cellphoneDB
setwd("./cellphonedb/")
#Loading networkV1 and cellphoneDB results
LR <- read.csv("./nichenet/LR_nichenet.csv",row.names = 1)
RT <- read.csv("./nichenet/RT_nichenet.csv",row.names = 1)
TFtarget <- read.csv("./nichenet/TFtarget_nichenet.csv",row.names = 1)
mypvals <- read.delim("pvalues.txt", check.names = FALSE)
mymeans <- read.delim("means.txt", check.names = FALSE)
significant_means <- read.delim("significant_means.txt")

ligands <- grep("CELSR1|COL18A1|CTGF|CXCL14|CXCL16|", #Ligand signalling in network (version=1)
                mymeans$interacting_pair,value = T)

cellphone_data(ligands,lr = "m_cluster1|Macrophage",rl = "Macrophage|m_cluster1")
LR_cellphoneDB <- read.csv("LR.csv")
cellphone_network(LR_cellphoneDB,Lr = "m_cluster1.Macrophage",Rl = "Macrophage.m_cluster1")

#step3 Construct for the network (version=3) by pySCENIC
setwd("./scenic/")
#Loading networkV2 and pySCENIC results
LR <- read.csv("./cellphonedb/LR_verified.csv",row.names = 1)
RT <- read.csv("./cellphonedb/RT_verified.csv",row.names = 1)
TT <- read.csv("./cellphonedb/TT_verified.csv",row.names = 1)
scenicLoomPath <- "./scenic/sample_SCENIC.loom"

SCENIC_network(scenicLoomPath)

#step4 Construct for the network (version=3) by DEGs
setwd("./degs/")
#Loading networkV3 and DEGs results
cancer_count_deg <- "./degs/counts_cancer_tlungnlung_markers.csv"
cancer_data_deg <- "./degs/data_cancer_tlungnlung_markers.csv"
mac_count_deg <- "./degs/counts_mac_tlungnlung_markers.csv"
mac_data_deg <- "./degs/data_mac_tlungnlung_markers.csv"

LR <- read.csv("./scenic/LR_scenic_tlung.csv",row.names = 1)
RT <- read.csv("./scenic/RT_scenic_tlung.csv")
TT <- read.csv("./scenic/TT_scenic_tlung.csv")
#Getting the final regulatory network
deg_network(cancer_count_deg,cancer_data_deg,mac_count_deg,mac_data_deg)




