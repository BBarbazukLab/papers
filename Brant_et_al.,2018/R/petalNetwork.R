# Install or Update petal
# library(devtools)
# install_github("julipetal/petalNet", subdir="src/petal")

library(petal)
library(igraph)

setwd("/Users/jasonbrant/Desktop/petal_data/ACO_VS_MUS_02-10-17")

graphHistQQ(expM, fileName = "Aco_Mus_QQ")
#write.table(expM, file="acomys_genes_counts_matrix_filtered.txt",quote = F, sep = '\t')

createSWSFnetFromFile(expMFile = "/Users/jasonbrant/Desktop/petal_data/ACO_VS_MUS_02-10-17/aco_vs_mus_counts.txt", metric = "SP")


# downstreamAnalysis(winnerT = "0.914",
#                    metric = "SP",
#                    outFile = "corematrisome_collagens_mm.txt_VN.txt",
#                    expMatrixFile = "/Users/jasonbrant/Desktop/petal_data/ACO_VS_MUS_02-10-17/aco_vs_mus_counts.txt",
#                    GoIFile = "../corematrisome_collagens_mm.txt")
# 
makeCytoFile(threshold = "0.914", "SP", orderedMM = MMo)

