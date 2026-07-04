#' ---
#' title: "Root cold stress network analysis"
#' subtitle: ""
#' author: "Elena van Zalen & Tuuli Aro & Alex Vergara"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    fig_width: 9
#'    fig_height: 6
#'    toc: true
#'    number_sections: true
#'    toc_depth: 3
#'    toc_float:
#'      collapsed: TRUE
#'      smooth_scroll: TRUE
#'    code_folding: hide
#'    theme: "flatly"
#'    highlight: pygments
#'    includes:
#'      before_body: header.html
#'      after_body: footer.html
#'    css: style.css
#' ---
#' 
#' Loading libraries
library(tidyverse)
theme_set(theme_classic())
theme_update(plot.title = element_text(face="bold"))
suppressPackageStartupMessages({
  library(arrow)
  library(cowplot)
  library(DESeq2)
  library(DT)
  library(gplots)
  library(igraph)
  library(matrixStats)
  library(purrr)
  library(RColorBrewer)
})


# Species name in the orthologs.paquet file
# "Arabidopsis_thaliana", "Betula_pendula", Picea_abies", "Pinus_sylvestris", "Populus_tremula"
species1_name <- "col0"
species2_name <- "ost0"
species3_name <- "aspen"
species4_name <- "birch"
species5_name <- "spruce"
species6_name <- "pine"

species1_expr_file <- "data/expression/col0_root_cold_stress_expression.txt"
species2_expr_file <- "data/expression/ost0_root_cold_stress_expression.txt"
species3_expr_file <- "data/expression/aspen_root_cold_stress_expression.txt"
species4_expr_file <- "data/expression/birch_root_cold_stress_expression.txt"
species5_expr_file <- "data/expression/spruce_root_cold_stress_expression.txt"
species6_expr_file <- "data/expression/pine_root_cold_stress_expression.txt"

dir.create("data/ComPlEx/networks", recursive = TRUE, showWarnings = FALSE)

# Read in ortholog groups
# =======================
# ortholog_group_file <- read_parquet("data/orthologs.parquet")
# 
# # Changing the species names of our species of interest for easy of interpretation later on
# species_rename_vector <- c("Arabidopsis_thaliana" = "col0", "Populus_tremula" = "aspen", 
#                            "Betula_pendula" = "birch", "Picea_abies" = "spruce", "Pinus_sylvestris" = "pine")
# 
# ortholog_group_RData <- paste0("data/ComPlEx/orthologs-ALL_SIX.RData")
# if (!file.exists(ortholog_group_RData)) {
#   ortho <- ortholog_group_file %>% 
#     mutate(species = dplyr::recode(species, !!!species_rename_vector)) %>% 
#     bind_rows(., filter(., grepl("col0", species)) %>%
#                 mutate(species = gsub("col0", "ost0", species)))
#   ortho <- ortho %>%
#     filter(species %in% c(species1_name, species2_name, species3_name,
#                           species4_name, species5_name, species6_name))
#   ortho_A <- ortho %>% 
#     dplyr::select(species,gene,Ortholog_Group) %>% 
#     filter(str_detect(species, species1_name)) 
#   ortho_B <- ortho %>% 
#     dplyr::select(species,gene,Ortholog_Group) %>% 
#     filter(str_detect(species, species2_name))
#   ortho_C <- ortho %>% 
#     dplyr::select(species,gene,Ortholog_Group) %>% 
#     filter(str_detect(species, species3_name))
#   ortho_D <- ortho %>% 
#     dplyr::select(species,gene,Ortholog_Group) %>% 
#     filter(str_detect(species, species4_name))
#   ortho_E <- ortho %>% 
#     dplyr::select(species,gene,Ortholog_Group) %>% 
#     filter(str_detect(species, species5_name))
#   ortho_F <- ortho %>% 
#     dplyr::select(species,gene,Ortholog_Group) %>% 
#     filter(str_detect(species, species6_name))
#   ortho_list <- list(ortho_A, ortho_B, ortho_C, ortho_D, ortho_E, ortho_F)
#   ortho <- ortho_list %>%  #purrr::reduce(list(x,y,z), dplyr::left_join, by = 'Flag')
#     purrr::reduce(list(ortho_A, ortho_B, ortho_C, ortho_D, ortho_E, ortho_F),
#                   dplyr::inner_join, by = "Ortholog_Group") %>% 
#     dplyr::select(-c(species.x, species.y)) %>% 
#     dplyr::select(Species1 = "gene.x", Species2 = "gene.y", OrthoGroup = "Ortholog_Group")
# 
#   save(ortho, file = ortholog_group_RData)
# } else {
#   load(file = ortholog_group_RData)
# }
# 
# 
# duplicated_keys <- ortho_A$Ortholog_Group[duplicated(ortho_A$Ortholog_Group)]
# ortho_A %>% filter(Ortholog_Group %in% duplicated_keys)
# # Add annotations from arabidopsis
# symbols <- read.delim("data/gene_aliases_20140331.txt", sep = "\t") %>%
#   dplyr::rename(Arabidopsis = locus_name, Symbol = symbol, Name = full_name)
# ortholog_annot_file <- "data/Orthogroups_130323_predefined_tree.tsv"

# Read in expression data
# =======================
# species1_expr <- read.delim(species1_expr_file, sep = "\t", header = TRUE)
# species2_expr <- read.delim(species2_expr_file, sep = "\t", header = TRUE)
# species3_expr <- read.delim(species3_expr_file, sep = "\t", header = TRUE)
# species4_expr <- read.delim(species4_expr_file, sep = "\t", header = TRUE)
# species5_expr <- read.delim(species5_expr_file, sep = "\t", header = TRUE)
# species6_expr <- read.delim(species6_expr_file, sep = "\t", header = TRUE)

# Filter
# ======
# Filtering the expression tables where they only contain the genes that have an ortholog in ortho
# skipped this filter for centrality calculation
# cat (length(unique(ortho$OrthoGroup)), " ortholog groups containing:\n",
#      " ", length(unique(ortho$Species1)), " ", species1_name, " genes\n",
#      " ", length(unique(ortho$Species2)), " ", species2_name, " genes\n\n",
#      length(unique(species1_expr$Genes)), " expressed ", species1_name, " genes\n",
#      length(unique(species2_expr$Genes)), " expressed ", species2_name, " genes\n",
#      sep = "")
# 
# ortho <- ortho %>%
#   filter(Species1 %in% species1_expr$Genes & Species2 %in% species2_expr$Genes)
# 
# species1_expr <- species1_expr[species1_expr$Genes %in% ortho$Species1,]
# species2_expr <- species2_expr[species2_expr$Genes %in% ortho$Species2,]
# 
# cat ("After filtering on expressed genes with ortholog:\n",
#      " ", length(unique(ortho$OrthoGroup)), " ortholog groups containing: \n",
#      "  ", length(unique(ortho$Species1)), " ", species1_name, " genes\n",
#      "  ", length(unique(ortho$Species2)), " ", species2_name, " genes\n",
#      sep = "")

# Read in expression without Genes column
species1_expr <- as.matrix(read.delim(species1_expr_file, sep = "\t", header = TRUE)[,-1])
species2_expr <- as.matrix(read.delim(species2_expr_file, sep = "\t", header = TRUE)[,-1])
species3_expr <- as.matrix(read.delim(species3_expr_file, sep = "\t", header = TRUE)[,-1])
species4_expr <- as.matrix(read.delim(species4_expr_file, sep = "\t", header = TRUE)[,-1])
species5_expr <- as.matrix(read.delim(species5_expr_file, sep = "\t", header = TRUE)[,-1])
species6_expr <- as.matrix(read.delim(species6_expr_file, sep = "\t", header = TRUE)[,-1])

# calculate rowsds
species1_expr <- species1_expr[which(rowSds(species1_expr) > 0),]
species2_expr <- species2_expr[which(rowSds(species2_expr) > 0),]
species3_expr <- species3_expr[which(rowSds(species3_expr) > 0),]
species4_expr <- species4_expr[which(rowSds(species4_expr) > 0),]
species5_expr <- species5_expr[which(rowSds(species5_expr) > 0),]
species6_expr <- species6_expr[which(rowSds(species6_expr) > 0),]

# add back Genes column
species1_expr <- data.frame(Genes = rownames(species1_expr), species1_expr, row.names = NULL)
species2_expr <- data.frame(Genes = rownames(species2_expr), species2_expr, row.names = NULL)
species3_expr <- data.frame(Genes = rownames(species3_expr), species3_expr, row.names = NULL)
species4_expr <- data.frame(Genes = rownames(species4_expr), species4_expr, row.names = NULL)
species5_expr <- data.frame(Genes = rownames(species5_expr), species5_expr, row.names = NULL)
species6_expr <- data.frame(Genes = rownames(species6_expr), species6_expr, row.names = NULL)

# Calculating co-expression networks
# ==================================
cor_method <- "pearson" # pearson spearman
cor_sign <- "" # abs
norm_method <- "MR" # CLR MR, you can find more conservation with MR
density_thr <- 0.03 # this defines what the neighbourhood is, you only consider the 3% highest links
randomize <- "" # randomisation, permuting name of genes, set yes or no, to check if what you see from real data is real

networkfile <- paste0("data/ComPlEx/networks/MR_network_", 
                      species1_name, ".RData")

if (!file.exists(networkfile)) {

if (randomize == "rand") {
    species1_expr$Genes <- sample(species1_expr$Genes, nrow(species1_expr), FALSE)
    species2_expr$Genes <- sample(species2_expr$Genes, nrow(species2_expr), FALSE)
  }
  ############################################################
  species1_net <- cor(t(species1_expr[,-1]), method = cor_method)
  dimnames(species1_net) <- list(species1_expr$Genes, species1_expr$Genes)

  species2_net <- cor(t(species2_expr[,-1]), method = cor_method)
  dimnames(species2_net) <- list(species2_expr$Genes, species2_expr$Genes)

  species3_net <- cor(t(species3_expr[,-1]), method = cor_method)
  dimnames(species3_net) <- list(species3_expr$Genes, species3_expr$Genes)

  species4_net <- cor(t(species4_expr[,-1]), method = cor_method)
  dimnames(species4_net) <- list(species4_expr$Genes, species4_expr$Genes)

  species5_net <- cor(t(species5_expr[,-1]), method = cor_method)
  dimnames(species5_net) <- list(species5_expr$Genes, species5_expr$Genes)

  species6_net <- cor(t(species6_expr[,-1]), method = cor_method)
  dimnames(species6_net) <- list(species6_expr$Genes, species6_expr$Genes)
  #############################################################
  if (cor_sign == "abs") {
    species1_net <- abs(species1_net)
    species2_net <- abs(species2_net)
  }
  
  if (norm_method == "CLR") {
    #                                  NOT INCLUDED
    z <- scale(species1_net)
    z[z < 0] <- 0
    species1_net <- sqrt(t(z)**2 + z**2)
    
    z <- scale(species2_net)
    z[z < 0] <- 0
    species2_net <- sqrt(t(z)**2 + z**2)
    ############################################################
  } else if (norm_method == "MR") {
    R <- t(apply(species1_net, 1, rank))
    species1_net <- sqrt(R * t(R))
    
    R <- t(apply(species2_net, 1, rank))
    species2_net <- sqrt(R * t(R))

    R <- t(apply(species3_net, 1, rank))
    species3_net <- sqrt(R * t(R))

    R <- t(apply(species4_net, 1, rank))
    species4_net <- sqrt(R * t(R))

    R <- t(apply(species5_net, 1, rank))
    species5_net <- sqrt(R * t(R))

    R <- t(apply(species6_net, 1, rank))
    species6_net <- sqrt(R * t(R))
  }############################################################
  
  diag(species1_net) <- 0
  diag(species2_net) <- 0
  diag(species3_net) <- 0
  diag(species4_net) <- 0
  diag(species5_net) <- 0
  diag(species6_net) <- 0
  
  R <- sort(species1_net[upper.tri(species1_net, diag = FALSE)], decreasing = TRUE)
  species1_thr <- R[round(density_thr*length(R))]
  #plot(density(R), xlab = paste0(species1_name, " correlations"), main = "")

  R <- sort(species2_net[upper.tri(species2_net, diag = FALSE)], decreasing = TRUE)
  species2_thr <- R[round(density_thr*length(R))]
  #plot(density(R), xlab = paste0(species2_name, " correlations"), main = "")

  R <- sort(species3_net[upper.tri(species3_net, diag = FALSE)], decreasing = TRUE)
  species3_thr <- R[round(density_thr*length(R))]

  R <- sort(species4_net[upper.tri(species4_net, diag = FALSE)], decreasing = TRUE)
  species4_thr <- R[round(density_thr*length(R))]

  R <- sort(species5_net[upper.tri(species5_net, diag = FALSE)], decreasing = TRUE)
  species5_thr <- R[round(density_thr*length(R))]

  R <- sort(species6_net[upper.tri(species6_net, diag = FALSE)], decreasing = TRUE)
  species6_thr <- R[round(density_thr*length(R))]
  
} else {
  load(file = comparison_RData)
}

cat(species1_name, "co-expr threshold", format(species1_thr, digits = 3) , "\n")
cat(species2_name, "co-expr threshold", format(species2_thr, digits = 3) , "\n")
cat(species3_name, "co-expr threshold", format(species3_thr, digits = 3) , "\n")
cat(species4_name, "co-expr threshold", format(species4_thr, digits = 3) , "\n")
cat(species5_name, "co-expr threshold", format(species5_thr, digits = 3) , "\n")
cat(species6_name, "co-expr threshold", format(species6_thr, digits = 3) , "\n")

#########################################################################################
# Threshold for neighborhood in module network
density_thr_neigh <- 0.01

R <- sort(species1_net[upper.tri(species1_net, diag = FALSE)], decreasing = TRUE)
species1_thr <- R[round(density_thr*length(R))]
species1_thr_neigh <- R[round(density_thr_neigh*length(R))]

R <- sort(species2_net[upper.tri(species2_net, diag = FALSE)], decreasing = TRUE)
species2_thr <- R[round(density_thr*length(R))]
species2_thr_neigh <- R[round(density_thr_neigh*length(R))]

R <- sort(species3_net[upper.tri(species3_net, diag = FALSE)], decreasing = TRUE)
species3_thr <- R[round(density_thr*length(R))]
species3_thr_neigh <- R[round(density_thr_neigh*length(R))]

R <- sort(species4_net[upper.tri(species4_net, diag = FALSE)], decreasing = TRUE)
species4_thr <- R[round(density_thr*length(R))]
species4_thr_neigh <- R[round(density_thr_neigh*length(R))]

R <- sort(species5_net[upper.tri(species5_net, diag = FALSE)], decreasing = TRUE)
species5_thr <- R[round(density_thr*length(R))]
species5_thr_neigh <- R[round(density_thr_neigh*length(R))]

R <- sort(species6_net[upper.tri(species6_net, diag = FALSE)], decreasing = TRUE)
species6_thr <- R[round(density_thr*length(R))]
species6_thr_neigh <- R[round(density_thr_neigh*length(R))]

#########################################################################################
# Apply the threshold to the adjacency matrix
species1_net_thresholded <- species1_net
species1_net_thresholded[species1_net < species1_thr] <- 0

species2_net_thresholded <- species2_net
species2_net_thresholded[species2_net < species2_thr] <- 0

species3_net_thresholded <- species3_net
species3_net_thresholded[species3_net < species3_thr] <- 0

species4_net_thresholded <- species4_net
species4_net_thresholded[species4_net < species4_thr] <- 0

species5_net_thresholded <- species5_net
species5_net_thresholded[species5_net < species5_thr] <- 0

species6_net_thresholded <- species6_net
species6_net_thresholded[species6_net < species6_thr] <- 0

# Find the edges that satisfy the threshold
# edge_indices1 <- which(species1_net_thresholded > 0, arr.ind = TRUE)
# 
# # Create an edgelist
# edgelist1 <- data.frame(
#   source = edge_indices1[, 1],
#   target = edge_indices1[, 2],
#   weight = species1_net_thresholded[edge_indices1]
# )
# 
# # Remove self-loops if needed (optional)
# edgelist1 <- edgelist1[edgelist1$source != edgelist1$target, ]

# Create a graph object from the adjacency matrix
graph1 <- graph_from_adjacency_matrix(species1_net_thresholded, mode = "undirected", weighted = TRUE, diag = FALSE)
graph2 <- graph_from_adjacency_matrix(species2_net_thresholded, mode = "undirected", weighted = TRUE, diag = FALSE)
graph3 <- graph_from_adjacency_matrix(species3_net_thresholded, mode = "undirected", weighted = TRUE, diag = FALSE)
graph4 <- graph_from_adjacency_matrix(species4_net_thresholded, mode = "undirected", weighted = TRUE, diag = FALSE)
graph5 <- graph_from_adjacency_matrix(species5_net_thresholded, mode = "undirected", weighted = TRUE, diag = FALSE)
graph6 <- graph_from_adjacency_matrix(species6_net_thresholded, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the edgelist with gene names
edgelist1 <- as_data_frame(graph1, what = "edges")
colnames(edgelist1) <- c("Source", "Target", "Weight")

edgelist2 <- as_data_frame(graph2, what = "edges")
colnames(edgelist2) <- c("Source", "Target", "Weight")

edgelist3 <- as_data_frame(graph3, what = "edges")
colnames(edgelist3) <- c("Source", "Target", "Weight")

edgelist4 <- as_data_frame(graph4, what = "edges")
colnames(edgelist4) <- c("Source", "Target", "Weight")

edgelist5 <- as_data_frame(graph5, what = "edges")
colnames(edgelist5) <- c("Source", "Target", "Weight")

edgelist6 <- as_data_frame(graph6, what = "edges")
colnames(edgelist6) <- c("Source", "Target", "Weight")

write.table(edgelist1, file = paste0("data/ComPlEx/networks/MR_network_", 
                                     species1_name,".tsv"), sep = "\t", dec = ".", col.names = TRUE)
write.table(edgelist2, file = paste0("data/ComPlEx/networks/MR_network_", 
                                     species2_name,".tsv"), sep = "\t", dec = ".", col.names = TRUE)
write.table(edgelist3, file = paste0("data/ComPlEx/networks/MR_network_", 
                                     species3_name,".tsv"), sep = "\t", dec = ".", col.names = TRUE)
write.table(edgelist4, file = paste0("data/ComPlEx/networks/MR_network_", 
                                     species4_name,".tsv"), sep = "\t", dec = ".", col.names = TRUE)
write.table(edgelist5, file = paste0("data/ComPlEx/networks/MR_network_", 
                                     species5_name,".tsv"), sep = "\t", dec = ".", col.names = TRUE)
write.table(edgelist6, file = paste0("data/ComPlEx/networks/MR_network_", 
                                     species6_name,".tsv"), sep = "\t", dec = ".", col.names = TRUE)



