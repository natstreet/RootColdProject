# ComPlEx pairwise co-expressolog comparison (Netotea et al. 2014)
# author: Elena van Zalen
#
# Builds per-species mutual-rank co-expression networks from the VST expression
# matrices and, for one ordered species pair, tests every orthologue pair with a
# bidirectional hypergeometric test for co-expression-neighbourhood overlap
# (Benjamini-Hochberg FDR). Writes one
#   data/ComPlEx/PairwiseComparisons/comparison-<sp1>-<sp2>-abspearsonMR0.03-table.RData
# Set species1_name / species2_name below and run once for each of the 15 pairs:
#   col0-ost0  col0-aspen  col0-birch  col0-spruce  col0-pine
#              ost0-aspen  ost0-birch  ost0-spruce  ost0-pine
#                          aspen-birch aspen-spruce aspen-pine
#                                      birch-spruce birch-pine
#                                                   spruce-pine

library(arrow)
library(dplyr)
library(gdata)
library(ggplot2)
library(matrixStats)
library(stringr)
library(tibble)
library(tidyr)

# ------------- Fill in below -------------

species1_name <- "col0"
#species1_name <- "ost0"
#species1_name <- "aspen"
#species1_name <- "birch"
#species1_name <- "spruce"

species2_name <- "ost0"
#species2_name <- "aspen"
#species2_name <- "birch"
#species2_name <- "spruce"
#species2_name <- "pine"

run_version <- "2"

# OrthoFinder orthogroups (raw input, not part of the deposit). The derived
# per-pair orthologue table is cached below as data/ComPlEx/orthologs-<sp1>-<sp2>.RData,
# which is shipped in the deposit, so this parquet is only read when regenerating.
ortholog_group_file <- read_parquet("data/orthologs/orthologs.parquet")

# Parameters
cor_method <- "pearson" # pearson spearman
cor_sign <- "abs" # abs
norm_method <- "MR" # CLR MR
density_thr <- 0.03
randomize <- "no" # yes no

# For test-runs
test_run <- "no" # "yes" "no"
# numb_of_cols <-
numb_of_rows_1 <- 3000
numb_of_rows_2 <- 3000

# ------------- Setup --------------
# Output directories (created on demand so a re-run does not fail)
for (d in c("data/ComPlEx/PairwiseComparisons", "data/ComPlEx/centrality",
            "data/ComPlEx/species_net_dens_thr", "data/ComPlEx/correlationPlots/densdist_filter")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Changing the species names of our species of interest for easy of interpretation later on
species_rename_vector <- c("Arabidopsis_thaliana" = "col0", "Populus_tremula" = "aspen", 
                           "Betula_pendula" = "birch", "Picea_abies" = "spruce", "Pinus_sylvestris" = "pine")

ortholog_group_RData <- paste0("data/ComPlEx/orthologs-", species1_name, "-", species2_name, ".RData")
if (!file.exists(ortholog_group_RData)) {
  ortho <- ortholog_group_file %>% 
    mutate(species = dplyr::recode(species, !!!species_rename_vector)) %>% 
    bind_rows(., filter(., grepl("col0", species)) %>%
                mutate(species = gsub("col0", "ost0", species)))
  ortho <- ortho %>%
    filter(species %in% c(species1_name, species2_name))
  ortho_A <- ortho %>% 
    dplyr::select(species,gene,Ortholog_Group) %>% 
    filter(str_detect(species, species1_name)) 
  ortho_B <- ortho %>%
    dplyr::select(species,gene,Ortholog_Group) %>% 
    filter(str_detect(species, species2_name))
  ortho <- inner_join(ortho_A, ortho_B, by = "Ortholog_Group", relationship = "many-to-many") %>% 
    dplyr::select(-c(species.x, species.y)) %>% 
    dplyr::select(Species1 = "gene.x", Species2 = "gene.y", OrthoGroup = "Ortholog_Group")
  save(ortho, file = ortholog_group_RData)
} else {
  load(file = ortholog_group_RData)
}

# Add annotations from arabidopsis
symbols <- read.delim("data/annotation/gene_aliases_20140331.txt", sep = "\t") %>%
  dplyr::rename(Arabidopsis = locus_name, Symbol = symbol, Name = full_name)
ortholog_annot_file <- "data/annotation/Orthogroups_130323_predefined_tree.tsv"

ortholog_annot_RData <- paste0("data/ComPlEx/annotation-", species1_name, "-", species2_name, ".RData")
if (!file.exists(ortholog_annot_RData)) {
  annot <- read.delim2(ortholog_annot_file, header = TRUE, sep = "\t") %>%
    dplyr::rename(Arabidopsis = Aratha.SHORT.pep, OrthoGroup = Orthogroup) %>%
    dplyr::select(Arabidopsis, OrthoGroup) %>%
    filter(Arabidopsis != "") %>%
    separate_rows(Arabidopsis, sep = ", ", convert = FALSE) %>%
    mutate(Arabidopsis = gsub("\\.\\d.p\\d$", "", Arabidopsis)) %>%
    mutate(Arabidopsis = gsub("Aratha_", "", Arabidopsis)) %>%
    left_join(symbols, by = "Arabidopsis") %>%
    group_by(OrthoGroup) %>%
    summarise(Arabidopsis = paste0(unique(Arabidopsis), collapse = "; "), 
              Symbol = paste0(unique(Symbol), collapse = "; "), 
              Name = paste0(unique(Name), collapse = "; ")) %>%
    mutate(Symbol = gsub("NA", "", Symbol),
           Name = gsub("NA", "", Name))
  save(ortho, annot, file = ortholog_annot_RData)
} else {
  load(file = ortholog_annot_RData)
}

species1_transcription_txt <- paste0("data/expression/", species1_name,"_root_cold_stress_expression.txt") 
species2_transcription_txt <- paste0("data/expression/", species2_name,"_root_cold_stress_expression.txt")

if(test_run == "yes") {
  # Only a subset of genes and samples is run
  species1_expr_unfilt <- read.delim(species1_transcription_txt, sep = "\t", header = TRUE)[1:numb_of_rows_1, ]
  species2_expr_unfilt <- read.delim(species2_transcription_txt, sep = "\t", header = TRUE)[1:numb_of_rows_2, ]
  
} else{
  species1_expr_unfilt <- read.delim(species1_transcription_txt, sep = "\t", header = TRUE)
  species2_expr_unfilt <- read.delim(species2_transcription_txt, sep = "\t", header = TRUE)
}

# Filter
# ======
# expression filter rowSds
species1_expr <- species1_expr_unfilt %>%
  {
    .x <- .                              # save original data
    rownames(.x) <- .x[[1]]              # set first column as rownames
    .x <- .x[, -1]                       # drop first column
    mat <- as.matrix(.x)                 # convert to matrix
    row_sd_nonzero <- rowSds(mat) > 0    # SD of each gene across all samples
    row_nonzero_count <- rowCounts(mat != 0)  # how many non-zero values per gene
    
    # Filter: keep genes with SD > 0 and at least 3 non-zero entries
    .x <- .x[row_sd_nonzero & row_nonzero_count >= 3, ]
    
    as_tibble(.x, rownames = "Genes")    # convert back to tibble with gene names
  }
species2_expr <- species2_expr_unfilt %>%
  {
    .x <- .                              # save original data
    rownames(.x) <- .x[[1]]              # set first column as rownames
    .x <- .x[, -1]                       # drop first column
    mat <- as.matrix(.x)                 # convert to matrix
    row_sd_nonzero <- rowSds(mat) > 0    # SD of each gene across all samples
    row_nonzero_count <- rowCounts(mat != 0)  # how many non-zero values per gene
    
    # Filter: keep genes with SD > 0 and at least 3 non-zero entries
    .x <- .x[row_sd_nonzero & row_nonzero_count >= 3, ]
    
    as_tibble(.x, rownames = "Genes")    # convert back to tibble with gene names
  }

cat("Number of genes before SD filtering for ", species1_name, " was ", length(species1_expr_unfilt$Genes), " genes\n", 
    " and after filtering it became ", length(species1_expr$Genes), " genes\n",
    "  for ", species2_name, "it was ", length(species2_expr_unfilt$Genes), " genes before filtering\n", 
    " and after filtering it became ", length(species2_expr$Genes), " genes\n", 
    sep = "")

# Filtering expression data with density distribution threshold
# Compute threshold from density of row means
threshold1 <- {
  tmp <- species1_expr
  expr_mat <- tmp %>%
    column_to_rownames("Genes") %>%
    as.matrix()
  
  d <- density(rowMeans(expr_mat)) 
  x_vals <- d$x
  y_vals <- d$y
  local_mins <- which(diff(sign(diff(y_vals))) == 2)
  
  x_vals[local_mins[1] + 1]# first local minimum
} 
threshold1 <- 1.0
pdf(paste0("data/ComPlEx/correlationPlots/densdist_filter/density_filter_plot_", species1_name, ".pdf"))
plot(d) #, xlim = c(0, 3)
abline(v = threshold1, col = "red", lty = 2)
dev.off()

threshold2 <- {
  tmp <- species2_expr
  expr_mat <- tmp %>%
    column_to_rownames("Genes") %>%
    as.matrix()
  
  d <- density(rowMeans(expr_mat)) 
  x_vals <- d$x
  y_vals <- d$y
  local_mins <- which(diff(sign(diff(y_vals))) == 2)
  
  x_vals[local_mins[1] + 1]
}
threshold2 <- 1.0
pdf(paste0("data/ComPlEx/correlationPlots/densdist_filter/density_filter_plot_", species2_name, ".pdf"))
plot(d)
abline(v = threshold2, col = "red", lty = 2)
dev.off()

# Apply filter and restore "Genes" column
species1_expr <- species1_expr %>%
  filter(rowMeans(across(-Genes)) > threshold1)
species2_expr <- species2_expr %>%
  filter(rowMeans(across(-Genes)) > threshold2)

cat("The density distribution threshold for ", species1_name, " was ", threshold1, " and for ", species2_name, " it was ", threshold2, "\n",
    "Gene count after this filter for ", species1_name, " was ", length(species1_expr$Genes), " genes\n",
    "and the gene count for ", species2_name, " was ",  length(species2_expr$Genes),  " genes\n",
    sep = "")

# Filtering the expression tables where they only contain the genes that have an ortholog in ortho
cat (length(unique(ortho$OrthoGroup)), " ortholog groups containing:\n",
     " ", length(unique(ortho$Species1)), " ", species1_name, " genes\n",
     " ", length(unique(ortho$Species2)), " ", species2_name, " genes\n\n",
     length(unique(species1_expr$Genes)), " expressed ", species1_name, " genes\n",
     length(unique(species2_expr$Genes)), " expressed ", species2_name, " genes\n",
     sep = "")

ortho <- ortho %>%
  filter(Species1 %in% species1_expr$Genes & Species2 %in% species2_expr$Genes)

species1_expr <- species1_expr[species1_expr$Genes %in% ortho$Species1,]
species2_expr <- species2_expr[species2_expr$Genes %in% ortho$Species2,]

cat ("After filtering on expressed genes with ortholog:\n",
     " ", length(unique(ortho$OrthoGroup)), " ortholog groups containing: \n",
     "  ", length(unique(ortho$Species1)), " ", species1_name, " genes\n",
     "  ", length(unique(ortho$Species2)), " ", species2_name, " genes\n",
     sep = "")

comparison_RData <- paste0("data/ComPlEx/PairwiseComparisons/comparison-", species1_name, "-", species2_name, "-", 
                           cor_sign, cor_method, norm_method, density_thr, "-table.RData")

species_net1_RData <- paste0("data/ComPlEx/species_net_dens_thr/species_net_", species1_name, ".RData")
species_net2_RData <- paste0("data/ComPlEx/species_net_dens_thr/species_net_", species2_name, ".RData")

# Capture the start time
start_time <- proc.time()

# ------------- Network calculation/comparison --------------

if (!file.exists(comparison_RData)) {
  if (randomize == "yes") {
    species1_expr$Genes <-
      sample(species1_expr$Genes, nrow(species1_expr), FALSE)
    species2_expr$Genes <-
      sample(species2_expr$Genes, nrow(species2_expr), FALSE)
  }
  # Correlate genes
  species1_net <- cor(t(species1_expr[, -1]), method = cor_method) 
  dimnames(species1_net) <- list(species1_expr$Genes, species1_expr$Genes)               
  
  species2_net <- cor(t(species2_expr[, -1]), method = cor_method)  
  dimnames(species2_net) <- list(species2_expr$Genes, species2_expr$Genes)                 
  
  
  if (cor_sign == "abs") {
    species1_net <- abs(species1_net)
    species2_net <- abs(species2_net)
  }
  #################################################################################
  if (norm_method == "CLR") {                          # currently skipped cause of method selection
    
    z <- scale(species1_net)
    z[z < 0] <- 0
    species1_net <- sqrt(t(z) ** 2 + z ** 2)
    
    z <- scale(species2_net)
    z[z < 0] <- 0
    species2_net <- sqrt(t(z) ** 2 + z ** 2)
    #################################################################################
    
  } else if (norm_method == "MR") {
    C_density_thr_neigh <- 0.01 # for centrality calculation
    
    R <- t(apply(species1_net, 1, rank)) # Apply rank to correlated genes
    species1_net <- sqrt(R * t(R)) # Geometric average 
    species1_thr_neigh <- R[round(C_density_thr_neigh*length(R))] # for centrality calculation
    
    R <- t(apply(species2_net, 1, rank))
    species2_net <- sqrt(R * t(R))
    species2_thr_neigh <- R[round(C_density_thr_neigh*length(R))]
  }
  
  diag(species1_net) <- 0
  diag(species2_net) <- 0
  
  if(!file.exists(species_net1_RData)) {
    save(species1_net, file = species_net1_RData)
  } else {
    print("File already exists")
  }
  
  if(!file.exists(species_net2_RData)) {
    save(species2_net, file = species_net2_RData)
  } else {
    print("File already exists")
  }
  
    ############################################# CENTRALITY ##########################################################
  C_species1_net <- species1_net
  C_species2_net <- species2_net
  
  centrality <- data.frame(Genes = rownames(C_species1_net),
                           Degree = c(NA))
  for (i in 1:nrow(centrality)) {
    neigh <- C_species1_net[centrality$Genes[i],]
    neigh <- names(neigh[neigh >= species1_thr_neigh])
    centrality$Degree[i] <- length(neigh)
  }
  save(centrality, file = paste0("data/ComPlEx/centrality/centrality_",
                                 species1_name, ".RData"))
  
  centrality <- data.frame(Genes = rownames(C_species2_net),
                           Degree = c(NA))
  for (i in 1:nrow(centrality)) {
    neigh <- C_species2_net[centrality$Genes[i],]
    neigh <- names(neigh[neigh >= species2_thr_neigh])
    centrality$Degree[i] <- length(neigh)
  }
  save(centrality, file = paste0("data/ComPlEx/centrality/centrality_",
                                 species2_name, ".RData"))
  
  ############################################# COMPARISON ##########################################################
  R <- sort(species1_net[upper.tri(species1_net, diag = FALSE)], decreasing = TRUE) 
  species1_thr <- R[round(density_thr * length(R))] 
  
  R <-sort(species2_net[upper.tri(species2_net, diag = FALSE)], decreasing = TRUE)
  species2_thr <- R[round(density_thr * length(R))]
  
  comparison <- ortho
  
  comparison$Species1.neigh <- c(NA)
  comparison$Species1.ortho.neigh <- c(NA)
  comparison$Species1.neigh.overlap <- c(NA)
  comparison$Species1.p.val <- c(NA)
  comparison$Species1.effect.size <- c(NA)
  
  comparison$Species2.neigh <- c(NA)
  comparison$Species2.ortho.neigh <- c(NA)
  comparison$Species2.neigh.overlap <- c(NA)
  comparison$Species2.p.val <- c(NA)
  comparison$Species2.effect.size <- c(NA)
  
  for (i in 1:nrow(ortho)) {
    
    if (i %% 100 == 0) {
      cat(i, "\n")
    }
    
    #i <- 18
    # Species 1 -> Species 2
    
    neigh <- species1_net[ortho$Species1[i],]  # Named numeric of all genes and their ranked correlations for gene i (entire co-expression network)
    neigh <- names(neigh[neigh >= species1_thr]) # Retain only the 3% top ranked genes
    
    ortho_neigh <- species2_net[ortho$Species2[i],] # Co-expression network for gene i in species 2 (the ortholog)
    ortho_neigh <- names(ortho_neigh[ortho_neigh >= species2_thr]) # Retain only the 3% top ranked genes
    ortho_neigh <- ortho$Species1[ortho$Species2 %in% ortho_neigh] # Overlapping the the networks, i.e. seeing how many of the species 1 genes have orthologs within ortho.neigh
    
    
    N <- nrow(species1_expr) # Number of all possible genes in S1
    m <- length(neigh) # Number of neighbours of gene i - white balls
    n <- N-m # Number of genes that are NOT neighbours - black balls
    k <- length(unique(ortho_neigh)) # Number of ortholog neighbours  - number of balls we draw
    x <- length(unique(intersect(neigh, ortho_neigh))) # Number of genes that are present in both networks. Must be at least 1.
    p_val <- 1
    effect_size <- 1
    if (x > 1) {
      p_val <- phyper(x-1, m, n, k, lower.tail = FALSE)
      effect_size <- (x/k) / (m/N)
    }
    
    comparison$Species1.neigh[i] <- m
    comparison$Species1.ortho.neigh[i] <- k
    comparison$Species1.neigh.overlap[i] <- x
    comparison$Species1.p.val[i] <- p_val
    comparison$Species1.effect.size[i] <- effect_size
    
    # Species 2 -> Species 1
    
    neigh <- species2_net[ortho$Species2[i],]
    neigh <- names(neigh[neigh >= species2_thr])
    
    ortho_neigh <- species1_net[ortho$Species1[i],]
    ortho_neigh <- names(ortho_neigh[ortho_neigh >= species1_thr])
    ortho_neigh <- ortho$Species2[ortho$Species1 %in% ortho_neigh] 
    N <- nrow(species2_expr)
    m <- length(neigh)
    n <- N-m
    k <- length(unique(ortho_neigh))
    x <- length(unique(intersect(neigh, ortho_neigh)))
    p_val <- 1
    effect_size <- 1
    
    if (x > 1) {
      p_val <- phyper(x-1, m, n, k, lower.tail = FALSE)
      effect_size <- (x/k) / (m/N)
    }
    
    comparison$Species2.neigh[i] <- m
    comparison$Species2.ortho.neigh[i] <- k
    comparison$Species2.neigh.overlap[i] <- x
    comparison$Species2.p.val[i] <- p_val
    comparison$Species2.effect.size[i] <- effect_size
  }
  
  # Remove gene pairs with no overlapping neighbours
  comparison <- comparison %>%
    filter(Species1.neigh.overlap > 0 & Species2.neigh.overlap > 0)
  
  # FDR correction
  comparison$Species1.p.val <- p.adjust(comparison$Species1.p.val, method = "fdr")
  comparison$Species2.p.val <- p.adjust(comparison$Species2.p.val, method = "fdr")
  
  
  comparison_table <- comparison %>%
    rowwise() %>%
    mutate(Max.p.val = max(Species1.p.val, Species2.p.val)) %>%
    mutate(Min.effect.size = min(Species1.effect.size, Species2.effect.size)) %>%
    left_join(annot, by = "OrthoGroup") %>%
    select(-c("Species1.neigh", "Species1.ortho.neigh", "Species2.neigh", "Species2.ortho.neigh")) %>%
    arrange(Max.p.val)
  
  comparison_table$Species1.p.val <- format(comparison_table$Species1.p.val, digits = 3, scientific = TRUE)
  comparison_table$Species2.p.val <- format(comparison_table$Species2.p.val, digits = 3, scientific = TRUE)
  comparison_table$Max.p.val <- format(comparison_table$Max.p.val, digits = 3, scientific = TRUE)
  
  save(comparison_table, file = comparison_RData)
} else{
  print("File already exists")
}

cat ("After filtering on gene pairs in the networks:\n",
     " ", length(unique(comparison$OrthoGroup)), " ortholog groups containing: \n",
     "  ", length(unique(comparison$Species1)), " ", species1_name, " genes\n",
     "  ", length(unique(comparison$Species2)), " ", species2_name, " genes\n",
     sep = "")

# Comparison of p-values of orthologs: species 1 -> species 2 vs species 2 -> species 1
R <- cor.test(-log10(comparison$Species1.p.val), -log10(comparison$Species2.p.val))

pdf(paste0("data/ComPlEx/correlationPlots/ortholog_correlation_plot_",
           species1_name, "-", species2_name,".pdf"))
data.frame(s1 = -log10(comparison$Species1.p.val),
           s2 = -log10(comparison$Species2.p.val)) %>%
  ggplot(aes(x = s1, y = s2)) +
  xlab(paste0(species1_name, " p-value (-log10)")) +
  ylab(paste0(species2_name, " p-value (-log10)")) +
  geom_point() + 
  geom_smooth(method=lm, formula = y ~ x, fill = "gainsboro") +
  ggtitle(paste0("Correlation = ", format(R$estimate, digits = 3))) 
# checking correlation of the p-values between the 2 species
dev.off()

# Capture the end time
end_time <- proc.time()

# Calculate the difference
execution_time <- end_time - start_time
print(execution_time)