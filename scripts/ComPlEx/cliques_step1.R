# Clique extraction, step 1 (author: Elena van Zalen)
#
# From the co-expressologs across the 15 pairwise comparisons, keep orthogroups
# that contain at least one clique and enumerate the maximal cliques per orthogroup,
# producing the filterable clique table used downstream. The deposit ships the final
# object as data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS; conserved-core
# selection (step 3) is in cliques_step2.R.

library(plyr)
library(dplyr)
library(tidyverse)
library(igraph)

dir.create("data/ComPlEx/cliques", recursive = TRUE, showWarnings = FALSE)

# **JUMP TO STEP 3 FOR IDENTIFICATION OF VARIOUS CLIQUES**
#
#   **Step 1: remove orthogroups with no cliques**
#   PRODUCES: CONS_OGS_RED which contains all orthogroups with at least one clique (3-6 members)

load("data/ComPlEx/all_expressed_genes.RData")
expressologs <- expr_genes %>%
  dplyr::rename(OrthoGroup = OrthologGroup) %>% 
  dplyr::rename(MaxpVal = MaxpVal) %>%
  dplyr::filter(MaxpVal < 0.05) #0.05 = 155186 rows, 0.9 = 664970 rows
#save(expressologs09,file="data/ComPlEx/expressologs09.RData") 

load("data/ComPlEx/ones_and_zeros_all.RData")

# Add summary columns and remove all orthogroups that no expressed gene pairs.
df_with_sums <- ones_and_zeros_all %>% 
  dplyr::mutate(Angiosperms = rowSums(. [1:6])) %>% 
  dplyr::mutate(Cross = rowSums(. [7:14]))%>%   
  dplyr::mutate(Gymnosperms = rowSums(. [15])) %>% 
  dplyr::mutate(Conserved = rowSums(. [1:15])) %>% 
  dplyr::filter(Conserved > 0)

cons_ogs <- df_with_sums %>%
  dplyr::filter(rownames(df_with_sums) %in% expressologs$OrthoGroup) %>%
  rownames()

cons_ogs <- cons_ogs[!(cons_ogs %in% c("OG0000000", "OG0000001", "OG0000003", "OG0000004", "OG0000006", "OG0000008",
                                        "OG0000013", "OG0000014"))] # I removed 5-9-10-15-107, see below
#OG0000000 149612, OG0000001 61531, OG0000003 46948,
# OG0000004 (28739), OG0000005 (14394), OG0000006 (23161), OG0000008 (31467), OG0000009 (11508), OG0000010 (16343)
# OG0000013 (19173), OG0000014 (25463), OG0000015 (10232), OG0000107 (10810)

#Loop over all orthogroups from above and count the number of max cliques that contain at least three members.

numb_of_cliques_per_og <- c()

for (g in cons_ogs) {
 print(g)
   # g <- "OG0001993"

 expressologs_g <- expressologs %>%
   dplyr::filter(OrthoGroup == g)

   expressologs_g <- expressologs_g %>%
     dplyr::mutate(UGeneSpecies1 = paste0(expressologs_g$Species1, "-", expressologs_g$GeneSpecies1),
                   UGeneSpecies2 = paste0(expressologs_g$Species2, "-", expressologs_g$GeneSpecies2))

   nodes <- data.frame(name = unique(c(expressologs_g$UGeneSpecies1, expressologs_g$UGeneSpecies2)))

   edges <- data.frame(from = expressologs_g$UGeneSpecies1,
                       to = expressologs_g$UGeneSpecies2)
   net <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

# largest_clique <- largest_weighted_cliques(net)
largest_clique <- max_cliques(net, min = 6) # 699 OG with at least one clique, 87 G with co-expressologs <0.05 filter
numb_of_cliques_per_og <- rbind(numb_of_cliques_per_og, data.frame(Orthogroup = g,
                                                                    numbCliques = length(largest_clique),
                                                                    minSize = min(lengths(largest_clique)),
                                                                    maxSize = max(lengths(largest_clique))))
}

# 2) Remove orthogroups with no cliques (min clique size = 3). The loop simply goes through all orthogroups and counts the number of cliques that contain at least three members. We therefore need to remove all orthogroups with zero cliques.
numb_of_cliques_per_og <- numb_of_cliques_per_og %>%
  dplyr::mutate(numbCliques = as.numeric(numbCliques)) %>%
  dplyr::filter(is.finite(numbCliques) & numbCliques > 0) %>% # added this line to account for Inf values
  dplyr::arrange(desc(numbCliques))

paste0("Number of orthologs with at least one clique (3 to 6 members): ", length(unique(numb_of_cliques_per_og$Orthogroup)))

# 3) Vectorise and save!
cons_og_red <- numb_of_cliques_per_og$Orthogroup # added "min6" for saving
save(cons_og_red, file = "data/ComPlEx/cliques/conserved_orthogroups_reduced_min6_coex005.RData")
save(numb_of_cliques_per_og, file = "data/ComPlEx/cliques/number_of_cliques_per_OG_min6_coex005.RData")

# **Step 2: The clique code**
#   PRODUCES: clique_genes_prefiltered (if saved) and clique_genes_filterable. The latter is the final file which can be used for identifying all clique types.
# 
#load("data/ComPlEx/number_of_cliques_per_OG_min6_coex005.RData") # consider dropping OG0000001
#load("data/ComPlEx/conserved_orthogroups_reduced_min6_coex005.RData") # orthogroups that contain at least one clique (size 3-6)
expressologs_from_max_cliques <- list()
max_numb_of_cliques <- 5500

for (g in cons_og_red) {
 print(g)
 #g <- "OG0000025"
 expressologs_g <- expressologs %>%
   dplyr::filter(OrthoGroup == g)
 expressologs_g <- expressologs_g %>%
   dplyr::mutate(UGeneSpecies1 = paste0(expressologs_g$Species1, "-", expressologs_g$GeneSpecies1),
                 UGeneSpecies2 = paste0(expressologs_g$Species2, "-", expressologs_g$GeneSpecies2))
 nodes <- data.frame(name = unique(c(expressologs_g$UGeneSpecies1, expressologs_g$UGeneSpecies2)))
 edges <- data.frame(from = expressologs_g$UGeneSpecies1,
                     to = expressologs_g$UGeneSpecies2)
 net <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
 largest_clique <- max_cliques(net, min = 6)
 if(length(largest_clique) > max_numb_of_cliques){
   largest_clique_red <- largest_clique[1:max_numb_of_cliques]
 } else{largest_clique_red <- largest_clique}
 for(c in 1:length(largest_clique_red)){
   #c <- 1 this only takes the clique annotated with "_1"
   if (c %% 500 == 0) {
     cat(c, "\n")
   }
   clique <- largest_clique[[c]]
   clique_name <- attr(clique, "names")
   clique_species <- str_split_fixed(clique_name, "-", n = 2)[,1]
   clique_genes <-  str_split_fixed(clique_name, "-", n = 2)[,2]

   expressologs_in_clique <- expressologs_g %>%
     dplyr::filter(UGeneSpecies1 %in% clique_name & UGeneSpecies2 %in% clique_name) %>%
     dplyr::select(c(1:5, 8)) %>%
     dplyr::mutate(cliqueID = paste0(g,"_" ,c))
   expressologs_from_max_cliques <- append(expressologs_from_max_cliques, list(expressologs_in_clique))
 }}
clique_genes_unlisted <- plyr::ldply(expressologs_from_max_cliques)
clique_genes_unlisted_prefilter <- clique_genes_unlisted %>%
 dplyr::group_by(cliqueID) %>%
 dplyr::mutate(MaxpValNegLog10 = -log10(MaxpVal)) %>%
 dplyr::mutate(NegLog10CliqueSum = sum(MaxpValNegLog10)) %>%
 dplyr::ungroup()

# # ------------- SAVING POINT -----------------
# # As the code above does take some time to run - it may be a good idea to save progress here.
save(clique_genes_unlisted_prefilter, file = "data/ComPlEx/cliques/clique_genes_prefiltered_5500_coex005.RData")
# # load("Data/DATA/clique_genes_prefiltered.RData")
#load("data/ComPlEx/clique_genes_prefiltered.RData")
weighted_max_cliques <- clique_genes_unlisted_prefilter %>%
  dplyr::mutate(species_pair = paste0(Species1, Species2))

speciesPairs <- tibble(SpeciesPair = weighted_max_cliques %>% distinct(species_pair) %>% pull(species_pair),
                       SpeciesPairClade = c("Angio", "Angio", "Angio", "Cross", "Cross", "Angio", "Angio", "Cross",
                                            "Cross", "Angio", "Cross", "Cross", "Cross", "Cross", "Gymno"))
#("Cross", "Cross", "Gymno", "Angio", "Cross", "Angio", "Cross", "Cross",
#  "Cross", "Angio", "Angio", "Angio", "Cross", "Angio", "Cross"))
weighted_max_cliques <-  left_join(weighted_max_cliques, speciesPairs, by = join_by("species_pair" == "SpeciesPair" ))

weighted_max_cliques_filterable <- weighted_max_cliques %>%
  dplyr::mutate(Clade = paste0(SpeciesPairClade, "-sum" )) %>%
  dplyr::group_by(cliqueID) %>%
  dplyr::mutate(OriginalCliqueSize = n()) %>%
  dplyr::mutate(OriginalCliqueSize = factor(OriginalCliqueSize)) %>%
  dplyr::ungroup() %>%
  pivot_wider(names_from = SpeciesPairClade, values_from = SpeciesPairClade)%>%
  pivot_wider(names_from = Clade, values_from = MaxpValNegLog10)

weighted_max_cliques_filterable <- replace_na(weighted_max_cliques_filterable, list(`Angio-sum` = 0, `Gymno-sum` = 0, `Cross-sum` = 0))

weighted_max_cliques_filterable <- weighted_max_cliques_filterable %>%
  dplyr::group_by(cliqueID) %>%
  dplyr::mutate(AngioSum = sum(`Angio-sum`)) %>%
  dplyr::mutate(GymnoSum = sum(`Gymno-sum`)) %>%
  dplyr::mutate(CrossSum = sum(`Cross-sum`))

weighted_max_cliques_filterable$Angio <- lapply(weighted_max_cliques_filterable$Angio, function(x) {
  if (is.null(x)) {
    return(0)  # Replace NULL with 0
  } else if (is.character(x)) {
    return(1)  # Replace character vectors with 1
  }
})
weighted_max_cliques_filterable$Gymno <- lapply(weighted_max_cliques_filterable$Gymno, function(x) {
  if (is.null(x)) {
    return(0)  # Replace NULL with 0
  } else if (is.character(x)) {
    return(1)  # Replace character vectors with 1
  }
})
weighted_max_cliques_filterable$Cross <- lapply(weighted_max_cliques_filterable$Cross, function(x) {
  if (is.null(x)) {
    return(0)  # Replace NULL with 0
  } else if (is.character(x)) {
    return(1)  # Replace character vectors with 1
  }
})

weighted_max_cliques_filterable <- weighted_max_cliques_filterable %>%
  dplyr::ungroup() %>%
  dplyr::mutate_at(c("Angio", "Gymno", "Cross"), as.numeric) %>%
  dplyr::mutate(UGeneSpecies1 = paste0(weighted_max_cliques_filterable$Species1, "-", weighted_max_cliques_filterable$GeneSpecies1),
                UGeneSpecies2 = paste0(weighted_max_cliques_filterable$Species2, "-", weighted_max_cliques_filterable$GeneSpecies2), .before = species_pair)

# ------------- SAVING POINT -----------------
# At this point the data set is ready to be used for identification of cliques.
save(weighted_max_cliques_filterable, file = "data/ComPlEx/cliques/clique_genes_filterable_coex005.RData")

#  