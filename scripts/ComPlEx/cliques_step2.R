# **Step 3: identifying cliques**
#   Note that the two previous steps can be skipped if already run (and hopefully saved).

weighted_max_cliques_filterable <- readRDS("data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS")
load("data/ComPlEx/ones_and_zeros_all.RData")

# --------------- DECISIONS ---------------
SELECTION_METHOD <- "SUM" # options: "SUM" = clique sum OR "PVAL" = p-value
TYPE_OF_GENE_SET <- "e" # options: "c" = conserved, "d" = differentiated, "e" = enhanced OR "u" = unique

# --------------- PARAMETERS ---------------
P_VAL <- 0.1 # P_VAL < x
C_SUM <- 15 # C_SUM >= x
CROSSED_PAIRS_PERMITTED <- 2 # CROSSED_PAIRS_PERMITTED <= x
#C_SUM_CLADES <- 3 # C_SUM_CLADES >= x
C_SUM_ANGIO <- 4 # C_SUM_CLADES >= x
C_SUM_GYMNO <- 2 # C_SUM_CLADES >= x

C_SUM_CROSS <- 6 # C_SUM_CROSS <= x
C_SUM_TARGET_CLADE <- 4 # C_SUM_TARGET_CLADES >= x
C_SUM_NON_TARGET_CLADE <- 2 # C_SUM_NON_TARGET_CLADES < x

# NB! When identifying enhanced genes - an additional parameter has been added (04.09.24) removing orthogroups that only have 3-membered cliques. This is to reduce the possibility of overlapping with unique genes.

# **Which parameters?**
#   
#   Conserved
# - p-value: P_VAL
# - clique sum: C_SUM
# 
# Differentiated
# - p-value: P_VAL + CROSSED_PAIRS_PERMITTED
# - clique sum: C_SUM_CLADES + C_SUM_CROSS
# 
# Enhanced
# - p-value: P_VAL
# - clique sum: C_SUM_TARGET_CLADE + C_SUM_NON_TARGET_CLADE + C_SUM_CROSS
# 
# Unique
# - p-value: P_VAL
# - clique sum: C_SUM_TARGET_CLADE

# ----------------CONSERVED GENES------------------------
if(TYPE_OF_GENE_SET == "c"){
  print(TYPE_OF_GENE_SET)
  
  if(SELECTION_METHOD == "PVAL"){
    print(SELECTION_METHOD)
    
    conserved_genes <- weighted_max_cliques_filterable %>%
      filter(MaxpVal < P_VAL) %>%
      group_by(cliqueID) %>%
      mutate(cliqueSize = n()) %>%
      filter(cliqueSize == 15)
    
    print(length(unique(conserved_genes$OrthoGroup)))
    file_name <- paste0("data/ComPlEx/conserved_genes_", SELECTION_METHOD, "_", P_VAL, ".RData")
  }
  
  if(SELECTION_METHOD == "SUM"){
    print(SELECTION_METHOD)
    conserved_genes <- weighted_max_cliques_filterable %>%
      group_by(cliqueID) %>%
      filter(NegLog10CliqueSum >= C_SUM) %>%
      mutate(cliqueSize = n()) %>%
      filter(cliqueSize == 15)
    
    print(length(unique(conserved_genes$OrthoGroup)))
    file_name <- paste0("data/ComPlEx/conserved_genes_", SELECTION_METHOD, "_", C_SUM, ".RData")
  }}
save(conserved_genes, file = file_name)


# ------------------DIFFERENTIATED GENES----------------------

if(TYPE_OF_GENE_SET == "d"){
  print(TYPE_OF_GENE_SET)
  
  if(SELECTION_METHOD == "PVAL"){
    print(SELECTION_METHOD)
    
    conserved_genes <- weighted_max_cliques_filterable %>%
      filter(MaxpVal < P_VAL) %>%
      group_by(cliqueID) %>%
      filter(sum(as.integer(Angio)) == 4 & sum(as.integer(Gymno)) == 2 & sum(as.integer(Cross)) <= 8) %>%
      filter(OriginalCliqueSize == 15)
    
    print(length(unique(conserved_genes$OrthoGroup)))
    file_name <- paste0("data/ComPlEx/differentiated_genes_", SELECTION_METHOD, "_", P_VAL,"_XPAIRS_",CROSSED_PAIRS_PERMITTED ,".RData")
  }
  
  if(SELECTION_METHOD == "SUM"){
    print(SELECTION_METHOD)
    
    conserved_genes <- weighted_max_cliques_filterable %>%
      group_by(cliqueID) %>%
      #filter(sum(Angio) == 4 & sum(Gymno) == 2 & sum(Cross) <= 6) %>%
      filter(AngioSum >= C_SUM_ANGIO & GymnoSum >= C_SUM_GYMNO & CrossSum < 2) %>%
      mutate(cliqueSize = n()) %>%
      filter(cliqueSize == 15) # look at expression profiles of these genes
    
    print(length(unique(conserved_genes$OrthoGroup)))
    file_name <- paste0("data/ComPlEx/differentiated_genes_", SELECTION_METHOD, "_C-X_", C_SUM_ANGIO, "-", C_SUM_GYMNO,"-", C_SUM_CROSS, ".RData")
  }}

# -----------------ENHANCED GENES----------------------- (clade diverged)
# clade specific cliques within the orthogroups. Look at expression profiles
# make arabidopsis and broadleaf set
if(TYPE_OF_GENE_SET == "e"){
  print(TYPE_OF_GENE_SET)
  
  if(SELECTION_METHOD == "PVAL"){
    print(SELECTION_METHOD)
    
    conserved_genes_angio <- weighted_max_cliques_filterable %>%
      filter(MaxpVal < P_VAL) %>%
      filter(!(OriginalCliqueSize == 3)) %>%  # Only orthogroups with cliques larger than 3.
      group_by(cliqueID) %>%
      filter(sum(Angio) == 4) %>%
      filter(sum(Gymno) == 0 & sum(Cross) == 0)
    
    print(length(unique(conserved_genes_angio$OrthoGroup)))
    
    conserved_genes_gymno <- weighted_max_cliques_filterable %>%
      filter(MaxpVal < P_VAL) %>%
      filter(!(OriginalCliqueSize == 3)) %>% # Only orthogroups with cliques larger than 3.
      group_by(cliqueID) %>%
      filter(sum(Gymno) == 1) %>%
      filter(sum(Angio) == 0 & sum(Cross) == 0)
    
    print(length(unique(conserved_genes_gymno$OrthoGroup)))
    
    file_name_angio <- paste0("data/ComPlEx/enhanced_angiosperm_genes_", SELECTION_METHOD, "_", P_VAL, ".RData")
    file_name_gymno <- paste0("data/ComPlEx/enhanced_gymnosperm_genes_", SELECTION_METHOD, "_", P_VAL, ".RData")
  }
  
  if(SELECTION_METHOD == "SUM"){
    print(SELECTION_METHOD)
    
    conserved_genes_angio <- weighted_max_cliques_filterable %>% 
      filter(GymnoSum < C_SUM_NON_TARGET_CLADE & CrossSum < C_SUM_CROSS)
    
    conserved_genes_angio <- weighted_max_cliques_filterable %>% 
      filter(AngioSum >= C_SUM_TARGET_CLADE) %>% 
      filter(GymnoSum < C_SUM_NON_TARGET_CLADE & CrossSum == 0) %>% 
      group_by(cliqueID) %>% 
      filter(sum(Angio) == 3) %>% 
      mutate(cliqueSize = n()) %>% 
      filter(cliqueSize > 3) # Removes orthogroups with only 3-membered cliques.
    
    print(length(unique(conserved_genes_angio$OrthoGroup)))
    
    conserved_genes_gymno <- weighted_max_cliques_filterable %>% 
      filter(GymnoSum >= C_SUM_NON_TARGET_CLADE) %>% 
      filter(AngioSum < C_SUM_NON_TARGET_CLADE & CrossSum == 0) %>% 
      group_by(cliqueID) %>% 
      filter(sum(Gymno) == 1) %>% 
      mutate(cliqueSize = n()) %>% 
      filter(cliqueSize > 3) # Removes orthogroups with only 3-membered cliques.
    
    print(length(unique(conserved_genes_gymno$OrthoGroup)))
    
    file_name_angio <- paste0("data/ComPlEx/enhanced_angiosperm_genes_", SELECTION_METHOD, "_T-NT-X_", C_SUM_TARGET_CLADE,"-",C_SUM_NON_TARGET_CLADE, "-", C_SUM_CROSS,".RData")
    file_name_gymno <- paste0("data/ComPlEx/enhanced_gymnosperm_genes_", SELECTION_METHOD, "_T-NT-X_", C_SUM_TARGET_CLADE,"-",C_SUM_NON_TARGET_CLADE, "-", C_SUM_CROSS, ".RData")
  }}

# ------------------UNIQUE GENES----------------------
if(TYPE_OF_GENE_SET == "u"){
  print(TYPE_OF_GENE_SET)
  
  if(SELECTION_METHOD == "PVAL"){
    print(SELECTION_METHOD)
    
    angio_unique_ogs <- ones_and_zeros_all %>% 
      mutate(Angiosperms = rowSums(. [1:6])) %>% 
      mutate(Cross = rowSums(. [7:14]))%>%  
      mutate(Gymnosperms = rowSums(. [15])) %>% 
      mutate(Conserved = rowSums(. [1:15])) %>% 
      filter(Angiosperms == 4 & Gymnosperms == 0 & Cross == 0)
    
    angio_unique_ogs_vctr <- rownames(angio_unique_ogs)
    
    conserved_genes_angio <- weighted_max_cliques_filterable %>%
      filter(OrthologGroup %in% angio_unique_ogs_vctr) %>%
      filter(MaxpVal < P_VAL) %>%
      group_by(cliqueID) %>%
      filter(sum(Angio) == 4)
    
    print(length(unique(conserved_genes_angio$OrthologGroup)))
    
    gymno_unique_ogs <- ones_and_zeros_all %>% 
      mutate(Angiosperms = rowSums(. [1:6])) %>% 
      mutate(Cross = rowSums(. [7:14]))%>%  
      mutate(Gymnosperms = rowSums(. [15])) %>% 
      mutate(Conserved = rowSums(. [1:15])) %>% 
      filter(Angiosperms == 0 & Gymnosperms == 2 & Cross == 0)
    
    gymno_unique_ogs_vctr <- rownames(gymno_unique_ogs)
    
    conserved_genes_gymno <- weighted_max_cliques_filterable %>%
      filter(OrthologGroup %in% gymno_unique_ogs_vctr) %>%
      filter(MaxpVal < P_VAL) %>%
      group_by(cliqueID) %>%
      filter(sum(Gymno) == 2)
    print(length(unique(conserved_genes_gymno$OrthoGroup)))
    
    file_name_angio <- paste0("data/ComPlEx/unique_angiosperm_genes_", SELECTION_METHOD, "_", P_VAL, ".RData")
    file_name_gymno <- paste0("data/ComPlEx/unique_gymnosperm_genes_", SELECTION_METHOD, "_", P_VAL, ".RData")
  }
  
  if(SELECTION_METHOD == "SUM"){
    print(SELECTION_METHOD)
    
    angio_unique_ogs <- ones_and_zeros_all %>% 
      mutate(Angiosperms = rowSums(. [1:6])) %>% 
      mutate(Cross = rowSums(. [7:14]))%>%  
      mutate(Gymnosperms = rowSums(. [15])) %>% 
      mutate(Conserved = rowSums(. [1:15])) %>% 
      filter(Angiosperms == 4 & Gymnosperms == 0 & Cross == 0)
    
    angio_unique_ogs_vctr <- rownames(angio_unique_ogs)
    
    conserved_genes_angio <- weighted_max_cliques_filterable %>%
      filter(OrthoGroup %in% angio_unique_ogs_vctr) %>%
      filter(AngioSum >= C_SUM_TARGET_CLADE)
    
    print(length(unique(conserved_genes_angio$OrthoGroup)))
    
    gymno_unique_ogs <- ones_and_zeros_all %>% 
      mutate(Angiosperms = rowSums(. [1:6])) %>% 
      mutate(Cross = rowSums(. [7:14]))%>%  
      mutate(Gymnosperms = rowSums(. [15])) %>% 
      mutate(Conserved = rowSums(. [1:15])) %>% 
      filter(Angiosperms == 0 & Gymnosperms == 2 & Cross == 0)
    
    gymno_unique_ogs_vctr <- rownames(gymno_unique_ogs)
    
    conserved_genes_gymno <- weighted_max_cliques_filterable %>%
      filter(OrthoGroup %in% gymno_unique_ogs_vctr) %>%
      filter(GymnoSum >= C_SUM_TARGET_CLADE)
    
    print(length(unique(conserved_genes_gymno$OrthoGroup)))
    
    file_name_angio <- paste0("data/ComPlEx/unique_angiosperm_genes_", SELECTION_METHOD,"_T_", C_SUM_TARGET_CLADE ,".RData")
    file_name_gymno <- paste0("data/ComPlEx/unique_gymnosperm_genes_", SELECTION_METHOD, "_T_", C_SUM_TARGET_CLADE ,".RData")
  }}

#SAVE
# add arabidopsis and broadleaf 
# broadleaf and gymnosperm comparison
save(conserved_genes, file = "data/ComPlEx/conserved_genes_step3.RData")
save(conserved_genes_angio, file = "data/ComPlEx/conserved_genes_angio_step3.RData")
save(conserved_genes_gymno, file = "data/ComPlEx/conserved_genes_gymno_step3.RData")

