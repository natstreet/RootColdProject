#' ---
#' title: "Root cold stress network setup"
#' subtitle: "Per-species variance-stabilising transformation of the DESeq2 objects
#'   into the co-expression input matrices (data/expression/<species>_root_cold_stress_expression.txt).
#'   Shown here for Scots pine (dds_Ps); the other species follow the same pattern
#'   (see the commented blocks below)."
#' author: "Elena van Zalen"
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
suppressPackageStartupMessages({
  library(arrow)
  library(cowplot)
  library(DESeq2)
  library(DT)
  library(gplots)
  library(here)
  library(matrixStats)
  library(RColorBrewer)
  library(tidyverse)
})
theme_set(theme_classic())
theme_update(plot.title = element_text(face="bold"))

#' Loading the data
#load(here("data/dds_col0_soil.rda"))
#load(here("data/dds_ost0_soil.rda"))
#load(here("data/dds_Pt_new.rda"))
#load(here("data/dds_Bp_new.rda"))
#load("data/dds/dds_Pa.rda")
load("data/dds/dds_Ps.rda")

# Generate expression files
#' # Variance Stabilising Transformation
#vsd_col0_s <- varianceStabilizingTransformation(dds_col0_soil, blind=TRUE)
#vst_col0_s <- assay(vsd_col0_s)
#vst_col0_s <- vst_col0_s - min(vst_col0_s)
#col0_s <- data.frame(Genes = row.names(vst_col0_s), vst_col0_s)
#write.table(col0_s, file = here("data/col0_soil_expression.txt"), sep = "\t", dec = ".", col.names = TRUE)

#vsd_ost0_s <- varianceStabilizingTransformation(dds_Ost0_soil, blind=TRUE)
#vst_ost0_s <- assay(vsd_ost0_s)
#vst_ost0_s <- vst_ost0_s - min(vst_ost0_s)
#ost0_s <- data.frame(Genes = row.names(vst_ost0_s), vst_ost0_s)
#write.table(ost0_s, file = here("data/ost0_soil_expression.txt"), sep = "\t", dec = ".", col.names = TRUE)

#vsd_Pt <- varianceStabilizingTransformation(dds_Pt, blind=TRUE)
#vst_Pt <- assay(vsd_Pt)
#vst_Pt <- vst_Pt - min(vst_Pt)
#Pt <- data.frame(Genes = row.names(vst_Pt), vst_Pt)
#write.table(Pt, file = here("data/Pt_expression.txt"), sep = "\t", dec = ".", col.names = TRUE)

#vsd_Bp <- varianceStabilizingTransformation(dds_Bp, blind=TRUE)
#vst_Bp <- assay(vsd_Bp)
#vst_Bp <- vst_Bp - min(vst_Bp)
#Bp <- data.frame(Genes = row.names(vst_Bp), vst_Bp)
#write.table(Bp, file = here("data/Bp_expression.txt"), sep = "\t", dec = ".", col.names = TRUE)

# vsd_Pa <- varianceStabilizingTransformation(dds_Pa, blind=TRUE)
# vst_Pa <- assay(vsd_Pa)
# vst_Pa <- vst_Pa - min(vst_Pa)
# save(vst_Pa, file ="data/vst_Pa.rda")
# Pa <- data.frame(Genes = row.names(vst_Pa), vst_Pa)
# write.table(Pa, file = "data/Pa_expression.txt", sep = "\t", dec = ".", col.names = TRUE)
# # 
vsd_Ps <- varianceStabilizingTransformation(dds_Ps, blind=TRUE)
vst_Ps <- assay(vsd_Ps)
vst_Ps <- vst_Ps - min(vst_Ps)
save(vst_Ps, file ="data/vst_Ps.rda")

Ps <- data.frame(Genes = row.names(vst_Ps), vst_Ps)
write.table(Ps, file = "data/expression/pine_root_cold_stress_expression.txt", sep = "\t", dec = ".", col.names = TRUE)
# 









