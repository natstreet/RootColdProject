#!/usr/bin/env Rscript
# check_deposit_consistency.R — assert the DEG-list invariants of the deposit and fail loudly.
#
# Run from a directory containing the unpacked deposit `data/`. Verifies, for all six
# species x four cold timepoints:
#   1. every gene in each _lfc0 / _lfc0585 / _l2fc1 per-timepoint list has |log2FC| >= 1
#      (all three families carry the post-hoc two-fold cut; the suffix is the results()
#      lfcThreshold, not the fold filter);
#   2. the per-species union of the _lfc0 lists equals the DEGs_<species>.RData object
#      (the DEG set used throughout the paper);
#   3. thresholding padj_<code>.RData at <= 0.05 reproduces the strict _l2fc1 lists
#      (documents that padj_*.RData is the lfcThreshold = 1 test, NOT the _lfc0 DEGs);
#   4. og_summary$n_species_with_DEGs == 6 sums to 270.
# Prints a compact PASS/FAIL table; exits non-zero on any FAIL.

ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
rl <- function(f){ x <- readLines(f); x[x != "" & x != "x"] }
SP <- list(col0=c(code="col0", dir="col0",  degs="DEGs_col0"),
           ost0=c(code="ost0", dir="ost0",  degs="DEGs_ost0"),
           aspen=c(code="pt",  dir="aspen", degs="DEGs_aspen"),
           birch=c(code="bp",  dir="birch", degs="DEGs_birch"),
           spruce=c(code="pa", dir="spruce",degs="DEGs_spruce"),
           pine=c(code="ps",   dir="pine",  degs="DEGs_pine"))
TP <- c("6h","24h","3d","10d"); FAM <- c("lfc0","lfc0585","l2fc1")
fail <- FALSE
say <- function(check, sp, ok, detail=""){
  lab <- if(is.na(ok)) "SKIP" else if(ok) "PASS" else "FAIL"
  cat(sprintf("  [%s] %-26s %-7s %s\n", lab, check, sp, detail))
  if(isFALSE(ok)) fail <<- TRUE
}

for(sp in names(SP)){ m <- SP[[sp]]
  lfc <- ld(sprintf("data/DEGs/log2FoldChange_%s.RData", m["code"]))
  bad <- 0; nfiles <- 0
  for(fam in FAM) for(tp in TP){
    f <- sprintf("data/DEGs/per_timepoint_DEG_lists/%s/DEGs_%s_C_vs_%s_p0.05_%s.txt", m["dir"], m["code"], tp, fam)
    if(!file.exists(f)) next
    nfiles <- nfiles + 1; g <- rl(f); v <- abs(lfc[g, paste0("log2FoldChange_", tp)])
    if(any(v < 1 - 1e-9, na.rm=TRUE)) bad <- bad + 1
  }
  say("1.lists |log2FC|>=1", sp, bad==0, sprintf("(%d files, %d violating)", nfiles, bad))

  u <- unique(unlist(lapply(TP, function(tp) rl(sprintf("data/DEGs/per_timepoint_DEG_lists/%s/DEGs_%s_C_vs_%s_p0.05_lfc0.txt", m["dir"], m["code"], tp)))))
  degf <- sprintf("data/DEGs/%s.RData", m["degs"])
  if(file.exists(degf)){ D <- ld(degf); say("2.union(lfc0)==DEGs", sp, setequal(u, D), sprintf("(union %d, obj %d)", length(u), length(D))) } else say("2.union(lfc0)==DEGs", sp, NA, "(no DEGs object)")

  pf <- sprintf("data/DEGs/padj_%s.RData", m["code"])
  if(file.exists(pf)){ padj <- ld(pf); okall <- TRUE
    for(tp in TP){ col <- paste0("padj_", tp); fromP <- rownames(padj)[!is.na(padj[[col]]) & padj[[col]] <= 0.05]
      l2 <- rl(sprintf("data/DEGs/per_timepoint_DEG_lists/%s/DEGs_%s_C_vs_%s_p0.05_l2fc1.txt", m["dir"], m["code"], tp))
      if(!setequal(fromP, l2)) okall <- FALSE }
    say("3.padj@0.05==l2fc1", sp, okall, "(padj_* = lfcThreshold=1 test)") } else say("3.padj@0.05==l2fc1", sp, NA, "(no padj object)")
}
og <- ld("data/DEGs/og_summary.RData"); n270 <- sum(og$n_species_with_DEGs == 6)
say("4.og_summary all6==270", "-", n270==270, sprintf("(=%d)", n270))
cat(if(fail) "\nRESULT: FAIL\n" else "\nRESULT: PASS\n")
quit(status = if(fail) 1 else 0)
