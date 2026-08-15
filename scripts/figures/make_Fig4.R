#!/usr/bin/env Rscript
# make_Fig4.R — Figure 4: conserved co-expression across the three species groups.
#   (A) Triple Venn of the co-expressolog orthogroups shared within each lineage pair —
#       Arabidopsis (col0|ost0), broadleaf trees (aspen|birch) and conifers (pine|spruce),
#       at MaxpVal < 0.1. Areas are computed from the data, not hard-coded.
#   (B) GO (biological-process) enrichment of three Venn regions (top terms): the shared core
#       ABC and the Arabidopsis–broadleaf AB set (Arabidopsis genes), and the broadleaf–conifer
#       BC set (aspen genes).
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs: data/ComPlEx/co_expressologs.RData, data/ComPlEx/all_go_results_clades.RData
# Output: Figure4.pdf, Figure4.png
# Usage:  Rscript scripts/figures/make_Fig4.R

suppressPackageStartupMessages({ library(VennDiagram); library(ggplot2); library(dplyr)
  library(cowplot); library(ggplotify); library(grid); library(futile.logger) })
flog.threshold(ERROR)                               # silence VennDiagram's log file chatter
ld <- function(f){ e<-new.env(); load(f,envir=e); get(ls(e)[1],envir=e) }

# ── panel A: Venn areas from the data ────────────────────────────────────────
co <- ld("data/ComPlEx/co_expressologs.RData")
co$pairkey <- apply(cbind(co$Species1, co$Species2), 1, function(v) paste(sort(v), collapse="|"))
ogset <- function(k, p=0.1) unique(co$OrthologGroup[co$pairkey==k & co$MaxpVal<p])
A <- ogset("col0|ost0"); B <- ogset("aspen|birch"); C <- ogset("pine|spruce")
vg <- draw.triple.venn(area1=length(A), area2=length(B), area3=length(C),
  n12=length(intersect(A,B)), n13=length(intersect(A,C)), n23=length(intersect(B,C)),
  n123=length(Reduce(intersect,list(A,B,C))),
  category=c("Arabidopsis","Broadleaf","Conifer"), fill=c("#7FBF9E","#F4C08A","#B9B9CE"),
  alpha=rep(0.55,3), lwd=1.4, col=rep("grey30",3), cex=rep(1.1,7), fontface="bold", fontfamily="sans",
  cat.col=c("#2E8B57","#D2691E","#6A5ACD"), cat.cex=rep(1.2,3), cat.fontface="bold",
  cat.pos=c(-20,20,180), cat.dist=c(0.06,0.06,0.045), ind=FALSE)
grid.newpage(); PA <- as.ggplot(grid::gTree(children=vg))

# ── panel B: GO dot-plots from the deposited clade GO results ─────────────────
go <- do.call(c, ld("data/ComPlEx/all_go_results_clades.RData"))   # flatten per-cell lists -> named list
pan <- c("ABC (Arabidopsis)"="Col0 ABC_DE", "AB (Arabidopsis)"="Col0 AB_DE", "BC (aspen)"="Aspen BC_DE")
subplot <- function(key, title, n=10, col){
  d <- go[[key]]
  d <- d[order(-d$scores), , drop=FALSE][seq_len(min(n, nrow(d))), ]
  d$Term <- factor(d$Term, levels=rev(d$Term))
  ggplot(d, aes(scores, Term)) +
    geom_segment(aes(x=0, xend=scores, yend=Term), colour="grey80", linewidth=0.4) +
    geom_point(aes(size=as.integer(Significant)), colour=col) +
    scale_size_continuous(range=c(2,6.5), name="Genes") +
    scale_x_continuous(expand=expansion(mult=c(0.01,0.08))) +
    labs(x=expression(-log[10]~italic(P)), y=NULL, title=title) +
    theme_bw(base_size=10) +
    theme(plot.title=element_text(face="bold", size=10), panel.grid.major.y=element_line(colour="grey93"))
}
cols <- c("#2E8B57","#B8860B","#6A5ACD")
Bs <- Map(function(k,t,c) subplot(k,t,col=c), pan, names(pan), cols)
PB <- plot_grid(plotlist=Bs, ncol=1, align="v")

fig <- plot_grid(PA, PB, ncol=1, labels=c("A","B"), label_size=16, rel_heights=c(0.85,1.7))
ggsave("Figure4.pdf", fig, width=8, height=12)
ggsave("Figure4.png", fig, width=8, height=12, dpi=110)
cat(sprintf("Venn areas: A=%d B=%d C=%d  ABC=%d\n", length(A), length(B), length(C),
            length(Reduce(intersect,list(A,B,C)))))
