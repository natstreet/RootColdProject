#!/usr/bin/env Rscript
# make_Fig3.R — Figure 3: characterisation of the four super clusters (SC_I–SC_IV).
#   (A) Clustering patterns: mean z-scored expression profile of each super cluster's member
#       genes over the cold time course (mean line, +/-1 SD ribbon).
#   (B) Cluster composition: number of member genes per species, split into transcription
#       factors (TFs) and the remaining DEGs.
#   (C) Conservation: distribution of each super cluster's orthogroups across the three
#       phylogenetic groups (Arabidopsis / Deciduous / Coniferous), for DEGs and TF-DEGs.
#
# Run from a directory containing the unpacked figshare deposit `data/`. Panels A and B use the
# per-species super cluster membership in data/superclusters/membership/; panel C uses the pooled
# super cluster gene lists and the orthogroup table (as in make_Fig3C.R).
#
# Output: Figure3.pdf, Figure3.png
# Usage:  Rscript scripts/figures/make_Fig3.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

SC      <- c(SC_I = 6, SC_II = 2, SC_III = 3, SC_IV = 1)
species <- c(col0="Col0", ost0="Ost0", birch="Birch", aspen="Aspen", pine="Pine", spruce="Spruce")
pre     <- c(col0="Col", ost0="Ost", aspen="Pt", birch="Bp", spruce="Pa", pine="Ps")
condlab <- c("01"="control","02"="6h","03"="24h","04"="3d","05"="10d")
ld <- function(f){ e<-new.env(); load(f,envir=e); get(ls(e)[1],envir=e) }
memb <- function(sp,n) read.csv(sprintf("data/superclusters/membership/%s/cluster_%d_members.csv",sp,n),
                                stringsAsFactors=FALSE)$gene_id
TFset <- unique(sub("\\.\\d+$","",unlist(lapply(Sys.glob("data/TranscriptionFactors/TF_*.RData"),
          function(f){x<-ld(f); if(is.data.frame(x)) x[[if("Gene_ID"%in%colnames(x))"Gene_ID" else 1]] else x}))))
EXP <- lapply(names(species), function(sp)
  read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt",sp),check.names=FALSE,row.names=1))
names(EXP) <- names(species)

# ── panel A: mean z-scored profile per super cluster ─────────────────────────
Arows <- list()
for (lab in names(SC)) { z <- list()
  for (sp in names(species)) {
    g <- intersect(memb(sp, SC[[lab]]), rownames(EXP[[sp]])); if (!length(g)) next
    ex <- EXP[[sp]][g, , drop=FALSE]; cc <- sub(paste0("^",pre[sp],"\\d+\\."),"",colnames(ex))
    cm <- t(apply(ex, 1, function(v) tapply(as.numeric(v), factor(cc,levels=names(condlab)), mean)))
    z[[sp]] <- t(apply(cm, 1, function(r) (r-mean(r))/sd(r)))
  }
  Z <- do.call(rbind, z)
  Arows[[lab]] <- data.frame(SC=lab, cond=names(condlab), m=colMeans(Z,na.rm=TRUE), s=apply(Z,2,sd,na.rm=TRUE))
}
A <- bind_rows(Arows); A$SC<-factor(A$SC,levels=names(SC)); A$cond<-factor(A$cond,levels=names(condlab))
pA <- ggplot(A, aes(cond, m, group=1)) +
  geom_ribbon(aes(ymin=m-s, ymax=m+s), fill="grey70", alpha=0.4) +
  geom_line() + geom_point(size=1.3) +
  facet_wrap(~SC, ncol=1, strip.position="right") +
  scale_x_discrete(labels=condlab) +
  labs(x=NULL, y="Mean z-score", title="(A) Clustering patterns") +
  theme_bw(base_size=10) + theme(panel.grid.minor=element_blank(),
    axis.text.x=element_text(angle=45,hjust=1))

# ── panel B: composition (TFs vs remaining DEGs) per species ──────────────────
Brows <- list()
for (lab in names(SC)) for (sp in names(species)) {
  g <- memb(sp, SC[[lab]]); nt <- sum(g %in% TFset)
  Brows[[paste(lab,sp)]] <- data.frame(SC=lab, Species=species[[sp]],
    part=c("DEGs","TFs"), n=c(length(g)-nt, nt))
}
B <- bind_rows(Brows); B$SC<-factor(B$SC,levels=names(SC))
B$Species<-factor(B$Species,levels=unname(species)); B$part<-factor(B$part,levels=c("DEGs","TFs"))
pB <- ggplot(B, aes(Species, n, fill=part)) + geom_col() +
  facet_wrap(~SC, ncol=1, scales="free_y", strip.position="right") +
  scale_fill_manual(values=c(DEGs="#5B6C9F", TFs="#E08214"), name=NULL) +
  labs(x=NULL, y="Count", title="(B) Cluster composition") +
  theme_bw(base_size=10) + theme(panel.grid.minor=element_blank(),
    axis.text.x=element_text(angle=45,hjust=1), legend.position="top")

# ── panel C: conservation pies (as in make_Fig3C.R) ──────────────────────────
og <- read.delim(Sys.glob("data/*/Orthogroups_20240613.tsv")[1], stringsAsFactors=FALSE, check.names=FALSE)
cols<-colnames(og); Ar<-"Arabidopsis_thaliana"; Bi<-"Betula_pendula"
Po<-grep("Populus_trem",cols,value=TRUE)[1]; Pa<-grep("Picea_abies",cols,value=TRUE)[1]; Ps<-grep("Pinus_sylv",cols,value=TRUE)[1]
splt<-function(x) unlist(strsplit(x,"[, ]+")); g2og<-new.env(hash=TRUE)
for(i in seq_len(nrow(og))) for(g in c(splt(og[i,Ar]),splt(og[i,Bi]),splt(og[i,Po]),splt(og[i,Pa]),splt(og[i,Ps])))
  if(nzchar(g)) assign(g,i,envir=g2og)
clade<-function(g) if(grepl("^AT",g)) "A" else if(grepl("^Bpev|^Potr",g)) "D" else "C"
mapping<-read.csv("data/superclusters/SC_mapping.csv",stringsAsFactors=FALSE)
CAT<-c("A","C","D","AC","AD","CD","ACD")
NAMES<-c(A="Arabidopsis only",C="Coniferous only",D="Deciduous only",AC="Arabidopsis & Coniferous",
         AD="Arabidopsis & Deciduous",CD="Coniferous & Deciduous",ACD="All three groups")
pie_of<-function(orig,tfonly){ sc<-sub("\\.\\d+$","",read.csv(Sys.glob(sprintf("data/superclusters/DEGs_TFs/%s/%s_gene_list.csv",orig,orig))[1],stringsAsFactors=FALSE)[,2])
  if(tfonly) sc<-sc[sc%in%TFset]
  ogi<-sapply(sc,function(g){i<-mget(g,envir=g2og,ifnotfound=NA)[[1]]; if(is.null(i)||is.na(i)) NA else i})
  cl<-sapply(sc,clade); ok<-!is.na(ogi); ogset<-tapply(cl[ok],ogi[ok],function(v) paste(sort(unique(v)),collapse=""))
  t<-table(factor(ogset,levels=CAT)); 100*t/sum(t) }
Crows<-list()
for(k in seq_len(nrow(mapping))) for(ty in c("DEGs","TF-DEGs")){
  v<-pie_of(mapping$Original_supercluster[k], ty=="TF-DEGs")
  Crows[[length(Crows)+1]]<-data.frame(SC=mapping$Main_text[k],type=ty,
    category=factor(NAMES[CAT],levels=NAMES[CAT]),pct=as.numeric(v)) }
Cd<-bind_rows(Crows); Cd$SC<-factor(Cd$SC,levels=names(SC)); Cd$type<-factor(Cd$type,levels=c("DEGs","TF-DEGs"))
pal<-c("Arabidopsis only"="#3F936A","Coniferous only"="#786C9C","Deciduous only"="#E49C24",
       "Arabidopsis & Coniferous"="#6CB490","Arabidopsis & Deciduous"="#E48460",
       "Coniferous & Deciduous"="#9090B4","All three groups"="#A8C048")
pC <- ggplot(Cd, aes(x="", y=pct, fill=category)) +
  geom_col(width=1,colour="white",linewidth=0.3) + coord_polar("y") +
  facet_grid(SC~type, switch="y") + scale_fill_manual(values=pal,name=NULL) +
  labs(title="(C) Conservation") +
  theme_void(base_size=10) + theme(legend.position="right", plot.title=element_text(hjust=0.5),
    strip.text=element_text(size=9))

fig <- (pA | pB | pC) + plot_layout(widths=c(1,1,1.4))
ggsave("Figure3.pdf", fig, width = 15, height = 9)
ggsave("Figure3.png", fig, width = 15, height = 9, dpi = 130)
cat("Panel B counts (should match Table S4 membership):\n")
print(B %>% group_by(SC,Species) %>% summarise(total=sum(n),.groups="drop") %>% pivot_wider(names_from=Species,values_from=total))
