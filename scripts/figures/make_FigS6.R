## Figure S6 — expression profiles of co-expressologs (P<0.1) directly linked to the
## Arabidopsis CBFs (CBF1 AT4G25490, CBF2 AT4G25470, CBF3 AT4G25480), plotted as median VST.
## Reproducible from the deposited data: data/ComPlEx/co_expressologs.RData + data/expression/.
## Run from the deposit root (the dir that contains data/).
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(ggplot2); library(ggrepel)
})

load("data/ComPlEx/co_expressologs.RData")            # -> co_expressologs (P<0.1 set)
prefix <- c(col0="Col", ost0="Ost", aspen="Pt", birch="Bp", spruce="Pa", pine="Ps")
rd <- function(sp) read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp),
                              header=TRUE, sep="\t", check.names=FALSE)
E <- lapply(names(prefix), rd); names(E) <- names(prefix)

## timepoint coding (verified): 01=Control 02=6h 03=24h 04=3d 05=10d
tp_levels <- c("Control","6h","24h","3d","10d")
tp_map    <- c("01"="Control","02"="6h","03"="24h","04"="3d","05"="10d")

## sanity: per-code sample counts per species (should match the Methods 25/24/25/25/15/18)
cat("== per-code sample counts (Control/6h/24h/3d/10d) ==\n")
for(sp in names(E)){
  cols <- setdiff(colnames(E[[sp]]), "Genes")
  code <- sub(sprintf("^%s[0-9]+\\.", prefix[sp]), "", cols)
  cnt <- table(factor(tp_map[code], levels=tp_levels))
  cat(sprintf("  %-7s %s  (total %d)\n", sp, paste(cnt, collapse="/"), sum(cnt)))
}

## co-expressologs directly linked to the col0 CBFs (verbatim hits() logic)
hits <- function(cbf){
  h <- co_expressologs %>% filter(GeneSpecies1==cbf, grepl("col0", Species1, ignore.case=TRUE)) %>%
         dplyr::select(GeneSpecies2, Species2)
  r <- co_expressologs %>% filter(GeneSpecies2==cbf, grepl("col0", Species1, ignore.case=TRUE),
                                  grepl("ost0", Species2, ignore.case=TRUE)) %>%
         dplyr::select(GeneSpecies1, Species1) %>%
         dplyr::rename(GeneSpecies2=GeneSpecies1, Species2=Species1)
  rbind(h, r)
}
CBFs <- c("AT4G25490","AT4G25470","AT4G25480")
HIT  <- bind_rows(lapply(CBFs, hits)) %>% distinct()

## build median-VST-per-timepoint long table for one species panel
panel <- function(sp){
  genes <- HIT %>% filter(Species2==sp) %>% pull(GeneSpecies2) %>% unique()
  df <- E[[sp]] %>% filter(Genes %in% genes)
  if(nrow(df)==0) return(NULL)
  long <- df %>% pivot_longer(-Genes, names_to="sample", values_to="expr") %>%
    mutate(code = sub(sprintf("^%s[0-9]+\\.", prefix[sp]), "", sample),
           Timepoint = factor(tp_map[code], levels=tp_levels)) %>%
    filter(!is.na(Timepoint)) %>%
    group_by(Genes, Timepoint) %>% summarise(med = median(expr, na.rm=TRUE), .groups="drop") %>%
    mutate(Species = sp)
  long
}
species_disp <- c(col0="Col-0", ost0="Ost-0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")
combined <- bind_rows(lapply(names(prefix), panel))
combined$Species <- factor(species_disp[combined$Species],
                           levels=c("Col-0","Ost-0","Aspen","Birch","Spruce","Pine"))
## label CBFs
combined <- combined %>% mutate(gene = case_when(
  Genes=="AT4G25490" ~ "AT4G25490_CBF1", Genes=="AT4G25470" ~ "AT4G25470_CBF2",
  Genes=="AT4G25480" ~ "AT4G25480_CBF3", TRUE ~ Genes))

cat("\n== genes plotted per species panel ==\n")
print(combined %>% distinct(Species, gene) %>% count(Species, name="n_genes"))

lab <- combined %>% group_by(Species, gene) %>%
  filter(as.numeric(Timepoint)==max(as.numeric(Timepoint))) %>% slice_tail(n=1) %>% ungroup()

p <- ggplot(combined, aes(Timepoint, med, group=gene, color=gene)) +
  geom_line(linewidth=0.9) +
  geom_text_repel(data=lab, aes(label=gene), hjust=0, size=3.0, direction="y",
                  segment.color="grey70", box.padding=0.2, max.overlaps=Inf,
                  show.legend=FALSE, xlim=c(5.2, NA)) +
  facet_wrap(~Species, scales="free_y", ncol=3) +
  scale_x_discrete(expand=expansion(add=c(0.1, 3.2))) +
  coord_cartesian(clip="off") + theme_minimal() +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(size=14),
        legend.position="none", plot.margin=margin(t=5,r=45,b=5,l=5),
        panel.spacing=unit(1.4,"lines"), plot.title=element_text(size=15,hjust=0.5),
        axis.title=element_text(size=13), axis.text=element_text(size=11)) +
  labs(title="Expression profiles of co-expressologs (P<0.1) directly linked to Arabidopsis CBFs",
       x="Cold time course", y="Expression (median VST)")

ggsave("Figure_S6.pdf", p, width=13.5, height=9)
ggsave("Figure_S6_preview.png", p, width=13.5, height=9, dpi=120, bg="white")
cat("\nDONE — wrote Figure_S6.pdf + Figure_S6_preview.png\n")
