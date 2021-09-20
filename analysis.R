## ------- get-libs --------

library(DESeq2)
library(tidyverse)
library(airway)
library(ggheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gt)

## ------- show-dataset -----

data("gse", package = "airway")

gse

## ------- cite-workflow ----

#citation("rnaseqGene")

## ------- run-deseq2 -----

dds <-  DESeqDataSet(gse, design = ~ donor + condition)

colnames(dds) <- paste(dds$donor,str_trunc(as.character(dds$condition),width = 3,ellipsis = ""), sep = "_")

keep <- rowSums(counts(dds)) > 1

dds <- dds[keep,]

vsd <- vst(dds, blind = FALSE)

sampleDists <- as.dist(1-cor(assay(vsd), method="spearman")) %>% as.matrix()

dds <- DESeq(dds)

res <- results(dds,tidy = T, alpha = 0.05, lfcThreshold = 1) %>% as_tibble()

top_gene_list <- filter(res,padj < 0.05) %>%
  group_by(dir=ifelse(sign(log2FoldChange)==-1,"down","up")) %>%
  slice_max(abs(log2FoldChange), n=100) %>%
  mutate(row = str_remove(row,"\\.\\d+")) %>%
  dplyr::select(dir,row) %>%
  split(.,.$dir) %>%
  map(pull,row)


ego_up <- enrichGO(gene = top_gene_list$up,
                   keyType = "ENSEMBL",
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 25,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05)

## ------- make-plots -----

g_sample_all_gene_hm <- ggheatmap(sampleDists,scale = "none",dist_method = NA,text_show_cols = F,tree_height_rows = T,
                                  legendName = "1 - Spearman's rho", text_position_rows = "left") +
  scale_fill_distiller(direction = -1) +
  theme(aspect.ratio = 1, legend.position = "bottom") +
  guides(fill = guide_legend(title.position = "bottom",
                             title.hjust = 0.5,
                             keyheight = rel(0.5),
                             keywidth = rel(0.75),
                            title.theme = element_text(face="italic"))) +
  theme(legend.title = element_text(size=rel(0.5)))


g_top_2_genes <- filter(res,padj < 0.05) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_min(padj, n=1) %>%
  slice_max(abs(log2FoldChange), n=1) %>%
  pull(row) %>%
  map_df(~tibble(plotCounts(dds, gene=.x,returnData=T),gene=.x)) %>%
  as_tibble() %>%
  ggplot(aes(condition,count)) +
  geom_jitter(width = 0.3) +
  facet_wrap(~gene, scales="free") +
  ylab("expression") +
  theme_classic() +
  xlab("") + theme(aspect.ratio = 1.5, 
                     strip.text.x = element_text(size=rel(0.5)),
                     axis.text.x = element_text(size=rel(0.5), angle=30, hjust=1))

g_volc <- ggplot(res,aes(log2FoldChange,-log10(pvalue))) +
  geom_point(data= . %>% filter(padj >= 0.05 | abs(log2FoldChange) <= 0.5), color="lightgray", alpha=0.5,size=rel(0.5)) +
  geom_point(data= . %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5), aes(color=log2FoldChange)) +
  coord_cartesian(xlim=c(-11,11)) +
  theme_classic() +
  scale_color_distiller(type="div",palette = 5) +
  guides(color="none")

grid_top10_heat <- filter(res,padj < 0.05) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_max(abs(log2FoldChange), n=10) %>%
  pull(row) %>%
  assay(vsd)[.,] %>%
  pheatmap::pheatmap(cellwidth = 12, cellheight = 12, fontsize = 8,treeheight_col = 8, treeheight_row = 0)
  
  #ggheatmap(show_cluster_cols = T,
  #          cluster_cols = T,
  #          cluster_rows=T,
  #          show_cluster_rows = T)

g_enrich_corticsteroid_up <- clusterProfiler::dotplot(object = ego_up) + coord_flip() + scale_size_continuous(breaks = c(5,7,9)) +
  theme(legend.title = element_text(size=rel(0.5)),legend.box = "horizontal",
        legend.text = element_text(size=rel(0.5)), legend.justification = "center",
        legend.key = element_rect(size=rel(0.2)), legend.key.size = unit(0.1,"in"),
        axis.text.x=element_text(angle=30, hjust=1, size=rel(0.5)))

gt_sample_table <- colData(dds) %>%
  as_tibble(rownames = "sample") %>%
  dplyr::rename(accession = "names") %>%
  gt::gt() %>%
  bstfun::as_ggplot()
