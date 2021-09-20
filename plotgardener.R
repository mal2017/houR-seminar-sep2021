source("analysis.R")

## ----------- create-pg-page -------------------

library(plotgardener)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationHub)

{{ pageCreate(width = 7, height = 4.25, default.units = "inches") }}

## ----------- create-plot-lab -------------------

plotgardener::pageGuideShow()

plotText(label = "A", fontsize = 12,
         x = 0.125, y = 0.125, just = "left", default.units = "inches")

## ----------- create-plot-a -------------------

a <- plotGG(g_sample_all_gene_hm,x = 0.5,y = 0.5,width = 3, height = 3)
