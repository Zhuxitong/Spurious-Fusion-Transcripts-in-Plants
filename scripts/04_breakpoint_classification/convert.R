library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(readr)

geneLoc <- import.bed("Omh63.geneLocus.bed")
names(geneLoc) <- geneLoc$name

genomeGtf <- import.gff("Omh63.gtf", format = "gtf")
fusionMap <- import.bed("mh63_flagleaf.pc.breakpoints.genomic.bed")

fusionMap_trans <- mapToTranscripts(x = fusionMap, transcripts = geneLoc)

fusionMap %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, name) %>%
  add_column(order = 1:length(fusionMap)) %>%
  merge(
    x = fusionMap_trans %>% as.data.frame() %>% dplyr::select(seqnames, start, end, strand, xHits),
    y = .,
    by.x = "xHits",
    by.y = "order"
  ) %>%
  write_tsv(file = "mh63_flagleaf.pc.breakpoints.trans.txt")