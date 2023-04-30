#!/usr/bin/env Rscript

### summarizes percent identity of consensus sequences to reference database
### top 3 hits




args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Three arguments must be supplied (.basecomp1.tsv, .blastn.tsv, output directory, samplename).", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  sprintf("arguments found. Running.")
}

library(dplyr)
library(data.table)

basecomp_dt <- fread(args[1], 
                     header = T, sep = "\t", 
                     col.names = c("consensusID", "length", "N", "ATGC"))


blastn_dt <- fread(args[2],
                   header = F, sep = "\t", 
                   col.names = c("consensusID", "refID", "pident", 
                                 "length", "mismatch", 
                                 "gapopen", "qstart", "qend", "sstart", "send", 
                                 "evalue", "bitscore", "nident")) %>%
  group_by(consensusID, refID) %>%
  summarize(total_length = sum(length),
            total_mismatch = sum(mismatch),
            total_identical = sum(nident),
            total_bitscore = sum(bitscore)) %>%
  ungroup()

merge_dt <- merge(blastn_dt, basecomp_dt, by = "consensusID") %>%
  mutate(perc_identical_ATGC = total_identical / ATGC) %>% 
  group_by(consensusID) %>%
  top_n(perc_identical_ATGC, n=3) %>%
  select(c(consensusID, refID, perc_identical_ATGC)) %>%
  mutate(hit_rank = round(rank(-perc_identical_ATGC))) %>%
  arrange(consensusID, desc(perc_identical_ATGC), refID)

write.table(merge_dt, 
            file = sprintf("%s/%s_consensus_seqs_vs_ref_seqs.tsv", args[3], args[4]),
            quote = F, row.names = F, sep = "\t")


