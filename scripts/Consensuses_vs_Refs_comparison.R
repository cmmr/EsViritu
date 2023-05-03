#!/usr/bin/env Rscript

### summarizes percent identity of consensus sequences to reference database
### top 5 hits




args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Four arguments must be supplied (virus_pathogen_database.all_metadata.tsv, .consensus.remove_NNs.VS.refs.anicalc.tsv, output directory, samplename).", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  sprintf("arguments found. Running.")
}

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

ani_removeNNs_dt <- fread(args[2], 
                          header = T, sep = "\t")

sum_ani_removeNNs_dt <- ani_removeNNs_dt %>%
  group_by(qname) %>%
  filter(qcov >= 70) %>%
  top_n((qcov * pid), n=5) %>%
  mutate(hit_rank = round(rank(-(qcov * pid)))) %>%
  arrange(qname, hit_rank)


seq_meta_dt <- fread(args[1], 
                     header = T, sep = "\t") %>%
  select(c(accession, Sequence_name, genus, species))


sum_meta_dt <- merge(sum_ani_removeNNs_dt, seq_meta_dt, 
                     by.x = "tname", by.y = "accession") %>%
  arrange(qname, hit_rank)

sum_meta_dt$sample_ID <- as.character(args[4])

setcolorder(sum_meta_dt, "qname")

write.table(sum_meta_dt, 
            file = sprintf("%s/%s_consensus_seqs_vs_ref_seqs.tsv", args[3], args[4]),
            quote = F, row.names = F, sep = "\t")


