#!/usr/bin/env Rscript

### makes a summary reactable for any number of samples (within reason) based on a summary info table + coverage table
### input summary files generated with: make_summary_batch_of_samples1.sh
suppressMessages(library(reactable))
suppressMessages(library(htmltools))
suppressMessages(library(dplyr))
suppressMessages(library(reactablefmtr))
suppressMessages(library(dataui))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
suppressMessages(library(scales))
suppressMessages(library(knitr))
suppressMessages(library(readxl))
suppressMessages(library(lubridate))
suppressMessages(library(stringr))

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Three arguments must be supplied (.combined.mean_cov.tsv, .detected_virus.combined.tax.tsv, then outname).", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  sprintf("arguments found. Running.")
}


coverage_data <- fread(args[1], sep = "\t", header = FALSE, col.names = c("sample_ID", "accession", "start_base", "end_base", "mean_depth"))

genome_data <- fread(args[2], sep = "\t", header = TRUE)


### format coverage data for reactable

coverage_data$mean_depth <- round(coverage_data$mean_depth)


sum_coverage <- coverage_data %>% 
  group_by(sample_ID, accession) %>%
  summarize(coverage = list(mean_depth)) 



### merge read depth and tax/abundance info on sample_ID and accession columns
#colnames(genome_data)
setkey(genome_data, sample_ID, accession)

sum_coverage <- setDT(sum_coverage)
#colnames(sum_coverage)
setkey(sum_coverage, sample_ID, accession)

combined_data <- genome_data[sum_coverage]
combined_data <- na.omit(combined_data)


combined_data$Percent_covered <- combined_data$covered_bases / combined_data$reference_length

combined_data <- subset(combined_data, 
                        select = c("sample_ID", "sequence_name", "accession", 
                                   "reference_length", "Percent_covered", 
                                   "RPKMF", "reads_aligned", "genus", "species", 
                                   "subspecies", "coverage"))

nice_table <- combined_data %>%
  reactable(
    .,
    pagination = TRUE,
    filterable = TRUE,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(10, 20, 100),
    defaultPageSize = 10,
    columns = list(
      Percent_covered = colDef(
        cell = data_bars(
          data = .,
          fill_color = viridis::magma(5, direction = -1),
          background = '#F1F1F1',
          min_value = 0,
          max_value = 1,
          round_edges = TRUE,
          text_position = 'outside-end',
          number_fmt = scales::percent
        )
      ),
      RPKMF = colDef(
        maxWidth = 100, 
        style = color_scales(., colors = c("grey", "gold", "maroon"), bias = 2), 
        format = colFormat(digits = 1)
      ),
      coverage = colDef(
        filterable = FALSE,
        width = 500,
        cell = react_sparkline(
          .,
          height = 50,
          decimals = 1,
          show_area = TRUE,
          area_color = "darkgreen",
          line_curve = "cardinal",
          highlight_points = highlight_points(max = "blue"),
          labels = c("max"),
          statline = "min",
          statline_color = "black"
          #bandline = "innerquartiles",
          #bandline_color = "darkgreen"
        )))) %>% 
  add_title(sprintf("%s Batch Detected Virus Summary", args[3])) %>%
  google_font(font_family = "Oswald")

nice_table %>% save_reactable_test(sprintf("%s.batch_detected_viruses.html", args[3]))

#write.table(genome_data, file = sprintf("%s.batch_detected_viruses.info.tsv", args[3]), quote = F, row.names = F, sep = "\t")


