#!/usr/bin/env Rscript

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

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Four arguments must be supplied (coverage_bed.tsv, then coverm.threshold.info.tsv, output directory, then sample ID).n", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  sprintf("arguments found. Running.")
}

coverage_data <- fread(args[1], sep = "\t", header = FALSE, col.names = c("accession", "start_base", "end_base", "mean_depth"))

genome_data <- fread(args[2], sep = "\t", header = TRUE)


#coverage_data

coverage_data$mean_depth <- round(coverage_data$mean_depth)

sum_coverage <- coverage_data %>% 
  #group_by(accession, grp = as.integer(gl(n(), 100, n()))) %>% 
  #summarise(mean = mean(coverage), n = n())
  group_by(accession) %>%
  summarize(coverage = list(mean_depth)) 



combined_data <- merge(genome_data, sum_coverage, by = "accession")
combined_data$Percent_covered <- combined_data$covered_bases / combined_data$reference_length

combined_data <- subset(combined_data, select = c("sequence_name", "accession", "reference_length", "Percent_covered", "RPKMF", "reads_aligned", "genus", "species", "subspecies", "coverage"))

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
  add_title(sprintf("%s Detected virus Summary", args[4])) %>%
  google_font(font_family = "Oswald")

nice_table %>% save_reactable_test(sprintf("%s/%s_EsViritu_reactable.html", args[3], args[4]))
