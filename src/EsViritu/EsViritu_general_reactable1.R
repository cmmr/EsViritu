#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(reactable)))
suppressMessages(suppressWarnings(library(htmltools)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(reactablefmtr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(viridis)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(knitr)))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Four arguments must be supplied:
    coverage windows tsv, 
    then main table tsv, 
    output directory, 
    sample_ID,
    reads_in_sample",
    call. = FALSE
  )
}

coverage_data <- fread(
  args[1],
  sep = "\t",
  header = TRUE,
)

genome_data <- fread(args[2], sep = "\t", header = TRUE)


#coverage_data
sum_coverage <- coverage_data %>%
  mutate(average_coverage = ceiling(average_coverage)) %>%
  group_by(Accession) %>%
  summarize(coverage = list(average_coverage))



combined_data <- merge(genome_data, sum_coverage, by = "Accession") %>%
  mutate(Percent_covered = covered_bases / Length) %>%
  select(c(
    "Name", "Accession", "Assembly", "Length",
    "Percent_covered", "RPKMF", "read_count",
    "genus", "species", "subspecies", "coverage"
  )
  ) %>%
  mutate(
    genus = gsub("^g__", "", genus),
    species = gsub("^s__", "", species),
    subspecies = gsub("^t__", "", subspecies)
  )

## check for dataui
is_dataui <- require(dataui)

if (is_dataui == TRUE) {
  suppressMessages(library(dataui))
  
  nice_table <- combined_data %>%
    reactable(
      .,
      pagination = TRUE,
      filterable = TRUE,
      showPageSizeOptions = TRUE,
      pageSizeOptions = c(10, 20, 100),
      defaultPageSize = 10,
      columns = list(
        read_count = colDef(
          name = "# of Aligned Reads"
        ),
        Percent_covered = colDef(
          name = "% Covered",
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
          name = "RPKMF\nReads per Kilobase\nper Million Filtered Reads",
          maxWidth = 100, 
          style = color_scales(., colors = c("grey", "gold", "maroon"), bias = 2), 
          format = colFormat(digits = 2)
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
          )))) %>% 
    add_title(sprintf("%s EsViritu Detected virus Summary", args[4])) %>%
    add_subtitle(
      sprintf(
        "Generated at %s | %s filtered reads in sample",
        format(Sys.time(), "%Y-%m-%d %H:%M"), args[5]
      )
    ) %>%
    google_font(font_family = "Oswald")
  
} else {
  nice_table <- combined_data %>%
    select(-coverage) %>%
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
          format = colFormat(digits = 2)
        )
      )) %>% 
    add_title(sprintf("%s EsViritu Detected virus Summary", args[4])) %>%
    add_subtitle(
      sprintf(
        "Generated at %s | %s filtered reads in sample",
        format(Sys.time(), "%Y-%m-%d %H:%M"), args[5]
      )
    ) %>%
    google_font(font_family = "Oswald")
  
}

nice_table %>% save_reactable_test(
  sprintf("%s/%s_EsViritu_reactable.html", args[3], args[4])
)
