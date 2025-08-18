#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(reactable)))
suppressMessages(suppressWarnings(library(htmltools)))
suppressMessages(suppressWarnings(library(reactablefmtr)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(magrittr)))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Three arguments must be supplied:
    coverage windows tsv, 
    main table tsv, 
    project_ID",
    call. = FALSE
  )
}
coverage_data <- read.table(
  args[1],
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)
genome_data <- read.table(
  args[2],
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)
coverage_data$average_coverage <- ceiling(coverage_data$average_coverage)
sum_coverage <- aggregate(
  average_coverage ~ Accession + sample_ID,
  data = coverage_data,
  FUN = function(x) list(x)
)
names(sum_coverage)[3] <- "coverage"
combined_data <- merge(genome_data, sum_coverage, by = c("sample_ID", "Accession"))
combined_data$Percent_covered <- combined_data$covered_bases / combined_data$Length
keep <- c(
  "sample_ID", "Name", "Accession", "Segment", "Assembly",
  "Length", "Percent_covered", "RPKMF",
  "read_count", "avg_read_identity", "Pi", "genus",
  "species", "subspecies", "coverage"
)
combined_data <- combined_data[, keep]
combined_data$genus <- sub("^g__", "", combined_data$genus)
combined_data$species <- sub("^s__", "", combined_data$species)
combined_data$subspecies <- sub("^t__", "", combined_data$subspecies)

magma_colors <- c("#F6E620", "#F98E09", "#E14E0E", "#8B0A50", "#000004")

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
            fill_color = magma_colors,
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
          maxWidth = 150,
          style = color_scales(., colors = c("grey", "gold", "maroon"), bias = 2),
          format = colFormat(digits = 2)
        ),
        avg_read_identity = colDef(
          name = "Average Read Identity",
          maxWidth = 80,
          style = color_scales(., colors = c("#e09c9c", "#93adc8"), bias = 2),
          format = colFormat(percent = TRUE, digits = 1)
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
    add_title(sprintf("%s EsViritu Detected virus Summary", args[3])) %>%
    add_subtitle(
      sprintf(
        "Generated at %s",
        format(Sys.time(), "%Y-%m-%d %H:%M")
      )
    ) %>%
    google_font(font_family = "Oswald")
  
} else {
  keep_nodataui <- c(
    "sample_ID", "Name", "Accession", "Segment", "Assembly",
    "Length", "Percent_covered", "RPKMF",
    "read_count", "avg_read_identity", "Pi", "genus",
    "species", "subspecies"
  )
  nice_table <- combined_data[, keep_nodataui] %>%
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
            fill_color = magma_colors,
            background = '#F1F1F1',
            min_value = 0,
            max_value = 1,
            round_edges = TRUE,
            text_position = 'outside-end',
            number_fmt = scales::percent
          )
        ),
        
        RPKMF = colDef(
          maxWidth = 150,
          style = color_scales(., colors = c("grey", "gold", "maroon"), bias = 2), 
          format = colFormat(digits = 2)
        ),
        avg_read_identity = colDef(
          name = "Average Read Identity",
          maxWidth = 80,
          style = color_scales(., colors = c("#e09c9c", "#93adc8"), bias = 2),
          format = colFormat(percent = TRUE, digits = 1)
        )
      )) %>% 
    add_title(sprintf("%s EsViritu Detected virus Summary", args[3])) %>%
    add_subtitle(
      sprintf(
        "Generated at %s",
        format(Sys.time(), "%Y-%m-%d %H:%M")
      )
    ) %>%
    google_font(font_family = "Oswald")
  
}

nice_table %>% save_reactable_test(
  sprintf("%s_EsViritu_project_reactable.html", args[3])
)
