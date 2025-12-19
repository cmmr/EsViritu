#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(reactable)))
suppressMessages(suppressWarnings(library(htmltools)))
suppressMessages(suppressWarnings(library(reactablefmtr)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(magrittr)))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Five arguments must be supplied:
    coverage windows tsv, 
    then main table tsv, 
    output directory, 
    sample_ID,
    reads_in_sample",
    call. = FALSE
  )
}
coverage_data <- read.table(
  args[1],
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
genome_data <- read.table(
  args[2],
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
coverage_data$average_coverage <- ceiling(coverage_data$average_coverage)
sum_coverage <- aggregate(
  average_coverage ~ Accession,
  data = coverage_data,
  FUN = function(x) list(x)
)
names(sum_coverage)[2] <- "coverage"
combined_data <- merge(genome_data, sum_coverage, by = "Accession")
combined_data$Percent_covered <- combined_data$covered_bases / combined_data$Length
keep <- c(
  "Name", "Accession", "Segment", "Assembly",
  "Length", "Percent_covered", "RPKMF",
  "read_count", "avg_read_identity", "Pi", "genus",
  "species", "subspecies", "coverage"
)
combined_data <- combined_data[, keep]
combined_data$genus <- sub("^g__", "", combined_data$genus)
combined_data$species <- sub("^s__", "", combined_data$species)
combined_data$subspecies <- sub("^t__", "", combined_data$subspecies)

magma_colors <- c("#F6E620", "#F98E09", "#E14E0E", "#8B0A50", "#000004")

safe_google_font <- function(tbl, font_family = "Oswald") {
  font_url <- sprintf(
    "https://fonts.googleapis.com/css2?family=%s&display=swap",
    utils::URLencode(font_family, reserved = TRUE)
  )

  can_fetch <- FALSE
  if (requireNamespace("curl", quietly = TRUE)) {
    can_fetch <- tryCatch({
      handle <- curl::new_handle(nobody = TRUE, connecttimeout = 2, timeout = 5)
      res <- curl::curl_fetch_memory(font_url, handle = handle)
      res$status_code >= 200 && res$status_code < 400
    }, error = function(e) FALSE)
  }

  if (!can_fetch) {
    warning(
      sprintf(
        "Skipping google_font('%s') because %s is unreachable.",
        font_family,
        font_url
      ),
      call. = FALSE
    )
    return(tbl)
  }

  tryCatch(
    google_font(tbl, font_family = font_family),
    error = function(e) {
      warning(
        sprintf(
          "Failed to apply google_font('%s'): %s. Continuing without custom font.",
          font_family,
          conditionMessage(e)
        ),
        call. = FALSE
      )
      tbl
    }
  )
}

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
    add_title(sprintf("%s EsViritu Detected virus Summary", args[4])) %>%
    add_subtitle(
      sprintf(
        "Generated at %s | %s filtered reads in sample",
        format(Sys.time(), "%Y-%m-%d %H:%M"), args[5]
      )
    ) %>%
    safe_google_font(font_family = "Oswald")
  
} else {
  keep_nodataui <- c(
    "Name", "Accession", "Segment", "Assembly",
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
    add_title(sprintf("%s EsViritu Detected virus Summary", args[4])) %>%
    add_subtitle(
      sprintf(
        "Generated at %s | %s filtered reads in sample",
        format(Sys.time(), "%Y-%m-%d %H:%M"), args[5]
      )
    ) %>%
    safe_google_font(font_family = "Oswald")
  
}

nice_table %>% save_reactable_test(
  sprintf("%s/%s_EsViritu_reactable.html", args[3], args[4])
)
