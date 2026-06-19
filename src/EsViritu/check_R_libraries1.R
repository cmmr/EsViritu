#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(if (!require(reactable)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(htmltools)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(reactablefmtr)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(scales)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(magrittr)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(base64enc)) quit(status=1)))

## install vendored dataui if not already available
if (!requireNamespace("dataui", quietly = TRUE)) {
  .this_script <- local({
    argv <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", argv, value = TRUE)
    if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]))
    else sys.frame(1)$ofile
  })
  vendor_path <- file.path(dirname(.this_script), "dataui_vendored")
  install.packages(vendor_path, repos = NULL, type = "source", quiet = TRUE)
}
suppressMessages(suppressWarnings(if (!require(dataui)) quit(status=1)))
