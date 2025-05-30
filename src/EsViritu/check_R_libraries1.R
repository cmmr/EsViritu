#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(if (!require(reactable)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(htmltools)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(dplyr)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(reactablefmtr)) quit(status=1)))
## not requiring dataui
#suppressMessages(suppressWarnings(if (!require(dataui)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(data.table)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(RColorBrewer)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(viridis)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(scales)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(knitr)) quit(status=1)))
suppressMessages(suppressWarnings(if (!require(readxl)) quit(status=1)))