#!/usr/bin/env Rscript

suppressMessages(if (!require(reactable)) quit(status=1))
suppressMessages(if (!require(htmltools)) quit(status=1))
suppressMessages(if (!require(dplyr)) quit(status=1))
suppressMessages(if (!require(reactablefmtr)) quit(status=1))
suppressMessages(if (!require(dataui)) quit(status=1))
suppressMessages(if (!require(data.table)) quit(status=1))
suppressMessages(if (!require(RColorBrewer)) quit(status=1))
suppressMessages(if (!require(viridis)) quit(status=1))
suppressMessages(if (!require(scales)) quit(status=1))
suppressMessages(if (!require(knitr)) quit(status=1))
suppressMessages(if (!require(readxl)) quit(status=1))