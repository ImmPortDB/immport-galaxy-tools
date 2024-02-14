#!/usr/bin/env Rscript
######################################################################
#                  Copyright (c) 2018 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1 - January 2018
# Author: Cristel Thomas
#
#

library(MetaCyto)

check_panel <- function(df, outfile = "") {
  report <- panelSummary(df, ".", cluster = F, width = 30, height = 20)
  markers <- data.frame("Markers" = row.names(report))
  s <- cbind(markers, report)
  write.table(s, file = outfile, quote = F, row.names = F, col.names = T, sep = "\t")
}

check_input_format <- function(infile = "", outfile = "") {
  df <- read.table(infile, sep = "\t", header = T, colClasses = "character")
  nm <- colnames(df)
  check_ab <- if ("antibodies" %in% nm) TRUE else FALSE
  check_sdy <- if ("study_id" %in% nm) TRUE else FALSE

  if (check_sdy && check_ab) {
    check_panel(df, outfile)
  } else {
    quit(save = "no", status = 10, runLast = FALSE)
  }
}

################################################################################
args <- commandArgs(trailingOnly = TRUE)

check_input_format(args[1], args[2])
