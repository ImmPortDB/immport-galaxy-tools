#!/usr/bin/env Rscript
# Module for Galaxy
# Compares groups of FCS to FlowSOM reference tree
# with FlowSOM
######################################################################
#                  Copyright (c) 2017 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1
# Cristel Thomas
#
#
library(FlowSOM)
library(flowCore)

check_files <- function(groups) {
  all_files <- unlist(groups)
  all_unique <- unique(all_files)
  if (length(all_unique) != length(all_files)) {
    quit(save = "no", status = 14, runLast = FALSE)
  }
}

compare_lists <- function(m1, m2) {
  list_check <- TRUE
  if (is.na(all(m1 == m2))) {
    mm1 <- is.na(m1)
    mm2 <- is.na(m2)
    if (all(mm1 == mm2)) {
      if (!all(m1 == m2, na.rm = TRUE)) {
        list_check <- FALSE
      }
    } else {
      list_check <- FALSE
    }
  } else if (!all(m1 == m2)) {
    list_check <- FALSE
  }
  return(list_check)
}

pretty_marker_names <- function(flow_frame) {
  n <- flow_frame@parameters@data[, "name"]
  d <- flow_frame@parameters@data[, "desc"]
  d[is.na(d)] <- n[is.na(d)]
  pretty_names <- list()
  if (any(grepl("#", d))) {
    # Support for hashtag notation:
    pretty_names <- gsub("#(.*)$", " (\\1)", d)
  } else {
    pretty_names <- paste(d, " <", n, ">", sep = "")
  }
  return(pretty_names)
}

compare_to_tree <- function(fst,
                            wilc_thresh = 0.05, output = "", plot = "", stats,
                            comp_groups, filenames) {
  group_res <- CountGroups(fst, groups = comp_groups, plot = FALSE)
  pdf(plot, useDingbats = FALSE, onefile = TRUE)
  tresh <- wilc_thresh
  pg <- PlotGroups(fst, group_res, p_tresh = tresh)
  dev.off()

  nb_nodes <- length(pg[[1]])
  nb_comp <- length(pg)
  m <- matrix(0, nrow = nb_nodes, ncol = nb_comp + 1)
  s <- seq_len(nb_nodes)
  m[, 1] <- as.character(s)
  for (i in 1:nb_comp) {
    m[s, i + 1] <- as.character(pg[[i]])
  }
  groupnames <- attr(comp_groups, "names")
  out_colnames <- paste(groupnames, collapse = "-")
  colnames(m) <- c("Node", out_colnames)
  write.table(m, file = output, quote = FALSE, row.names = FALSE,
              col.names = TRUE, sep = "\t",
              append = FALSE)

  ## get filenames
  filepaths <- unlist(comp_groups)
  fnames <- unlist(filenames)
  nb_files <- length(filepaths)
  comp_files <- list()
  for (i in seq_along(filepaths)) {
    comp_files[[filepaths[[i]]]] <- fnames[[i]]
  }

  group_list <- list()
  for (grp in attr(comp_groups, "names")) {
    for (f in comp_groups[[grp]]) {
      group_list[[f]] <- grp
    }
  }
  out_stats <- attr(stats, "names")
  if ("counts" %in% out_stats) {
    gp_counts <- as.matrix(group_res$counts)
    tpc <- matrix("", nrow = nb_files, ncol = 2)
    tpc[, 1] <- as.character(
                             lapply(rownames(gp_counts),
                                    function(x) comp_files[[x]]))
    tpc[, 2] <- as.character(
                             lapply(rownames(gp_counts),
                                    function(x) group_list[[x]]))
    gp_counts <- cbind(tpc, gp_counts)
    colnames(gp_counts)[[1]] <- "Filename"
    colnames(gp_counts)[[2]] <- "Group"
    t_gp_counts <- t(gp_counts)
    write.table(t_gp_counts,
                file = stats[["counts"]],
                quote = FALSE,
                row.names = TRUE,
                col.names = FALSE,
                sep = "\t",
                append = FALSE)
  }
  if ("pctgs" %in% out_stats) {
    gp_prop <- as.matrix(group_res$pctgs)
    tpp <- matrix("", nrow = nb_files, ncol = 2)
    tpp[, 1] <- as.character(
                             lapply(rownames(gp_prop),
                                    function(x) comp_files[[x]]))
    tpp[, 2] <- as.character(
                             lapply(rownames(gp_prop),
                                    function(x) group_list[[x]]))
    gp_prop <- cbind(tpp, gp_prop)
    colnames(gp_prop)[[1]] <- "Filename"
    colnames(gp_prop)[[2]] <- "Group"
    t_gp_prop <- t(gp_prop)
    write.table(t_gp_prop,
                file = stats[["pctgs"]],
                quote = FALSE,
                row.names = TRUE,
                col.names = FALSE,
                sep = "\t",
                append = FALSE)
  }
  if ("means" %in% out_stats) {
    gp_mean <- as.matrix(group_res$means)
    t_gp_mean <- t(gp_mean)
    tpm <- matrix(0, nrow = nb_nodes, ncol = 1)
    tpm[, 1] <- seq_len(nb_nodes)
    t_gp_mean <- cbind(tpm, t_gp_mean)
    colnames(t_gp_mean)[[1]] <- "Nodes"
    write.table(t_gp_mean,
                file = stats[["means"]],
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t",
                append = FALSE)
  }
  if ("medians" %in% out_stats) {
    gp_med <- as.matrix(group_res$medians)
    t_gp_med <- t(gp_med)
    tpd <- matrix(0, nrow = nb_nodes, ncol = 1)
    tpd[, 1] <- seq_len(nb_nodes)
    t_gp_med <- cbind(tpd, t_gp_med)
    colnames(t_gp_med)[[1]] <- "Nodes"
    write.table(t_gp_med,
                file = stats[["medians"]],
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t",
                append = FALSE)
  }
}

check_fcs <- function(tree,
                      output = "", plot = "", thresh = 0.05, stats, groups,
                      filenames) {

  fcsfiles <- unlist(groups)
  tree_valid <- FALSE
  marker_check <- TRUE
  tryCatch({
    fsomtree <- readRDS(tree)
    tree_valid <- TRUE
  }, error = function(ex) {
    print(paste(ex))
  })

  fst <- if (length(fsomtree) == 2) fsomtree[[1]] else fsomtree

  if (tree_valid) {
    tree_markers <- as.vector(fst$prettyColnames)
    if (length(tree_markers) < 1) {
      quit(save = "no", status = 11, runLast = FALSE)
    }
  } else {
    quit(save = "no", status = 11, runLast = FALSE)
  }

  for (i in seq_along(fcsfiles)) {
    tryCatch({
      fcs <- read.FCS(fcsfiles[i], transformation = FALSE)
    }, error = function(ex) {
      print(paste(ex))
    })
    if (i  ==  1) {
      m1 <- as.vector(pData(parameters(fcs))$desc)
      c1 <- colnames(fcs)
      # compare to tree markers
      pm <- pretty_marker_names(fcs)
      if (!all(tree_markers %in% pm)) {
        quit(save = "no", status = 13, runLast = FALSE)
      }
    } else {
      m2 <- as.vector(pData(parameters(fcs))$desc)
      c2 <- colnames(fcs)
      marker_check <- compare_lists(m1, m2)
      marker_channel <- compare_lists(c1, c2)
    }
  }
  if (marker_check && marker_channel) {
    compare_to_tree(fst, thresh, output, plot, stats, groups, filenames)
  } else {
    quit(save = "no", status = 12, runLast = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)

first_g1 <- 5
tot_args <- length(args)
g <- list()
tmplist <- c("counts", "means", "medians", "pctgs")

for (i in 5:13) {
  if (args[i] %in% tmplist) {
    first_g1 <- first_g1 + 2
    g[[args[i]]] <- args[i + 1]
  }
}

tmpargs <- paste(args[first_g1:tot_args], collapse = "=%=")
tmpgroups <- strsplit(tmpargs, "=%=DONE=%=")

groups <- list()
filenames <- list()
for (gps in tmpgroups[[1]]) {
  tmpgroup <- strsplit(gps, "=%=")
  nb_files <- (length(tmpgroup[[1]]) - 1) / 2
  tmplist <- character(nb_files)
  tmpnames <- character(nb_files)
  j <- 1
  for (i in 2:length(tmpgroup[[1]])) {
    if (!i %% 2) {
      tmplist[[j]] <- tmpgroup[[1]][i]
      tmpnames[[j]] <- tmpgroup[[1]][i + 1]
      j <- j + 1
    }
  }
  groups[[tmpgroup[[1]][1]]] <- tmplist
  filenames[[tmpgroup[[1]][1]]] <- tmpnames
}

check_files(groups)
check_fcs(args[1], args[2], args[3], args[4], g, groups, filenames)
