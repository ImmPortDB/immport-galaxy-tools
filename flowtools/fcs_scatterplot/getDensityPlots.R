#
# Density Plot Module for Galaxy
# FlowDensity
# 
# Version 1
# Cristel Thomas
#
#

library(flowCore)
library(flowDensity)

generateGraph <- function(input, channels, output, plot_default) {
  fcs <- read.FCS(input, transformation=F)
  ## marker names
  markers <- colnames(fcs)
  print_markers <- as.vector(pData(parameters(fcs))$desc)
  # Update print_markers if the $P?S not in the FCS file
  for (i in 1:length(print_markers)) {
    if (is.na(print_markers[i])) {
      print_markers[i] <- markers[i]
    }
  }

  if (plot_default) {
    channels <- c(grep(colnames(fcs), pattern="FSC"), 
                  grep(colnames(fcs), pattern="SSC"))
    if (length(channels) > 2) {
      #get first FSC and corresponding SSC
      channels <- c(grep(colnames(fcs), pattern="FSC-A"), 
                    grep(colnames(fcs), pattern="SSC-A"))
      if (length(channels) == 0) {
        channels <- c(grep(colnames(fcs), pattern="FSC-H"), 
                      grep(colnames(fcs), pattern="SSC-H"))
        if (length(channels) == 0) {
          channels <- c(grep(colnames(fcs), pattern="FSC-W"), 
                        grep(colnames(fcs), pattern="SSC-W"))
        }
      }
    }
    if (length(channels) == 0) {
      warning('No forward/side scatter channels found, no plots will be generated.')
	  quit(save = "no", status = 10, runLast = FALSE)
    }
  }

  nb_markers <- length(channels)
  pdf(output, useDingbats=FALSE, onefile=TRUE) 
  par(mfrow=c(2,2))
  for (m in 1:(nb_markers - 1)) {
    for (n in (m+1):nb_markers) {
      plotDens(fcs, c(channels[m],channels[n]), xlab = print_markers[channels[m]], ylab = print_markers[channels[n]])
    }
  }
  dev.off()
}

checkFCS <- function(input_file, channels, output_file, plot_default) {
  isValid <- F
  # Check file beginning matches FCS standard
  tryCatch({
    isValid = isFCSfile(input_file)
  }, error = function(ex) {
    print (paste("    ! Error in isFCSfile", ex))
  })

  if (isValid) {
    generateGraph(input_file, channels, output_file, plot_default)
  } else {
    print (paste(input_file, "does not meet FCS standard"))
  }
}

args <- commandArgs(trailingOnly = TRUE)
channels <- ""
flag_default <- FALSE

if (args[3]=="None") {
  flag_default <- TRUE
} else {
  if (args[3] == "i.e.:1,3,4"){
  	flag_default <- TRUE
  } else {
    channels <- as.numeric(strsplit(args[3], ",")[[1]])
    for (channel in channels){
	  if (is.na(channel)){
	    quit(save = "no", status = 11, runLast = FALSE)
	  } 
    }
	if (length(channels) == 1){
	  warning('Please indicate more than one marker to plot.')
	  quit(save = "no", status = 10, runLast = FALSE)
	}
  }
}

checkFCS(args[2], channels, args[4], flag_default)
