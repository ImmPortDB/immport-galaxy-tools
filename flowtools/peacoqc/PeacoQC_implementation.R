#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("PeacoQC")

library(PeacoQC)

# Determine channels on which quality control should be done
channels <- c(1, 3, 5:14, 18, 21)

# Remove margins
ff <- RemoveMargins(ff=ff, channels=channels, output="frame")

# Compensate and transform the data
ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
ff <- flowCore::transform(ff, 
                          flowCore::estimateLogicle(ff, colnames(flowCore::keyword(ff)$SPILL)))


#Recursively read-in FCS
for (file in files){
  ff <- flowCore::read.FCS(file)
  
  # Remove margins
  ff <- RemoveMargins(ff=ff, channels=channels, output="frame")
  
  # Compensate and transform the data
  ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
  ff <- flowCore::transform(ff, 
                            flowCore::estimateLogicle(ff, 
                                                      colnames(flowCore::keyword(ff)$SPILL)))
  peacoqc_res <- PeacoQC(
    ff,
    channels,
    determine_good_cells="all",
    IT_limit=0.6,
    save_fcs=T,
    plot=T)
}


# PeacoQCHeatmap

# Find the path to the report that was created by using the PeacoQC function
location <- system.file("extdata", "PeacoQC_report.txt", package="PeacoQC")

# Make heatmap overview of the quality control run
PeacoQCHeatmap(report_location=location, show_values = FALSE,
               show_row_names = FALSE)

# Make heatmap with only the runs of the last test
PeacoQCHeatmap(report_location=location, show_values = FALSE, 
               latest_tests=TRUE, show_row_names = FALSE)

# Make heatmap with row annotation
PeacoQCHeatmap(report_location=location, show_values = FALSE,
               show_row_names = FALSE,
               row_split=c(rep("r1",7), rep("r2", 55)))






