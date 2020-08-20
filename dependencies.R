### dependencies.R for multiEditR

###########################################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

# run these lines one at a time, and be careful to note if there are 
# any errors put out during installation

# CRAN packages
install.packages("shiny")
install.packages("magrittr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("plotly")
install.packages("grr")
install.packages("printr")
install.packages("plyr")
install.packages("readr")
install.packages("rmarkdown")
install.packages("readr")
install.packages("gamlss")
install.packages("BiocManager")
install.packages("DT")
install.packages("grDevices")

# Bioconductor packages
# Updated 4.7.19 due to error with Bioconductor packages
# Updated again 4.30.19
# Run this code chunk in terminal first
#
  # options(repos = BiocManager::repositories())
  # source("https://bioconductor.org/biocLite.R")
#

BiocManager::install("Biostrings")
BiocManager::install("sangerseqR")
