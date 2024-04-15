### global.R for multiEditR

###########################################################################################
# Copyright (C) 2020-2021 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

### Options
# turn on universal error message:
# "Error: An error has occurred. Check your logs or contact the app author for clarification."
options(shiny.sanitize.errors = TRUE)

### Libraries
library(sangerseqR)
library(tidyverse)
library(magrittr)
library(plyr)
library(gamlss)
library(readr)
library(shiny)
library(plotly)
library(shinythemes)
library(berryFunctions)
library(shinycssloaders)
library(markdown)
library(multiEditR)


### Data
bases = c("A", "C", "G", "T")
ACGT = bases

phred_scores = data.frame(stringsAsFactors=FALSE,
                          phred = c("!", "“", "$", "%", "&", "‘", "(", ")", "*", "+", ",", "–",
                                    ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                    ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F",
                                    "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                                    "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_",
                                    "`", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
                                    "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y",
                                    "z", "{", "|", "}", "~"),
                          prob = c(1, 0.794328235, 0.501187234, 0.398107171, 0.316227766,
                                   0.251188643, 0.199526232, 0.158489319, 0.125892541, 0.1,
                                   0.079432824, 0.063095734, 0.050118723, 0.039810717, 0.031622777,
                                   0.025118864, 0.019952623, 0.015848932, 0.012589254, 0.01,
                                   0.007943282, 0.006309573, 0.005011872, 0.003981072, 0.003162278,
                                   0.002511886, 0.001995262, 0.001584893, 0.001258925, 0.001, 0.000794328,
                                   0.000630957, 0.000501187, 0.000398107, 0.000316228, 0.000251189,
                                   0.000199526, 0.000158489, 0.000125893, 1e-04, 7.94328e-05,
                                   6.30957e-05, 5.01187e-05, 3.98107e-05, 3.16228e-05, 2.51189e-05,
                                   1.99526e-05, 1.58489e-05, 1.25893e-05, 1e-05, 7.9433e-06,
                                   6.3096e-06, 5.0119e-06, 3.9811e-06, 3.1623e-06, 2.5119e-06, 1.9953e-06,
                                   1.5849e-06, 1.2589e-06, 1e-06, 7.943e-07, 6.31e-07, 5.012e-07,
                                   3.981e-07, 3.162e-07, 2.512e-07, 1.995e-07, 1.585e-07, 1.259e-07,
                                   1e-07, 7.94e-08, 6.31e-08, 5.01e-08, 3.98e-08, 3.16e-08,
                                   2.51e-08, 2e-08, 1.58e-08, 1.26e-08, 1e-08, 7.9e-09, 6.3e-09, 5e-09,
                                   4e-09, 3.2e-09, 2.5e-09, 2e-09, 1.6e-09, 1.3e-09, 1e-09, 8e-10,
                                   6e-10, 5e-10)
)

### Functions
# Check that a string is IUPAC notation
checkIUPAC = function(x){all(strsplit(x, "")[[1]] %in% c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"))}

# convert nucleotide to a factor
nucleotide_factor = function(x){factor(x, levels = c("A", "C", "G", "T"))}

# Function to analyze an individual .fasta file
# Converts to neccesary form to be used in later operations compatabile with a sanger derived control
make_ctrl_fasta_df = function(fasta_file){
  init_ctrl_seq = read_lines(fasta_file)[2]
  ctrl_df = data.frame(max_base = init_ctrl_seq %>% strsplit(., split = "") %>% unlist(),
                       base_call = init_ctrl_seq %>% strsplit(., split = "") %>% unlist()) %>%
    mutate(index = 1:NROW(max_base)) %>%
    as_tibble()
  return(ctrl_df)
}

# Function to analyze an individual .ab1 file
# Need to make sure the file ends in .ab1 for this function
make_ctrl_sanger_df = function(sanger_file){
  options(warn=-1)
  base_calls = makeBaseCalls(sanger_file)
  sanger_df = base_calls %>% peakAmpMatrix %>% data.frame() %>% as_tibble()
  colnames(sanger_df) = c("A_area","C_area","G_area","T_area")
  sanger_df %<>% 
    mutate(., max_base = {apply(., 1, which.max) %>% bases[.]}) %<>%
    mutate(Tot.Area = A_area + C_area + G_area + T_area,
           A_perc = 100*A_area / Tot.Area,
           C_perc = 100*C_area / Tot.Area,
           G_perc = 100*G_area / Tot.Area,
           T_perc = 100*T_area / Tot.Area) %<>%
    mutate(base_call = strsplit(x = toString(base_calls@primarySeq), split = "") %>% unlist) %<>%
    mutate(index = 1:NROW(.)) %<>%
    mutate(max_base_height = {ifelse(max_base == "A", A_area,
                                     ifelse(max_base == "C", C_area,
                                            ifelse(max_base == "G", G_area,
                                                   ifelse(max_base == "T", T_area,NA))))}) %<>%
    mutate(max_base_perc = {ifelse(max_base == "A", A_perc,
                                   ifelse(max_base == "C", C_perc,
                                          ifelse(max_base == "G", G_perc,
                                                 ifelse(max_base == "T", T_perc,NA))))}) %<>%
    mutate(pre_5 = lag(base_call, n = 5),
           pre_4 = lag(base_call, n = 4),
           pre_3 = lag(base_call, n = 3),
           pre_2 = lag(base_call, n = 2),
           pre_1 = lag(base_call, n = 1),
           pre_pentanucleotide = paste0(pre_5, pre_4, pre_3, pre_2, pre_1)) %<>%
    mutate(post_5 = lead(base_call, n = 5),
           post_4 = lead(base_call, n = 4),
           post_3 = lead(base_call, n = 3),
           post_2 = lead(base_call, n = 2),
           post_1 = lead(base_call, n = 1),
           post_pentanucleotide = paste0(post_1, post_2, post_3, post_4, post_5)) %<>%
    mutate(pre_height_5 = lag(max_base_height, n = 5),
           pre_height_4 = lag(max_base_height, n = 4),
           pre_height_3 = lag(max_base_height, n = 3),
           pre_height_2 = lag(max_base_height, n = 2),
           pre_height_1 = lag(max_base_height, n = 1)) %<>%
    mutate(post_height_5 = lead(max_base_height, n = 5),
           post_height_4 = lead(max_base_height, n = 4),
           post_height_3 = lead(max_base_height, n = 3),
           post_height_2 = lead(max_base_height, n = 2),
           post_height_1 = lead(max_base_height, n = 1))
}


### Makes a df for the edited sample that accounts for the misalignment of peaks
### Employs the secondary basecalls
### Not aligned for peaks that are shifted by a complete basecall
make_samp_sanger_df = function(samp_sanger, ctrl_seq){
  
  ### Align the phase of the primary and secondary basecalls in the sample sequence to that of the control
  ### This changes where and what the bases are called as
  ### Returns an object of class sangerseq
  phased_samp_sanger = samp_sanger %>%
    makeBaseCalls(.) %>%
    setAllelePhase(., ctrl_seq) # Does not require an actual chromatogram, so it is ammenable to just using a fasta file
  
  ### Return a data frame from the phased sample sangerseq object with the position of each 
  ### This method appears to return higher intensities for the noise, thus we're going to trust it more for detecting noise for modelling
  samp_peakAmpDF = phased_samp_sanger@peakPosMatrix %>%
    as.data.frame() %>%
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    na_if(., 0) %>%
    dplyr::mutate(primary_base_call = primarySeq(phased_samp_sanger) %>% as.character() %>% str_split(., pattern = "") %>% unlist(),
                  secondary_base_call = secondarySeq(phased_samp_sanger) %>% as.character() %>% str_split(., pattern = "") %>% unlist()) %>%
    dplyr::mutate(identical = {primary_base_call == secondary_base_call}) %>%
    dplyr::mutate(row_sum = as.integer(!is.na(A)) +
                    as.integer(!is.na(C)) +
                    as.integer(!is.na(G)) +
                    as.integer(!is.na(`T`))
    ) %>% 
    dplyr::mutate(peak_pos = {ifelse(grepl("A|C|G|T", primary_base_call),
                                     ifelse(primary_base_call == "A", A, 
                                            ifelse(primary_base_call == "C", C, 
                                                   ifelse(primary_base_call == "G", G, 
                                                          ifelse(primary_base_call == "T", `T`, "Error, stop!")))),
                                     pmin(A, C, G, `T`, na.rm = TRUE))}) %>%
    dplyr::mutate(A = ifelse(is.na(A), peak_pos, A) %>% as.numeric(),
                  C = ifelse(is.na(C), peak_pos, C) %>% as.numeric(),
                  G = ifelse(is.na(G), peak_pos, G) %>% as.numeric(),
                  `T` = ifelse(is.na(`T`), peak_pos, `T`) %>% as.numeric())
  
  ### Reformat df to be identical to output of make_ctrl_sanger_df()
  samp_peakAmpDF %<>%
    dplyr::mutate(A_area = samp_sanger@traceMatrix[samp_peakAmpDF$A, 1],
                  C_area = samp_sanger@traceMatrix[samp_peakAmpDF$C, 2],
                  G_area = samp_sanger@traceMatrix[samp_peakAmpDF$G, 3],
                  T_area = samp_sanger@traceMatrix[samp_peakAmpDF$`T`, 4]) %<>%
    dplyr::select(A_area:T_area, primary_base_call, secondary_base_call) %<>%
    mutate(., max_base = {apply(., 1, which.max) %>% bases[.]}) %<>%
    mutate(max_base = factor(max_base, levels = ACGT)) %<>%
    mutate(Tot.Area = A_area + C_area + G_area + T_area,
           A_perc = 100*A_area / Tot.Area,
           C_perc = 100*C_area / Tot.Area,
           G_perc = 100*G_area / Tot.Area,
           T_perc = 100*T_area / Tot.Area) %>%
    mutate(index = 1:NROW(Tot.Area)) %<>%
    mutate(max_base_height = {ifelse(max_base == "A", A_area,
                                     ifelse(max_base == "C", C_area,
                                            ifelse(max_base == "G", G_area,
                                                   ifelse(max_base == "T", T_area,NA))))}) %<>%
    mutate(max_base_perc = {ifelse(max_base == "A", A_perc,
                                   ifelse(max_base == "C", C_perc,
                                          ifelse(max_base == "G", G_perc,
                                                 ifelse(max_base == "T", T_perc,NA))))}) %<>%
    mutate(pre_5 = lag(max_base, n = 5),
           pre_4 = lag(max_base, n = 4),
           pre_3 = lag(max_base, n = 3),
           pre_2 = lag(max_base, n = 2),
           pre_1 = lag(max_base, n = 1),
           post_5 = lead(max_base, n = 5),
           post_4 = lead(max_base, n = 4),
           post_3 = lead(max_base, n = 3),
           post_2 = lead(max_base, n = 2),
           post_1 = lead(max_base, n = 1),
           pre_height_5 = lag(max_base_height, n = 5),
           pre_height_4 = lag(max_base_height, n = 4),
           pre_height_3 = lag(max_base_height, n = 3),
           pre_height_2 = lag(max_base_height, n = 2),
           pre_height_1 = lag(max_base_height, n = 1),
           post_height_5 = lead(max_base_height, n = 5),
           post_height_4 = lead(max_base_height, n = 4),
           post_height_3 = lead(max_base_height, n = 3),
           post_height_2 = lead(max_base_height, n = 2),
           post_height_1 = lead(max_base_height, n = 1)
    )
  
  return(samp_peakAmpDF)
}

# The alignment index and subsequent filtering is currently off -- need to adapt to take the longest region of alignment, and then only filter our those sequences.
align_sanger_dfs = function(control_df, sample_df){
  
  ctrl_seq = control_df$max_base %>% paste0(., collapse = "") %>% DNAString()
  sample_seq = sample_df$max_base %>% paste0(., collapse = "") %>% DNAString()
  
  if(ctrl_seq@length > sample_seq@length){
    alignment = pairwiseAlignment(pattern = sample_seq, subject = ctrl_seq)
    control_df = control_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@subject@range@start) %>%
      filter(align_index < (alignment@subject@range@start+alignment@subject@range@width))
    sample_df = sample_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@pattern@range@start) %>%
      filter(align_index < alignment@pattern@range@start+alignment@pattern@range@width)
  } else
  {
    alignment = pairwiseAlignment(pattern = ctrl_seq, subject = sample_seq)
    control_df = control_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@pattern@range@start) %>%
      filter(align_index < alignment@pattern@range@start+alignment@pattern@range@width)
    sample_df = sample_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@subject@range@start) %>%
      filter(align_index < (alignment@subject@range@start+alignment@subject@range@width))
  }
  return(list(control_df, sample_df))
}

### Convert the abif to fastq
abif_to_fastq = function (seqname = "sample", path, trim = TRUE, cutoff = 1, 
                          min_seq_len = 20, offset = 33, recall = FALSE) 
{
  sangerseqr <- requireNamespace("sangerseqR")
  stopifnot(isTRUE(sangerseqr))
  abif <- sangerseqR::read.abif(path)
  if (is.null(abif@data$PCON.2)) {
    message(sprintf("failed on %s", seqname))
    return()
  }
  nucseq <- substring(abif@data$PBAS.2, 1, length(abif@data$PLOC.2))
  if (!typeof(abif@data$PCON.2) == "integer") {
    num_quals <- utf8ToInt(abif@data$PCON.2)[1:length(abif@data$PLOC.2)]
  }
  else {
    num_quals <- abif@data$PCON.2[1:length(abif@data$PLOC.2)]
  }
  if (isTRUE(recall)) {
    recalled <- sangerseqR::makeBaseCalls(sangerseqR::sangerseq(abif))
    nucseq <- sangerseqR::primarySeq(recalled, string = TRUE)
    if (nchar(nucseq) != length(num_quals)) {
      trim <- FALSE
      num_quals <- rep(60, nchar(nucseq))
      warning("Length of quality scores does not equal length of\n              re-called base sequence, ignoring quality scores")
    }
  }
  if (trim == FALSE) {
    tmp1 = list(seqname = seqname, seq = nucseq, 
                quals = rawToChar(as.raw(num_quals + offset)))
    return(tmp1)
  }
  trim_msg <- "Sequence %s can not be trimmed because it is shorter than the trim\n               segment size"
  if (nchar(nucseq) <= min_seq_len) {
    warning(sprintf(trim_msg, seqname))
    return()
  }
  scores = cutoff - 10^(num_quals/-10)
  running_sum <- rep(0, length(scores) + 1)
  for (i in 1:length(scores)) {
    num <- scores[i] + running_sum[i]
    running_sum[i + 1] <- ifelse(num < 0, 0, num)
  }
  trim_start <- min(which(running_sum > 0)) - 1
  trim_finish <- which.max(running_sum) - 2
  if (trim_finish - trim_start < min_seq_len - 1) {
    warning(sprintf(trim_msg, seqname))
    return()
  }
  tmp2 = list(seqname = seqname, seq = substring(nucseq, 
                                                 trim_start, trim_finish), quals = rawToChar(as.raw(num_quals[trim_start:trim_finish] + 
                                                                                                      offset)))
  return(tmp2)
}

### Takes sequence and runs through a while loop until they are properly trimmed
### Appears for CTSS 11 sample
align_and_trim = function(pattern_seq, subject_seq, min_continuity = 15){
  
  gap_length = min_continuity - 1
  raw_pattern = lapply(FUN = rep, X = 1:gap_length, x = "[A-Z]") %>%
    lapply(FUN = paste0, X = ., collapse = "") %>%
    unlist()
  gsub_pattern = raw_pattern %>%
    paste0("^", ., "-|") %>%
    paste0(collapse = "") %>%
    paste0(., raw_pattern %>%
             paste0("-", ., "$|") %>%
             paste0(collapse = ""),
           collapse = "") %>%
    gsub('.{1}$','', .)
  
  input_pattern = pattern_seq
  input_subject = subject_seq
  
  alignment = pairwiseAlignment(pattern = pattern_seq, subject = subject_seq)
  
  output_pattern = alignment@pattern %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
  output_subject = alignment@subject %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
  
  while(input_pattern != output_pattern | input_subject != output_subject){
    alignment = pairwiseAlignment(pattern = output_pattern, subject = output_subject)
    
    old_output_pattern = output_pattern
    old_output_subject = output_subject
    
    new_output_pattern = alignment@pattern %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
    new_output_subject = alignment@subject %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
    
    output_pattern = new_output_pattern
    output_subject = new_output_subject
    
    input_pattern = old_output_pattern
    input_subject = old_output_subject
  }
  
  return(list("pattern" = alignment@pattern %>% as.character() %>% gsub("-", "", .),
              "subject" = alignment@subject %>% as.character() %>% gsub("-", "", .),
              "alignment" = alignment))
}

### Substitute multiple string positions simultaneously
subchar <- function(string, pos, char) { 
  for(i in pos) { 
    substr(string, i, i) = char
  } 
  string 
} 

### Make a dataframe with all of the ZAGA information
### Backup ZAGA vector
replacement_zaga = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
                     0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
                     0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)

make_ZAGA_df = function(sanger_df, p_adjust){
  nvals <- list()
  nvals$A = filter(sanger_df, max_base != "A")$A_area 
  nvals$C = filter(sanger_df, max_base != "C")$C_area 
  nvals$G = filter(sanger_df, max_base != "G")$G_area 
  nvals$T = filter(sanger_df, max_base != "T")$T_area 
  
  n_models =   n_models <-lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
    {x = replacement_zaga; message("Replacement vector used for low noise.")} # add noise if all 0s, or all 0s and one other value.
    tryCatch(gamlss((x)~1, family = ZAGA), error=function(e) # Progressively step up the mu.start if it fails
      tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 1), error=function(e) 
        tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 2), error=function(e) 
          tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 3), error=function(e) # additional step added.
            gamlss((x)~1, family = ZAGA, mu.start = mean(x))
          )
        )
      )
    )
    # throws errors when a completely 0 vector
  })
  
  null_m_params = lapply(n_models, FUN = function(x){
    mu <- exp(x$mu.coefficients[[1]])
    sigma <- exp(x$sigma.coefficients[[1]])
    nu.logit <- x$nu.coefficients[[1]]
    nu <- exp(nu.logit)/(1+exp(nu.logit))
    fillibens <-cor(as.data.frame(qqnorm(x$residuals, plot = FALSE)))[1,2]
    crit = qZAGA(p = 1-p_adjust, mu = mu, nu = nu, sigma = sigma)
    
    return(data.frame(μ = mu, σ = sigma, n = nu, crit = crit, `Fillibens` = fillibens))
  })
  
  null_m_params %>%
    plyr::ldply(., "data.frame") %>%
    dplyr::rename(Base = `.id`) %>%
    return()
}

### convert negatives to zeros
neg_to_zero = function(x){ifelse(x < 0, 0, x)}

### Adjust the height and percent values based on GBM and the significance
gbm_adjust = function(sanger_df, wt, boi, motif, sample_file, critical_values){
  sanger_df %>%
    # only bases of interest
    dplyr::filter(grepl(wt, ctrl_max_base)) %>% # To only pull out the bases of interest from the motif
    dplyr::select(index, ctrl_index, motif_id, ctrl_max_base, ctrl_max_base_perc, max_base, A_area:T_area, A_eta:T_eta) %>%
    # If the channel is not found in the potential edits, then it is not applicable or NA, otherwise
    # If the height of the channel is greater than the critical value for that channel,
    # return TRUE, otherwise, return FALSE
    mutate(A_sig = if(!grepl("A", boi)) {FALSE} else {ifelse(A_area >= critical_values[1], TRUE, FALSE)},
           C_sig = if(!grepl("C", boi)) {FALSE} else {ifelse(C_area >= critical_values[2], TRUE, FALSE)},
           G_sig = if(!grepl("G", boi)) {FALSE} else {ifelse(G_area >= critical_values[3], TRUE, FALSE)},
           T_sig = if(!grepl("T", boi)) {FALSE} else {ifelse(T_area >= critical_values[4], TRUE, FALSE)}) %>%
    mutate(A_eta = ifelse(A_sig, A_eta, 0),
           C_eta = ifelse(C_sig, C_eta, 0), 
           G_eta = ifelse(G_sig, T_eta, 0), 
           T_eta = ifelse(T_sig, T_eta, 0)) %>%
    mutate(N_sig = A_sig + C_sig + G_sig + T_sig) %>%
    #mutate(Tot_eta = A_eta + C_eta + G_eta + T_eta) %>%
    #mutate(A_eta = ifelse(A_eta != 0, A_eta, NA),
    #       C_eta = ifelse(C_eta != 0, C_eta, NA), 
    #       G_eta = ifelse(G_eta != 0, G_eta, NA), 
    #       T_eta = ifelse(T_eta != 0, T_eta, NA)) %>%
    #mutate(A_adj = A_area*ifelse(A_sig,( 1 + (1 - ( (N_sig*A_eta) / Tot_eta ) ) ), 1),
    #       C_adj = C_area*ifelse(C_sig,( 1 + (1 - ( (N_sig*C_eta) / Tot_eta ) ) ), 1),
    #       G_adj = G_area*ifelse(G_sig,( 1 + (1 - ( (N_sig*G_eta) / Tot_eta ) ) ), 1),
    #       T_adj = T_area*ifelse(T_sig,( 1 + (1 - ( (N_sig*T_eta) / Tot_eta ) ) ), 1)) %>%
    #mutate(Tot_adj = A_adj + C_adj + G_adj + T_adj,
    #       A_adj_perc =  A_adj / Tot_adj, 
  #       C_adj_perc =  C_adj / Tot_adj,
  #       G_adj_perc =  G_adj / Tot_adj,
  #       T_adj_perc =  T_adj / Tot_adj) %>%
  mutate(Tot_area = A_area + C_area + G_area + T_area,
         A_perc =  A_area / Tot_area, 
         C_perc =  C_area / Tot_area,
         G_perc =  G_area / Tot_area,
         T_perc =  T_area / Tot_area) %>%
    mutate(wt_base = wt) %>%
    mutate(sig = {ifelse(wt_base == "A", (C_sig + G_sig + T_sig),
                         ifelse(wt_base == "C", (A_sig + G_sig + T_sig),
                                ifelse(wt_base == "G", (A_sig + C_sig + T_sig),
                                       ifelse(wt_base == "T", (A_sig + C_sig + G_sig), 0))))}) %>%
    mutate(sig = {ifelse(sig > 0, TRUE, FALSE)}) %>%
    #mutate(sig = {ifelse(wt == "A", 1,
    #                     ifelse(wt == "C", 2,
    #                             ifelse(wt == "G", 3,
    #                                   ifelse(wt == "T", 4, 0))))}) %>%
    dplyr::select(index, ctrl_index, motif_id, ctrl_max_base, ctrl_max_base_perc, max_base, A_area:T_area, A_perc:T_perc,  A_sig:T_sig, sig) %>%
    mutate(motif = motif, sample_file = sample_file)
}


### Functions for plotting

# Function for plotting raw sample data
plotRawSample = function(raw_sample_df, sample_alt, pre_cross_align_sample_df){
  raw_sample_df %>%
    ggplot(aes(x = index, y = max_base_perc)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0)) +
    scale_x_continuous(limits = c(0, max(raw_sample_df$index)), expand = c(0,0)) +
    geom_bar(data = sample_alt, aes(x = index, y = 100), stat = "identity", fill = "#53BCC2", color = "#53BCC2") +
    geom_rect(xmin = min(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(pre_cross_align_sample_df$index), ymax = 100,
              fill = "white", color = "black", alpha = 0) +
    xlab("Position in sample file") +
    ylab("Percent signal of basecall") +
    geom_hline(yintercept = mean(pre_cross_align_sample_df$max_base_perc), color = "darkred", size = 1) +
    geom_line() +
    geom_rect(xmin = 0, ymin = 0, xmax = min(pre_cross_align_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    geom_rect(xmin = max(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(raw_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    geom_rect(xmin = min(pre_cross_align_sample_df$index) +
                (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025 - 5,
              ymin = 3,
              xmax = min(pre_cross_align_sample_df$index) +
                (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025 + 215,
              ymax = 20,
              fill = "white", color = "black") +
    annotate(geom = "text",
             label = paste0("Average percent signal (", round(mean(pre_cross_align_sample_df$max_base_perc), 1),"%)") ,
             color = "darkred",
             x = min(pre_cross_align_sample_df$index) +
               (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025,
             y = 17,
             hjust = 0,
             size = 6) +
    annotate(geom = "text", label = "Low phred trimmed regions", color = "grey30",
             x = min(pre_cross_align_sample_df$index) +
               (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025,
             y = 12,
             hjust = 0,
             size = 6) +
    annotate(geom = "text", label = "Motif of interest", color = "#53BCC2",
             x = min(pre_cross_align_sample_df$index) +
               (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025,
             y = 7,
             hjust = 0,
             size = 6) +
    theme_classic(base_size = 18)
}

# Function for plotting trimmed sample data
plotTrimmedSample = function(sample_df, pre_cross_align_sample_df, output_sample_alt, raw_sample_df, sample_alt){
    sample_df %>%
      mutate(perc = (100-ctrl_max_base_perc) %>% round(., 0)) %>%
      ggplot(aes(x = index, y = perc)) +
      geom_line() +
      scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0)) +
      scale_x_continuous(limits = c(0, max(raw_sample_df$index)), expand = c(0,0)) +
      # geom_bar(data = sample_alt, aes(x = index, y = 100), stat = "identity", fill = "#53BCC2", color = "#53BCC2", alpha = 0.1) +
      # geom_bar(data = output_sample_alt, aes(x = index, y = 100), stat = "identity", fill = "lightblue", color = "lightblue") +
      geom_rect(xmin = min(pre_cross_align_sample_df$index), ymin = 0,
                xmax = max(pre_cross_align_sample_df$index), ymax = 100,
                fill = "white", color = "black", alpha = 0) +
      xlab("Position in sample file") +
      ylab("Percent noise in WT basecall") +
      geom_hline(yintercept = mean(pre_cross_align_sample_df$max_base_height), color = "darkred", size = 1) +
      geom_line() +
      geom_point(data = output_sample_alt %>%
                   mutate(sig = {ifelse(sig, "Signficant", "Non-significant")}) %>%
                   mutate(perc = (100-ctrl_max_base_perc) %>% round(., 0)),
                 aes(fill = motif_id, x = index, y = perc, alpha = sig), pch = 21, size = 2, color = "black") +
      geom_rect(xmin = 0, ymin = 0, xmax = min(pre_cross_align_sample_df$index), ymax = 100,
                fill = "grey", color = "black", alpha = 0.01) +
      geom_rect(xmin = max(pre_cross_align_sample_df$index), ymin = 0,
                xmax = max(raw_sample_df$index), ymax = 100,
                fill = "grey", color = "black", alpha = 0.01) +
      theme_classic(base_size = 18) +
      scale_alpha_manual(values = c(0.3, 1)) +
      theme(legend.position = "none")
}

# Function for producing the table data on editing
calculateEditingData = function(output_sample_alt, edit) {
  output_sample_alt %>%
    dplyr::select(index, ctrl_index, A_perc:T_perc) %>%
    gather(base, perc, A_perc:`T_perc`) %>%
    inner_join(., output_sample_alt %>%
                 dplyr::select(index, A_sig:T_sig) %>%
                 gather(base_sig, sig, A_sig:T_sig)
    ) %>%
    mutate(base = gsub("_perc", "", base), base_sig = gsub("_sig", "", base_sig)) %>%
    mutate(perc = perc*100, sig = {ifelse(sig, "Significant", "Non-significant")}) %>%
    filter(base == base_sig) %>%
    filter(grepl(pattern = edit, x = base)) %>%
    mutate(tmp = 1) %>%
    group_by(base) %>%
    mutate(tally = cumsum(tmp)) %>%
    dplyr::select(-tmp, -base_sig)
}

plotEditingData = function(editing_data) {
  editing_data %>%
    ggplot(aes(x = sig, y = perc, color = sig)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
    geom_bar(aes(fill = sig), stat = "summary", fun.y = "mean",
             color = "black", alpha = 0.3, show.legend = F) +
    geom_jitter(size = 2, alpha = 0.7) +
    ylab("Percent height") +
    xlab("") +
    labs(color = "Potential edit") +
    guides(fill = NULL) +
    theme_classic(base_size = 18) +
    theme(aspect.ratio = 4/1,
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +#sum(grepl(edit, bases))) +
    facet_wrap(.~base, nrow = 1, scales = "free_x")
}

tableEditingData = function(editing_data) {
  editing_data %>%
      group_by(base, sig) %>%
      dplyr::summarize(Mean = mean(perc), Max = max(perc), Min = min(perc), SD = sd(perc), N = length(perc)) %>%
      dplyr::rename(Base = base, Significance = sig)
}

# Define geom_chromatogram function
# input the sanger object, start index, end index and output chromatogram with underlying base percents
geom_chromatogram = function(sanger, start_index, end_index){
  
  # define bases
  bases = c("A", "C", "G", "T")
  
  # define the basecalls  in the sanger object initially to set the proper phasing to the data
  sanger = makeBaseCalls(sanger)
  
  # make rawtraces from the sanger trace matrix
  rawtraces = sanger@traceMatrix %>%
    # the trace matrix is the scatter of points that make up the chromatogram
    # convert the large matrix of traces into a tibble object
    as_tibble() %>%
    # name the columns by base channel
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    # create a column 1 through n trace measurements taken
    rowid_to_column("x")
  
  # create the peakPosMatrix giving the x coordinate of each basecall
  peakPosMatrix = sanger@peakPosMatrix %>% as_tibble()
  # create column names for each channel
  colnames(peakPosMatrix) = c("A", "C", "G", "T")
  # format peakPosMatrix
  peakPosMatrix = peakPosMatrix %>%
    # add the basecall for each position
    mutate(basecall = sanger@primarySeq %>% as.character %>% strsplit(split = "") %>% unlist()) %>%
    rowid_to_column("index") %>%
    gather(base, x, A:`T`) %>%
    filter(basecall == base) %>%
    arrange(index) %>%
    dplyr::select(index, basecall, x)
  
  # Make the traces data frame
  traces = peakPosMatrix %>% 
    # filter only the trace information at positions that are registered as the position of the peak height
    inner_join(., rawtraces) %>%
    # define the previous base index
    mutate(x_1 = lag(x, default = 1)) %>%
    # start the chromatogram for each base as the position between the two bases divided by two
    # start the chromatogram as the x coordinate bewteen the first  base of interest and the n-1 base
    mutate(start = floor((x + x_1) / 2)) %>%
    # end the chromatogram for each base as the position
    # empirically derived I believe
    mutate(end = lead(start, default = (max(x) + 6))) %>%
    # define the total area for each basecall
    mutate(Tot.Area = A + C + G + `T`) %>%
    # calculate percent heights for each base from the total area
    mutate(A_perc = round(A*100/Tot.Area),
           C_perc = round(C*100/Tot.Area),
           G_perc = round(G*100/Tot.Area),
           T_perc = round(`T`*100/Tot.Area))
  
  # The position list is for every basecall, these are the corresponding x-coordinates
  # This is a list object, where the primary index of the list represents the basecall index
  # The vector element under each index represents the x-coordinates of the peakAmpMatrix for each basecall
  position_list = mapply(FUN = seq,
                         from = c(1, traces$end),
                         to = c(traces$end, max(traces$end)+10),
                         by = 1) %>%
    head(., -1)
  
  # make the names of the position_list indentical to the index of the position list
  names(position_list) = c(1:length(position_list))
  
  indices = ldply(position_list, "data.frame") %>%
    dplyr::rename(index = ".id", x = "data.frame") %>%
    as_tibble() %>%
    mutate(index = as.numeric(index))
  
  # This joins the basecall index to the rawtrace x coordinate based data
  rawtraces %<>% inner_join(., indices)
  
  # Enter the plotting based on the indices
  # defines the basecall_indices as the start and end indices of interest for making the chromatogram
  basecall_indices = c(start_index, end_index)
  
  # determine the start and end x coordinates of the rawtraces to filter on for plotting
  index_filters = rawtraces %>%
    # filter the rawtraces data to just include the indices that are within the region of interest
    # some of the bleedthrough is still observed at this point
    filter(index >= basecall_indices[1] & index <= basecall_indices[2]) %>%
    # the x coordinates corresponding to the region of interest
    .$x %>%
    # takes the min and max of the peakAmpMatrix x cooridinates
    quantile(., c(0, 1)) %>%
    # adds an x coordinate cushion
    {. + c(-6, 6)}
  
  # gather the rawtraces data for ggplot format
  # filter in only traces that fall within the x coordinates rawtraces index filter
  # could merge this line of code with the preceeding index_filters chunk
  plotTraces = rawtraces %>%
    gather(base, height, A:`T`) %>%
    filter(x >= index_filters[1] & x <= index_filters[2]) %>%
    # This join operation is used to creat a base code that later allows a manual gradient to be employed
    inner_join(.,
               data.frame(stringsAsFactors = FALSE,
                          base = c("A", "C", "G", "T"),
                          base_code = c("1000", "2000", "3000", "4000"))
    ) %>%
    # calculate the min and max height for geom_ribbon
    mutate(max_height = pmax(height), min_height = 0)
  
  # calculate the y coordinates for the chromatogram
  y_bases = max(plotTraces$height)*(1/2) %>%
  {.*-c(2/11,4/11,6/11,8/11,10/11)}
  
  # calculate the number of bases involved
  n_bases = basecall_indices[2]-basecall_indices[1]+1
  
  # calculate the x coordinates for the chromatogram
  x_bases = (index_filters[2]-index_filters[1]) %>%
  {.*(seq(2, 2*n_bases+1, 2)/(2*n_bases + 2))} %>%
  {.+index_filters[1]}
  
  tile_plot = traces %>%
    filter(index >= basecall_indices[1] & index <= basecall_indices[2]) %>%
    dplyr::select(basecall, A_perc, C_perc, G_perc, T_perc) %>%
    gather(col, labels, basecall:T_perc) %>%
    mutate(y = rep(y_bases, each = n_bases),
           x = rep(x_bases, times = 5))
  
  base_annotation = data.frame(label = bases,
                               x = rep(min(plotTraces$x), 4),
                               y = tail(rev(sort(unique(tile_plot$y))), -1))
  
  ## Create functions for the fill color palettes
  colors1 = colorRampPalette(colors = c("white", "white"))
  colors2 = colorRampPalette(colors = c("white", "#e41a1c"))
  colors3 = colorRampPalette(colors = c("#e41a1c", "#e41a1c"))
  colors4 = colorRampPalette(colors = c("#e41a1c", "#377eb8"))
  
  # establish a fill_key vector
  fill_key = c(0:6, 7:19, 20:49, 50:100, 1000, 2000, 3000, 4000)
  
  # establish the colors in the fill_key
  fill_colors = c(colors1(length(0:6)),
                  colors2(length(7:19)),
                  colors3(length(20:49)),
                  colors4(length(50:100)),
                  "#32CD32", "#4169E1", "#121212", "#FC4338")
  
  # tie together the fill_key and fill_colors using the names() function
  names(fill_colors) = fill_key
  
  chromatogram = plotTraces %>%
    mutate(max_height = pmax(height), min_height = 0) %>%
    ggplot(data = ., aes(x = x, y = height)) +
    geom_ribbon(aes(ymin = min_height, ymax = max_height, color = base, fill = base_code), alpha = 0.1) +
    scale_color_manual(values = c("A" = "#4daf4a", "C" = "#377eb8", "G" = "black", "T" = "#e41a1c")) +
    geom_tile(data = filter(tile_plot, col != "basecall"),
              aes(y = y, x = x, fill = as.character(as.numeric(labels))), color = "black") +
    scale_fill_manual(values = fill_colors) +
    geom_text(data = tile_plot, aes(x = x, y = y, label = labels), color = "black") +
    geom_text(data = base_annotation, aes(x = x, y = y, label = label), color = "black") +
    theme_void(base_size = 36) +
    theme(legend.position = "none",
          aspect.ratio = 1/1.5)
  
  return(chromatogram)
}

geom_coxplot = function(editingData){
  sigEditingData = editingData %>% filter(sig == "Significant")
  proportion = 2/3
  n_points = max(editingData$tally)
  editingData %>%
    ggplot(aes(x = reorder(as.character(tally), tally), y = 1, fill = perc, color = reorder(sig, sig))) +
    geom_tile(width = 0.9, height = 0.9, lwd  = 1) +
    geom_point(data = sigEditingData, aes(x = n_points + (perc/n_points)*proportion, y = 1)) +
    scale_x_discrete(limits = c(0, n_points/proportion)) +
    #geom_tile(aes(data = sigEditingData, x = reorder(as.character(tally), tally), y = 1, fill = perc, color = sig)) +
    # scale_x_discrete(breaks = unique(editingData$tally), labels = unique(editingData$tally)) +
    ylab("") +
    xlab("") +
    facet_wrap(.~base, nrow = length(unique(editingData$base))) +
    scale_fill_gradient(low = "white", high = "#650C15") +
    scale_color_manual(values = c("Non-significant" = "grey", "Significant" = "black")) +
    theme_minimal() +
    theme(aspect.ratio = 1/length(unique(editingData$ctrl_index)),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
}

# Define logError function
# This function is used to log application errors for troubleshooting
logError = function(expr, step){
  require(berryFunctions)
  if(!berryFunctions::is.error(expr))
  {paste0(step, " is successful.")} else 
  {paste0("*Error: ", step, " failed.*")}
}

# Define which.closest
# A function to calcualte which value in a string is the closest
which.closest = function (vec, x) 
{
  which.min(abs(vec - x))
}
