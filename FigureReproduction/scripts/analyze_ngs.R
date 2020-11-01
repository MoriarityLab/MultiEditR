### Analyze EditR NGS comparison

###########################################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

##### Dependencies
### All libraries are loaded from the run_for_ngs.R and the dependencies.R
source("/Users/kluesner/Desktop/Research/EditR/multiEditR/program/working_branch/run_for_ngs.R")

##### Functions
# Calculates the R-squared value of a relationship
rsq = function(x, y){
  y_bar = mean(y)
  y_hat = x
  ss_total = sum((y - y_bar)^2)
  ss_residual = sum((y - y_hat)^2)
  1-(ss_residual/ss_total)
}

# Converts a character string to its complement 
complement = function(x){as.character(reverseComplement(reverse(DNAString(x))))}

# Converts a character string to its reverse complement
revcom = function(x){as.character(reverseComplement(DNAString(x)))}

# Returns the GC content of a fasta file
fastaGC = function(fasta){
  sequence = unlist(strsplit(read_lines(fasta)[2], split = ""))
  sum(sequence == "C" | sequence == "G")/length(sequence)
}

# Determine which is the maxbase for a data frame
maxBase = function(data, columns = c("A", "C", "G", "T")){
  mapply(FUN = function(i){names(which.max(data[i,columns]))},
         i = 1:NROW(data))
}

viewSamp = function(data){data %>% dplyr::select(ctrl_index, A_perc:T_perc, ctrl_base_call) %>% View}


##### Parameters
# Base colors
base_colors = c("A" = "forestgreen", "C" = "royalblue", "G" = "black", "T" = "firebrick")

# Issue with 6
# Issue with Rev samples
# CTSB looking at wrong strand?
# WT1 Rab7 is a great example

# analyzeNGS = function(
# parameters = readxl::read_excel(paste0(directory, "/parameters.xls")) %>%
#   mutate(sample_file = {ifelse(fwd,
#                                paste0(directory, "/fwd/", sample_file),
#                                paste0(directory, "/rev/", sample_file))}) %>%
#   mutate(ctrl_file = paste0(directory, "/fasta/", ctrl_file)) %>%
#   mutate(ngs_file = paste0(directory, "/xlsx/", ngs_file))
# 
# i = 1
#   directory = "/Users/kluesner/Desktop/Research/EditR/multiEditR/NGS/final"
#   sample_file = parameters$sample_file[i]
#   ctrl_file = parameters$ctrl_file[i]
#   ngs_file = parameters$ngs_file[i]
#   fwd = parameters$fwd[i]
#   motif =parameters$motif[i]
#   wt = parameters$wt[i]
#   edit = parameters$edit[i]
#   phred_cutoff = parameters$phred_cutoff[i]
#   p_value = 0.01
#   include = parameters$include[i]
# # ){

analyzeNGS = function(
  directory,
  sample_file,# = parameters$sample_file[i],
  ctrl_file,# = parameters$ctrl_file[i],
  ngs_file,# = parameters$ngs_file[i],
  fwd,# = parameters$fwd[i],
  motif,# =parameters$motif[i],
  wt,# = parameters$wt[i],
  edit,# = parameters$edit[i],
  phred_cutoff,#= parameters$phred_cutoff[i],
  p_value = 0.01,#
  include # = parameters$include[i]
){
# Handling to reverse complement the base information based on the read orientation of the .ab1 file
init_orient = fwd
if(fwd) {} else {motif = revcom(motif); wt = revcom(wt); edit = revcom(edit)}

##### Generate EditR data
### Generates a list with alternative df, a null df, and the control sequence for reference alignment
sample_data = extract_df(sample_file = sample_file,
                         ctrl_file = ctrl_file,
                         p_value = p_value,
                         motif = motif,
                         wt = wt,
                         edit = edit,
                         phred_cutoff = phred_cutoff,
                         fwd_orientation = fwd)

### Assign parts of data to individual objects to make it easier to work with
samp_df = sample_data[[1]] #%>% dplyr::rename(ctrl_index = ctrl_post_aligned_index)
alt_df = sample_data[[2]] # This should have all the data needed for the comparison
ctrl_seq = sample_data[[3]]
ctrl_df = sample_data[[4]]
zaga_parameters = sample_data[[5]]

# 1. Load in the data
### Create NGS sequence from ngs dataframe
# Read in the xlsx file Ricca gave to me
ngs_sample =readxl::read_excel(ngs_file) %>%
  # Replace "[" and "]" with "" to separate easier
  mutate(`BaseCount[A,C,G,T]` = gsub('[[]|[]]', "", `BaseCount[A,C,G,T]`)) %>%
  # separate on the ", " to make individual columns for each base
  separate(`BaseCount[A,C,G,T]`, into = c("A", "C", "G", "T"), sep = ", ", convert = TRUE) %>%
  # Replace "[" and "]" with "" to separate easier for KO data
  mutate(`gBaseCount[A,C,G,T]` = gsub('[[]|[]]', "", `gBaseCount[A,C,G,T]`)) %>%
  # separate on the ", " to make individual columns for each base for KO data
  separate(`gBaseCount[A,C,G,T]`, into = c("gA", "gC", "gG", "gT"), sep = ", ", convert = TRUE) %>%
  mutate(gA = {ifelse(gA == "-", NA, gA)}) %>%
  # Convert read columns to numeric %>%
  mutate_at(., c("A", "C", "G", "T", "gA", "gC", "gG", "gT"), as.integer) %>%
  # Rename the Coverage column as Reads
  dplyr::rename(Reads = `Coverage-q25`) %>%
  dplyr::mutate(samp = gsub(directory, "", ngs_file) %>%
                  gsub("/xlsx/|[.]xlsx", "", .)) %$%
  # Establish the KO data for fisher test
  inner_join(., 
             
             ## Establish the null probabilities for Fisher exact test
             
             # Only include rows that have read values
             filter(., !is.na(gA)) %>%
               
               # Group the data by the WT base in the reference genome
               dplyr::group_by(., Reference) %>% 
               
               # For each reference base establish the total number from each base channel
               dplyr::summarise(sum_reference_reads_A = sum(gA, na.rm = T),
                                sum_reference_reads_C = sum(gC, na.rm = T),
                                sum_reference_reads_G = sum(gG, na.rm = T),
                                sum_reference_reads_T = sum(gT, na.rm = T)) %>%
               
               # reshape the data
               gather(base, reads, sum_reference_reads_A:sum_reference_reads_T) %>%
               mutate(base = gsub("sum_reference_reads_", "", base)) %>%
               
               # pair the reference bases with their edited base identity
               inner_join(., tibble(Reference = c("A", "C", "G", "T"), Alternative = c("G", "T", "A", "C"))) %>%
               
               # only keep rows where the information is either for the Reference (Unedited) or Alternative (Edited) base
               filter(base == Reference | base == Alternative) %>%
               
               # reshape the data
               mutate(Type = {ifelse(base == Reference, "gRef", "gAlt")}) %>%
               dplyr::select(Reference, Alternative, reads, Type) %>%
               tidyr::spread(., value = "reads", key = "Type")
  ) %$%
  inner_join(.,
             mutate(., a = A, c = C, g = G, t = `T`) %>%
               gather(., A:`T`, key = "base", value = "reads") %>% 
               filter(base == Reference | base == Alternative) %>% 
               mutate(Type = {ifelse(base == Reference, "Ref", "Alt")}) %>%
               dplyr::select(., Position, reads, Type) %>%
               tidyr::spread(., value = "reads", key = "Type")
  )

  # # 02.02.2020
  # # Establish the editing indices for each sample
  # # group by the reference base call
  # dplyr::group_by(Reference) %>%
  # # dplyr::mutate(AEI_ngs = sum(G)/(sum(G) + sum(A)),test = sum(G))
  # # Calculate the editing index for each transition mutation
  # # Need to specify dplyr:: over plyr::, otherwise it doesn't behave with groupings.
  # dplyr::mutate(AEI_ngs = (sum(G) / (sum(G) + sum(A))),
  #        CEI_ngs = (sum(`T`) / (sum(C) + sum(`T`))),
  #        GEI_ngs = (sum(A) / (sum(A) + sum(G))),
  #        TEI_ngs = (sum(C) / (sum(`T`) + sum(C)))) %>%
  # # ungroup
  # ungroup() %>%
  # # Keep only the EI_index for each reference base
  # gather(EI_base, EI_ngs, AEI_ngs:TEI_ngs) %>%
  # mutate(EI_base = gsub("EI_ngs", "", EI_base)) %>%
  # filter(EI_base == Reference) %>%
  # dplyr::select(-EI_base)

# 2. Arrange the NGS in the correct order
# IF the NGS is on the positive strand, THEN arrange in ascending order, ELSE arrange in descending order
if(names(sort(table(ngs_sample$Strand), T)[1]) == "1")
  {ngs_sample = arrange(ngs_sample, Position) %>%
    mutate(ngs_index  = 1:NROW(.))} else
    {ngs_sample = arrange(ngs_sample, desc(Position)) %>%
                            mutate(ngs_index  = 1:NROW(.))}

# 3. Assign NGS reference sequence
ngs_ref_seq = ngs_sample$Reference %>% paste0(., collapse = "")

# If the alignment requires a reverse complement then reverse complement the data frame
alignment = matchPattern(ctrl_seq, ngs_ref_seq, max.mismatch =  143)

if(length(alignment@ranges@start))
  {ngs_sample = ngs_sample %>% mutate(alignment_index = ngs_index)} else
  {alignment = matchPattern(ctrl_seq,
                            reverseComplement(DNAString(ngs_ref_seq)),
                            max.mismatch =  143)
   ngs_sample = ngs_sample %>%
     arrange(desc(ngs_index)) %>%
     mutate(alignment_index = 1:NROW(.)) %>%
     mutate(A_tmp = A, C_tmp = C, G_tmp = G, T_tmp = `T`) %>%
     dplyr::select(-A, -C, -G, -`T`) %>%
     dplyr::rename(`T` = A_tmp, G = C_tmp, C = G_tmp, A = T_tmp) %>%
     mutate(max_base = maxBase(.)) %>%
     inner_join(., data.frame(Reference = c("A", "C", "G", "T"),
                              reference = c("T", "G", "C", "A"))) %>%
     dplyr::select(-Reference) %>%
     dplyr::rename(Reference = reference)
   
   
  }

# Filter and assign control index
# also add editing index
ngs_sample =  ngs_sample %>%
  filter(alignment_index >= alignment@ranges@start,
         alignment_index < (alignment@ranges@start + alignment@ranges@width)) %>%
  mutate(ctrl_index = 1:NROW(.)) %>%
  mutate(A_perc = A/Reads, C_perc = C/Reads, G_perc = G/Reads, T_perc = `T`/Reads) %>%
  # 02.02.2020
  # Establish the editing indices for each sample
  # group by the reference base call
  dplyr::group_by(Reference) %>%
  # dplyr::mutate(AEI_ngs = sum(G)/(sum(G) + sum(A)),test = sum(G))
  # Calculate the editing index for each transition mutation
  # Need to specify dplyr:: over plyr::, otherwise it doesn't behave with groupings.
  dplyr::mutate(AEI_ngs = (sum(G) / (sum(G) + sum(A))),
                CEI_ngs = (sum(`T`) / (sum(C) + sum(`T`))),
                GEI_ngs = (sum(A) / (sum(A) + sum(G))),
                TEI_ngs = (sum(C) / (sum(`T`) + sum(C)))) %>%
  # ungroup
  ungroup() %>%
  # Keep only the EI_index for each reference base
  gather(EI_base, EI_ngs, AEI_ngs:TEI_ngs) %>%
  mutate(EI_base = gsub("EI_ngs", "", EI_base)) %>%
  filter(EI_base == Reference) %>%
  dplyr::select(-EI_base)
  

sample_name = gsub("_", " ", unique(ngs_sample$samp))

# Its easy to judge people, but it doesn't pay off in the long run

# NGS data in region of interest
ngs_plot = ngs_sample %>%
  dplyr::select(-A, -C, -G, -`T`) %>%
  dplyr::rename(A = A_perc, C = C_perc, G = G_perc, `T` = T_perc) %>%
  gather(base, perc, A:`T`) %>%
  ggplot(aes(x = ctrl_index, y = perc, color = base)) +
  geom_point(alpha = 0.8) +
  labs(color = "Base") +
  scale_y_continuous(labels = scales::percent(seq(0,1, 0.1)), limits = c(0,1), breaks = seq(0,1, 0.1),
                     sec.axis = sec_axis(~.*5000, name = "Read Depth")
                     ) +
  scale_x_continuous(limits = c(0,600), breaks = seq(0,600,100)) +
  xlab("Position in Control Sequence") +
  ylab("Percent Base Composition") +
  ggtitle(paste0(sample_name, "\nRNASeq - Reditools")) +
  scale_color_manual(values = base_colors) +
  theme_bw(base_size = 18) +
  geom_line(aes(x = ctrl_index, y = Reads/5000), color = "#386cb0", alpha = 0.4)

# Sanger data in region of interest
# 1.17.19 would be good to modify to have ALL base information and either 0.05 alpha, or grey out non-sig or non-of-interest regions
sanger_plot = alt_df %>%
  dplyr::select(ctrl_index, A_perc:T_perc, A_sig:T_sig) %>% #, EI_sanger, A_sig_adjust:T_sig_adjust, A_pvalue:T_p_adjust
  dplyr::rename(A = A_perc, C = C_perc, G = G_perc, `T` = T_perc) %>%
  gather(base, perc, A:`T`) %>%
  dplyr::rename(A = A_sig, C = C_sig, G = G_sig, `T` = T_sig) %>%
  gather(sig_base, sig, A:`T`) %>%
  filter(base == sig_base) %>%
  dplyr::select(-sig_base) %>%
  ggplot(aes(x = ctrl_index, y = perc, color = base, alpha = sig)) +
  geom_point() +
  scale_y_continuous(labels = scales::percent(seq(0,1, 0.1)), limits = c(0,1), breaks = seq(0,1, 0.1)) +
  scale_x_continuous(limits = c(0,600), breaks = seq(0,600,100)) +
  xlab("Position in Control Sequence") +
  ylab("Percent Base Composition") +
  ggtitle(paste0(sample_name, "\nSanger sequencing - mEditR")) +
  labs(alpha = "Signficant?", color = "Base") +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_color_manual(values = base_colors) +
  theme_bw(base_size = 18)

##### Plot direct comparison
# Determine which ctrl reference positions are significant in the sanger file
#if(init_orient == "fwd"){test_edit = edit} else {test_edit = revcom(edit)}
test_edit = edit
# if(test_edit == "A") {sig_index = alt_df %>% filter(A_sig == TRUE) %>% .$ctrl_index} else
#   if(test_edit  == "C") {sig_index = alt_df %>% filter(C_sig == TRUE) %>% .$ctrl_index} else
#     if(test_edit  == "G") {sig_index = alt_df %>% filter(G_sig == TRUE) %>% .$ctrl_index} else
#       if(test_edit  == "T") {sig_index = alt_df %>% filter(T_sig == TRUE) %>% .$ctrl_index}

sig_index = alt_df$ctrl_index


# Make a filtered df from the sanger data
sanger_sample = alt_df %>%
  filter(ctrl_index %in% sig_index) %>%
  #dplyr::select(ctrl_index, A_adj_perc:T_adj_perc) %>%
  #dplyr::rename(A = A_adj_perc, C = C_adj_perc, G = G_adj_perc, `T` = T_adj_perc) %>%
  dplyr::select(ctrl_index, A_perc:T_perc, A_area:T_area, EI_sanger, A_pvalue:T_p_adjust) %>%
  dplyr::rename(A = A_perc, C = C_perc, G = G_perc, `T` = T_perc) %>%
  gather(base, sanger_perc, A:`T`)

# Make a filtered df from the NGS data
final_ngs_sample = ngs_sample %>%
  filter(ctrl_index %in% sig_index) %>%
  dplyr::select(Region, Position, Reference, Strand, ctrl_index, A_perc:T_perc, Reads, EI_ngs, gRef, gAlt, Ref, Alt) %>%
  dplyr::rename(A = A_perc, C = C_perc, G = G_perc, `T` = T_perc) %>%
  gather(base, ngs_perc, A:`T`)

# join filtered sanger and NGS data
final_data = inner_join(sanger_sample, final_ngs_sample) %>%
  filter(., base == wt | base == edit) %>%
  mutate(sample_file = sample_file) %>%
  mutate(phred = phred_cutoff, p_value = p_value, include = include)

# Plot the NGS data as a function of sanger data
comparison_plot = final_data %>%
  ggplot(aes(x = sanger_perc, y = ngs_perc)) +
  scale_x_continuous(labels = scales::percent(seq(0,1, 0.1)), limits = c(0,1), breaks = seq(0,1, 0.1)) +
  scale_y_continuous(labels = scales::percent(seq(0,1, 0.1)), limits = c(0,1), breaks = seq(0,1, 0.1)) +
  ggtitle(paste0(sample_name, "\nNGS - EditR Comparison")) +
  xlab("Sanger sequencing - mEditR") +
  ylab("RNASeq - Reditools") + 
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 18)

### Statistical analysis of the Sanger vs. NGS data
# Difference distributed about 0?
t_test = t.test(x = final_data$sanger_perc, y = final_data$ngs_perc, paired = TRUE)

# How well can sanger values be predicted by NGS values?
lm = lm(sanger_perc ~ ngs_perc + 0, data = final_data) %>% summary()

# How well can NGS data be predicted by Sanger data assuming Sanger data is the true value?
# I.e. how well does the purported truth predict the actual truth?
Rsq = rsq(y = final_data$ngs_perc, x = final_data$sanger_perc)

# produce a data frame with all the needed data for subsequent analyses

output_data = ctrl_df %>%
  dplyr::select(index, pre_pentanucleotide) %>%
  dplyr::rename(ctrl_index = index) %>%
  inner_join(alt_df, .) %>%
  mutate(tot_height = A_area + C_area + G_area + T_area) %>%
  dplyr::select(ctrl_index, index, ctrl_max_base, pre_pentanucleotide, tot_height, A_perc:T_perc, A_area:T_area, A_sig:T_sig, A_pvalue:T_p_adjust, motif, EI_sanger) %>%
  dplyr::mutate(edit_perc = !!sym(paste0(edit, "_perc"))) %>%
  dplyr::mutate(edit_sig = !!sym(paste0(edit, "_sig"))) %>%
  join(.,
       ngs_sample %>%
         dplyr::select(Region, Position, Reference, Strand, ctrl_index, ngs_index, Reads, A_perc, C_perc, G_perc, T_perc, EI_ngs, gRef, gAlt, Ref, Alt) %>%
         dplyr::rename(A_ngs_perc = A_perc, C_ngs_perc = C_perc, G_ngs_perc = G_perc, T_ngs_perc = T_perc)
         ) %>%
  as_tibble() %>%
  dplyr::mutate(edit_ngs_perc = !!sym(paste0(edit, "_ngs_perc"))) %>%
  mutate(wt = wt, edit = edit, fwd = fwd, p_value = p_value, phred_cutoff = phred_cutoff, sample_file = sample_file, ctrl_file = ctrl_file, ngs_file = ngs_file)

# Test alignment
test_alignment = pairwiseAlignment(samp_df %>% arrange(ctrl_index) %>%
                                .$ctrl_base_call %>%
                                paste0(., collapse = "") %>%
                                DNAString(),
                              ngs_sample %>%
                                arrange(ctrl_index) %>%
                                .$Reference %>%
                                paste0(., collapse = "") %>%
                                DNAString()
                              )
if(include) {return(
  list("ngs_plot" = ngs_plot,
       "sanger_plot" = sanger_plot,
       "comparison_plot" = comparison_plot,
       "rsq" = Rsq,
       "t_test" = t_test,
       "data" = final_data,
       "alignment" = test_alignment,
       "output_data" = output_data,
       "zaga_parameters" = zaga_parameters)
)} else {return(NULL)}
}

# Version for when sample gets skipped
AnalyzeNGS = function(
  directory,
  sample_file,# = parameters$sample_file[i],
  ctrl_file,# = parameters$ctrl_file[i],
  ngs_file,# = parameters$ngs_file[i],
  fwd,# = parameters$fwd[i],
  motif,# =parameters$motif[i],
  wt,# = parameters$wt[i], 
  edit,# = parameters$edit[i], 
  phred_cutoff,#= parameters$phred_cutoff[i],
  p_value = 0.01,#
  include # = parameters$include[i]
){
  tryCatch(analyzeNGS(directory = directory,
               sample_file = sample_file,# = parameters$sample_file[i],
               ctrl_file = ctrl_file,# = parameters$ctrl_file[i],
               ngs_file = ngs_file,# = parameters$ngs_file[i],
               fwd = fwd,# = parameters$fwd[i],
               motif = motif,# =parameters$motif[i],
               wt = wt,# = parameters$wt[i], 
               edit = edit,# = parameters$edit[i], 
               phred_cutoff = phred_cutoff,#= parameters$phred_cutoff[i],
               p_value = p_value,#
               include = include),
    error=function(e){message(paste0(sample_file, " failed to analyze."))})
  }
