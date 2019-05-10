### global.R for multiEditR

###########################################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

# Load sourcing scripts
# Don't need to run a sourcing line for an app, as it is all loaded into the global environment if it's in the app directory (?)
#source("/Users/kluesner/Desktop/Research/EditR/multiEditR/program/working_branch/dependencies.R")

runEditR = function(
  # Set parameters
  
  sample_file,
  ctrl_file,
  
  motif = "NAN", # Use IUPAC notation
  wt = "A", # Enter wt bases of interest with | separation
  edit = "G", # Enter edit of interest with | separation
  boi =  paste0(wt, "|", edit),
  phred_cutoff = 0.00001, ### 0.0001 seems good.
  trim = TRUE,
  p_value = 0.01,
  adjust_p = FALSE,
  bases = c("A", "C", "G", "T"),
  use_ctrl_seq = FALSE
  
  
  # Load sequencing files
  #sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_wt.ab1"
  #ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_ko.ab1"
  
  #sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/ME1/ME1_CmAG/CmAG_012_RP008_2018\ Sep\ 27.ab1"
  #ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/ME1/ME1_CmAG/CmAG/CmAG_001_RP008_2018 Sep 27.ab1"
  
  #sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_wt.ab1"
  #ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_ko.ab1"
){
  start = Sys.time()
  message("initialized.")
  # Make sangerseq objects
  # Need to flesh out the TRUE statement branch
  if(use_ctrl_seq)
  { #input_seq = abif_to_fastq(path = ctrl_file, cutoff = phred_cutoff)$seq
    input_seq = read_lines(ctrl_file)[2]
    init_ctrl_seq = input_seq
    ctrl_fastq = list()
    ctrl_fastq$seq = input_seq
    ctrl_df = data.frame(max_base = init_ctrl_seq %>% strsplit(., split = "") %>% unlist(),
                         base_call = init_ctrl_seq %>% strsplit(., split = "") %>% unlist()) %>%
      mutate(index = 1:NROW(max_base))
  } else 
  {
    # Generate ctrl sanger data frame
    # Generate ctrl primary basecalls
    ctrl_sanger = readsangerseq(ctrl_file)
    ctrl_df = make_ctrl_sanger_df(ctrl_sanger)
    init_ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
    ctrl_fastq = abif_to_fastq(path = ctrl_file, cutoff = phred_cutoff)
  }
  # set seed for reproducibility
  set.seed(666)
  
  # Make sangerseq object
  # Generate samp sanger data frame
  # Generate samp primary basecalls
  sample_sanger = readsangerseq(sample_file)
  sample_df = make_samp_sanger_df(sample_sanger, init_ctrl_seq)
  init_sample_seq = sample_df$primary_base_call %>% paste0(., collapse = "")
  # Genereate phred scores for ctrl and samp, trimming is built in using mott's algorithm
  sample_fastq = abif_to_fastq(path = sample_file, cutoff = phred_cutoff)
  
  # Align the both the ctrl and samp to their fastq filtered sequences
  # reasonable to assume phred scored sequence will always be smaller than the primary seq
  # use high gap penalty to force SNP alignment
  sample_alignment = pairwiseAlignment(pattern = sample_fastq$seq, subject = init_sample_seq)
  ctrl_alignment = pairwiseAlignment(pattern = ctrl_fastq$seq, subject = init_ctrl_seq)
  
  # Save unfiltered dataframes
  raw_sample_df = sample_df
  raw_ctrl_df = ctrl_df
  
  message("samples loaded.")
  
  # Filter dfs on high phred sequence
  sample_df %<>% 
    filter(index >= sample_alignment@subject@range@start) %<>%
    filter(index <= sample_alignment@subject@range@start + sample_alignment@subject@range@width - 1) %<>%
    mutate(post_filter_index = 1:NROW(index))
  ctrl_df %<>% 
    filter(index >= ctrl_alignment@subject@range@start) %<>%
    filter(index <= ctrl_alignment@subject@range@start + ctrl_alignment@subject@range@width - 1) %<>%
    mutate(post_filter_index = 1:NROW(index))
  
  # Regenerate primary basecalls
  ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
  sample_seq = sample_df$max_base %>% paste0(., collapse = "")
  
  # Save the pre-cross alignment dataframes
  pre_cross_align_sample_df = sample_df
  pre_cross_align_ctrl_df = ctrl_df
  
  message("samples filtered.")
  
  
  ### Bring samples together ###
  ### 01.07.19, if a ctrl sequence is used instead of a ctrl sequence, it would enter here.
  ### To use the context correction, would you still need to have the ctrl sequence to apply the GBM to ?
  # Align sample_seq to ctrl_seq
  trimmed_alignment = align_and_trim(sample_seq, ctrl_seq, min_continuity = 15)
  samp_alignment_seq = trimmed_alignment$alignment@pattern %>% as.character()
  ctrl_alignment_seq = trimmed_alignment$alignment@subject %>% as.character()
  
  # Align the trimmed sequences to the sequences from the data frame
  sample_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$pattern, subject = sample_seq, gapOpening = 1000, gapExtension = 1000, type = "local")
  ctrl_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$subject, subject = ctrl_seq, gapOpening = 1000, gapExtension = 1000, type = "local")
  
  ### Add predicted values from GBM
  sample_df = sample_df %>%
    dplyr::select(-(A_area:max_base_height), -post_filter_index, index) %>%
    bind_rows(., ., ., .) %>%
    mutate(base = rep(ACGT, each = NROW(index)/4), max_base = base) %>%
    mutate_if(is.character, nucleotide_factor) %>%
    mutate(eta = 1) %>%
    dplyr::select(-base) %>%
    spread(max_base, eta) %>%
    inner_join(sample_df, .) %>%
    dplyr::rename(A_eta = A, C_eta = C, G_eta = G, T_eta = `T`) %>%
    mutate(pred_height = {ifelse(max_base == "A", A_eta,
                                 ifelse(max_base == "C", C_eta,
                                        ifelse(max_base == "G", G_eta,
                                               ifelse(max_base == "T", `T_eta`,NA))))}) %>%
    dplyr::select(A_area:T_perc,
                  max_base, Tot.Area, index, max_base_height, post_filter_index, A_eta:T_eta, pred_height)
  
  
  message("GBM model applied.")
  
  # Filter dfs to aligned sequences
  sample_df %<>%
    filter(post_filter_index >= sample_trimmed_alignment@subject@range@start) %<>%
    filter(post_filter_index <= sample_trimmed_alignment@subject@range@start + sample_trimmed_alignment@subject@range@width - 1) %>%
    mutate(post_aligned_index = 1:NROW(index))
  ctrl_df %<>%
    filter(post_filter_index >= ctrl_trimmed_alignment@subject@range@start) %<>%
    filter(post_filter_index <= ctrl_trimmed_alignment@subject@range@start + ctrl_trimmed_alignment@subject@range@width - 1) %>%
    mutate(post_aligned_index = 1:NROW(index))
  
  message("sample aligned.")
  
  tmp = sample_df
  
  
  ### Generate post-filter, post-aligned sequence
  ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
  samp_seq = sample_df$max_base %>% paste0(., collapse = "")
  
  ### Save a df for NGS analysis
  pre_aligned_sample_df = sample_df
  
  ### Filter out any base positons that have an indel
  samp_indel = samp_alignment_seq %>% gregexpr("-", .) %>% unlist
  ctrl_indel = ctrl_alignment_seq %>% gregexpr("-", .) %>% unlist
  
  ### Code added 10.27.18, as when there is no indels in one sample it returned a -1 instead of 0, which introduced an error
  if(samp_indel == -1){samp_indel = 0} else {}
  if(ctrl_indel == -1){ctrl_indel = 0} else {}
  
  ctrl_df = ctrl_df %>%
    mutate(.,
           indel_filter = ctrl_alignment_seq %>%
             subchar(., samp_indel, "_") %>%
             gsub("-", "", .) %>%
             strsplit(x = ., split = "") %>%
             .[[1]]
    ) %>%
    filter(indel_filter != "_") %>%
    dplyr::select(-indel_filter) %>%
    mutate(filtered_index = 1:NROW(max_base))
  
  sample_df = sample_df %>%
    mutate(.,
           indel_filter = samp_alignment_seq %>%
             subchar(., ctrl_indel, "_") %>%
             gsub("-", "", .) %>%
             strsplit(x = ., split = "") %>%
             .[[1]]
    ) %>%
    filter(indel_filter != "_") %>%
    dplyr::select(-indel_filter) %>%
    mutate(filtered_index = 1:NROW(max_base)) %>%
    mutate(ctrl_post_aligned_index = ctrl_df$post_aligned_index)
  
  message("Indels removed.")
  ### join the initial ctrl index to the sample to give a reference
  sample_df = ctrl_df %>%
    dplyr::select(post_aligned_index, index) %>%
    dplyr::rename(ctrl_post_aligned_index = post_aligned_index, ctrl_index = index) %>%
    inner_join(., sample_df) #%>%
  #dplyr::select(everything(), ctrl_index, ctrl_post_aligned_index, pre_trinucleotide, post_trinucleotide)
  
  ### Assign the sample and ctrl file names
  sample_df %<>% mutate(sample_file = sample_file, ctrl_file = ctrl_file)
  ctrl_df %<>% mutate(sample_file = sample_file, ctrl_file = ctrl_file)
  
  ### reassign the sequences post filtering
  samp_df_seq = sample_df$max_base %>% paste0(., collapse = "")
  ctrl_df_seq = ctrl_df$max_base %>% paste0(., collapse = "")
  
  message("Indices adjusted.")
  ### ENTER MOTIF ISOLATION ###
  # Will want to be compare this to post alignment seq and index, but before indel removal
  # ctrl_seq will have the same indexing as post_aligned_index, but the indexing is not the same across the two
  # Will need to figure out how to be able to pull the same indexing out of both of them
  # Could take the post_aligned_index from ctrl and apply it to the sample
  
  # Align the motif of interest to the ctrl_seq
  motif_alignment = matchPattern(pattern = DNAString(motif), subject = DNAString(ctrl_seq), fixed = FALSE)
  n_alignments = motif_alignment@ranges %>% length()
  
  motif_positions = mapply(FUN = seq,
                           from = motif_alignment@ranges@start,
                           to = (motif_alignment@ranges@start + nchar(motif) - 1)) %>% as.vector()
  names(motif_positions) = rep(x = c(1:n_alignments), each = nchar(motif))
  
  # Append the sequences from the ctrl df to the sample df
  sample_df = sample_df %>%
    mutate(ctrl_max_base = ctrl_df$max_base, ctrl_base_call = ctrl_df$base_call) %>%
    # Add a column for the percent base for the ctrl basecall
    mutate(ctrl_max_base_perc = {ifelse(ctrl_max_base == "A", A_perc,
                                        ifelse(ctrl_max_base == "C", C_perc,
                                               ifelse(ctrl_max_base == "G", G_perc,
                                                      ifelse(ctrl_max_base == "T", T_perc,0))))})
  
  # Generate null and alternative samples for distribution
  # Perform for both the sample df
  sample_null = sample_df %>% filter(!(ctrl_post_aligned_index %in% motif_positions)) # This dataframe consists of all bases in the sample where the motif is not found in the ctrl sequence
  sample_alt = sample_df %>% filter(ctrl_post_aligned_index %in% motif_positions) # This dataframe consists of all bases in the sample where the motif is found in the ctrl sequence
  
  # Find all potential events of significant noise
  filtered_sample_alt = sample_alt %>%
    filter(grepl(wt, ctrl_max_base)) %>% # Use the ctrl wt base for determining data
    dplyr::rename(A = A_area, C = C_area, G = G_area, `T` = T_area) %>%
    gather(base, height, A:`T`) %>%
    filter(grepl(edit, base)) # Filter out hypothetical mutations that are not of interest
  
  message("Motifs of interest mapped and subsetted.")
  
  # Adjust p-value
  n_comparisons = NROW(filtered_sample_alt)
  # Holm-sidak correction, may need to use the smirnov correction for family-wise error rates
  # Will need to read original paper and cite appropriately
  if(adjust_p == TRUE) {p_adjust = 1-((1-p_value)^(1/n_comparisons))} else {p_adjust = p_value}
  
  # Generate zG models for each base
  # uses the sample_null to calculate
  zaga_parameters = make_ZAGA_df(sample_null, p_adjust = p_adjust)
  critical_values = zaga_parameters$crit
  
  ### Find significant edits and then apply the GBM adjustments
  # Determine which values are significant
  # Keep significant values and replace all n.s. values with NA
  # Use the significant values to calculated an adjusted height using formula 1
  # Find all significant edits
  output_sample_alt = gbm_adjust(sample_alt, wt, boi, motif, sample_file, critical_values)
  output_sample_null = sample_null %>%
    dplyr::select(ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc) %>%
    mutate(motif = motif, sample_file = sample_file)
  
  ### Make plots
  plot1 = raw_sample_df %>%
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
  
  plot2 = sample_df %>%
    ggplot(aes(x = index, y = (100-ctrl_max_base_perc))) +
    geom_line() +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0)) +
    scale_x_continuous(limits = c(0, max(raw_sample_df$index)), expand = c(0,0)) +
    #geom_bar(data = sample_alt, aes(x = index, y = 100), stat = "identity", fill = "lightblue", color = "lightblue") +
    geom_rect(xmin = min(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(pre_cross_align_sample_df$index), ymax = 100,
              fill = "white", color = "black", alpha = 0) +
    xlab("Position in sample file") +
    ylab("Percent noise in WT basecall") +
    geom_hline(yintercept = mean(pre_cross_align_sample_df$max_base_height), color = "darkred", size = 1) +
    geom_line() +
    geom_point(data = output_sample_alt %>% mutate(sig = {ifelse(sig, "Signficant", "Non-significant")}),
               aes(x = index, y = 100-ctrl_max_base_perc, color = sig), size = 3) +
    geom_rect(xmin = 0, ymin = 0, xmax = min(pre_cross_align_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    geom_rect(xmin = max(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(raw_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    theme_classic(base_size = 18) +
    labs(color = "Potential edits") +
    theme(legend.position = c(1/6, 5/6))
  
  editing_data = output_sample_alt %>%
    dplyr::select(index, A_perc:T_perc) %>%
    gather(base, perc, A_perc:`T_perc`) %>%
    inner_join(., output_sample_alt %>%
                 dplyr::select(index, A_sig:T_sig) %>%
                 gather(base_sig, sig, A_sig:T_sig)
    ) %>%
    mutate(base = gsub("_perc", "", base), base_sig = gsub("_sig", "", base_sig)) %>%
    mutate(perc = perc*100, sig = {ifelse(sig, "Significant", "Non-significant")}) %>%
    filter(base == base_sig) %>%
    dplyr::select(-base_sig) %>%
    filter(grepl(pattern = edit, x = base))
  
  plot3 = editing_data %>%
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
  
  output_editing_data = editing_data %>%
    group_by(base, sig) %>%
    dplyr::summarize(Mean = mean(perc), Max = max(perc), Min = min(perc), SD = sd(perc), N = length(perc)) %>%
    dplyr::rename(Base = base, Significance = sig)
  
  output_zaga_parameters = zaga_parameters %>% filter(grepl(edit, Base))
  output_sample_alt = output_sample_alt %>%
    dplyr::select(-A_area:-T_area)
  # Reset start directory
  message(paste0(round(Sys.time() - start, 2), " seconds elapsed."))
  
  return(list(output_editing_data, plot1, plot2, plot3, output_zaga_parameters, output_sample_alt))
  
}
