### Working script for multiEditR

###########################################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

# Load sourcing scripts

source("/Users/kluesner/Desktop/Research/EditR/multiEditR/program/working_branch/dependencies.R")

extract_df = function(
# Set parameters
  
sample_file,
ctrl_file,

motif = "A", # Use IUPAC notation
wt = "A", # Enter wt bases of interest with | separation
edit = "G", # Enter edit of interest with | separation
boi = paste0(wt, "|", edit),
phred_cutoff = 0.001, ### 0.0001 seems good.
trim = TRUE,
p_value = 0.01,
bases = c("A", "C", "G", "T"),
use_ctrl_seq = T,
fwd_orientation = TRUE

# Load sequencing files
# sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_wt.ab1"
# ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_ko.ab1"

#sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/ME1/ME1_CmAG/CmAG_012_RP008_2018\ Sep\ 27.ab1"
#ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/ME1/ME1_CmAG/CmAG/CmAG_001_RP008_2018 Sep 27.ab1"

#sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_wt.ab1"
#ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_ko.ab1"
){
message("initialized.")
start = Sys.time()
# Make sangerseq objects
# Need to flesh out the TRUE statement branch
if(use_ctrl_seq)
{ #input_seq = abif_to_fastq(path = ctrl_file, cutoff = phred_cutoff)$seq
  input_seq = read_lines(ctrl_file)[2]
  init_ctrl_seq = input_seq
  if(fwd_orientation){} else {init_ctrl_seq = revcom(init_ctrl_seq)}
  ctrl_fastq = list()
  ctrl_fastq$seq = input_seq
  ctrl_df = data.frame(max_base = init_ctrl_seq %>% base::strsplit(., split = "") %>% unlist(),
                       base_call = init_ctrl_seq %>% base::strsplit(., split = "") %>% unlist()) %>%
    mutate(index = 1:NROW(max_base)) %>% 
    mutate(pre_5 = lag(base_call, n = 5),
           pre_4 = lag(base_call, n = 4),
           pre_3 = lag(base_call, n = 3),
           pre_2 = lag(base_call, n = 2),
           pre_1 = lag(base_call, n = 1),
           pre_pentanucleotide = paste0(pre_5, pre_4, pre_3, pre_2, pre_1)) 
} else 
{
  # Generate ctrl sanger data frame
  # Generate ctrl primary basecalls
  ctrl_sanger = readsangerseq(ctrl_file)
  ctrl_df = make_ctrl_sanger_df(ctrl_sanger)
  init_ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
  if(fwd_orientation){} else {init_ctrl_seq = revcom(init_ctrl_seq)}
  ctrl_fastq = abif_to_fastq(path = ctrl_file, cutoff = phred_cutoff)
}
  
# IF the sample needs to be reversed, then reverse complement the ctrl_seq

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
# sample_df = sample_df %>%
#   dplyr::select(A_area:T_perc,
#                 max_base, Tot.Area, index, max_base_height, post_filter_index)


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
sample_df %<>% mutate(ctrl_max_base = ctrl_df$max_base, ctrl_base_call = ctrl_df$base_call)

# calculate an editing index

# sample_df %<>%
#   group_by(ctrl_max_base) %<>%
#   mutate(EI = sum(!! sym(paste0(edit, "_perc"))) / sum(!! sym(edit)) + sum(!! sym(wt)))

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
# if(adjust_p) {p_adjust = p.adjust(p_value, method = "holm", n = n_comparisons)} else {p_adjust = p_value}

# Generate zG models for each base
# uses the sample_null to calculate
zaga_parameters = make_ZAGA_df(sample_null, p_adjust = p_value) %>%
  mutate(sample_file = sample_file)
critical_values = zaga_parameters$crit

### Find significant edits and then apply the GBM adjustments
# Determine which values are significant
# Keep significant values and replace all n.s. values with NA
# Use the significant values to calculated an adjusted height using formula 1
# Find all significant edits
output_sample_alt = pvalue_adjust(sample_alt, wt, boi, motif, sample_file, critical_values, zaga_parameters, p_value)

# Create Editing index for output_sample_alt
# 02.02.2020
output_sample_alt = output_sample_alt %>%
  dplyr::group_by(ctrl_max_base) %>% # technically shouldn't change anything as it's the same across samples
  dplyr::mutate(AEI_sanger = (sum(G_perc) / (sum(A_perc) + sum(G_perc))),
         CEI_sanger = (sum(T_perc) / (sum(C_perc) + sum(T_perc))),
         GEI_sanger = (sum(A_perc) / (sum(G_perc) + sum(A_perc))),
         TEI_sanger = (sum(C_perc) / (sum(T_perc) + sum(C_perc)))
  ) %>%
  ungroup() %>%
  # Keep only the EI_index for each reference base
  gather(EI_base, EI_sanger, AEI_sanger:TEI_sanger) %>%
  mutate(EI_base = gsub("EI_sanger", "", EI_base)) %>%
  filter(EI_base == ctrl_max_base) %>%
  dplyr::select(-EI_base)


output_sample_null = sample_null %>%
  dplyr::select(ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc) %>%
  mutate(motif = motif, sample_file = sample_file)

output = list(sample_df, output_sample_alt, init_ctrl_seq, ctrl_df, zaga_parameters)

# Reset start directory
message(paste0(round(Sys.time() - start, 2), " seconds elapsed."))

return(output)

}

#tmp = output_sample_alt %>%
#  dplyr::select(ctrl_index, A_perc:T_perc, A_adj_perc:T_adj_perc) %>%
#  dplyr::rename(A = A_perc, C = C_perc, G = G_perc, `T` = T_perc) %>%
#  gather(base, perc, A:`T`) %>%
 # dplyr::rename(A = A_adj_perc, C = C_adj_perc, G = G_adj_perc, `T` = T_adj_perc) %>%
#  gather(adj_base, adj_perc, A:`T`) %>%
#  filter(base == adj_base) %>%
#  dplyr::select(-adj_base) %>%
#  gather(measure, perc, perc:adj_perc)

#tmp %>%
 # ggplot(aes(x = ctrl_index, y = perc, color = base)) +
 # geom_point() +
 # theme_bw(base_size = 12) +
 # facet_wrap(.~measure)

### Algorithm
# load the sanger sequencing files
# Extract primary basecalls
# Use sang_to_df() to make a df
#   Future work would normalize sanger sequences -- maybe look at TIDE?
#   Trim the df
#   Possibly trim on phred or something?
# Align primary basecalls
#   Use some level of mismatch tolerance
# Subset aligned region of dfs
# Delete the ctrl df from the sample df

