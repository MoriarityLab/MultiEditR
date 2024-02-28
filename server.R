### server.R for multiEditR

###########################################################################################
# Copyright (C) 2020-2021 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

##### Example data ----------------------------------------------------------------
example_sample = "RP272_cdna_wt.ab1"
example_control = "RP272_cdna_ko.ab1"
example_control_fasta = "RP272_cdna_ko.fasta"
example_motif = "YAR"
example_wt = "A"
example_edit = "G"
example_use_ctrl_seq = FALSE
example_phred = 0.0001

##### Server function -------------------------------------------------------------
shinyServer(
  function(input, output, session) {
    
    ##### Establish reactive input functions --------------------------------------
    # These reactive functions are used as the inputs for the analysis
    
    sample.Reactive = reactive({
      # if else things to handle example data loading
      if(input$use_example) {
        return(example_sample)
      } else if(!is.null(input$sample_file)) {
        return(input$sample_file$datapath)
      } else return(validate(
        need(input$sample_file, "Please upload your sanger sequence sample file")
      ))
    })
    
    ctrl.Reactive = reactive({
      # if else things to handle example data loading
      if(input$use_example) {
        if(input$use_fasta){
          return(example_control_fasta)
        } else
          return(example_control)
      } else if(!is.null(input$ctrl_file)) {
        return(input$ctrl_file$datapath)
      } else return(validate(
        need(input$ctrl_file, "Please upload your control sequence file")
      ))
    })

    motif.Reactive = reactive({
      # if else things to handle example data loading
      if(input$use_example & (input$motif == "")) {
        return(example_motif)
      } else if(input$motif != "") {
        return(input$motif)
      } else return(validate(
        need(input$motif, "Please enter a motif of interest")
      ))
    })
    
    wt.Reactive = reactive({
      # if else things to handle example data loading
      if(input$use_example & (input$wt == "")) {
        return(example_wt)
      } else if(input$wt != "") {
        return(input$wt)
      } else return(validate(
        need(input$wt, "Please enter a wt base of interest")
      ))
    })
    
    edit.Reactive = reactive({
      # if else things to handle example data loading
      if(input$use_example & (input$wt == "")) {
        return(example_edit)
      } else if(input$edit != "") {
        return(input$edit)
      } else return(validate(
        need(input$edit, "Please enter an edited base of interest")
      ))
    })
    
    ##### Begin computations ----------------------------------------------------
    
    #### Establish control sample information -----------------------------------
    # Reactive ctrl sanger object, only if using a fasta is NOT specified
    ctrl_sanger.Reactive = reactive(
      if(!input$use_fasta)
      {readsangerseq(ctrl.Reactive())} else
      {NULL})
    
    # Reactive ctrl df object, handled differently for a fasta vs sanger input
    ctrl_df.Reactive = reactive({
      if(input$use_fasta)
      {make_ctrl_fasta_df(ctrl.Reactive())} else
      {make_ctrl_sanger_df(ctrl_sanger.Reactive())}
    })
    
    # Reactive ctrl initial sequence
    init_ctrl_seq.Reactive = reactive({
      ctrl_df.Reactive()$base_call %>% paste0(., collapse = "")
    })
    
    # Reactive ctrl fastq
    ctrl_fastq.Reactive = reactive({
      if(input$use_fasta)
      {ctrl_fastq = list(); ctrl_fastq$seq = init_ctrl_seq.Reactive(); return(ctrl_fastq)} else
        # Genereate phred scores for ctrl and samp, trimming is built in using mott's algorithm
      {abif_to_fastq(path = ctrl.Reactive(), cutoff = input$phred)}
    })
    
    # Reactive ctrl alignment
    ctrl_alignment.Reactive = reactive({
      pairwiseAlignment(pattern = ctrl_fastq.Reactive()$seq, subject = init_ctrl_seq.Reactive())
    })
    
    #### Establish treatment sample information -----------------------------------
    # Reactive sample sanger object
    sample_sanger.Reactive = reactive({readsangerseq(sample.Reactive())})
    
    # Reactive sample df object
    sample_df.Reactive = reactive({
      make_samp_sanger_df(sample_sanger.Reactive(), init_ctrl_seq.Reactive())
    })
    
    # Reactive sample initial sequence
    init_sample_seq.Reactive = reactive({
      sample_df.Reactive()$primary_base_call %>% paste0(., collapse = "")
    })
    
    # Reactive sample fastq object
    # Genereate phred scores for ctrl and samp, trimming is built in using mott's algorithm
    sample_fastq.Reactive = reactive({
      abif_to_fastq(path = sample.Reactive(), cutoff = input$phred)
    })
    
    # Reactive sample alignment
    sample_alignment.Reactive = reactive({
      pairwiseAlignment(pattern = sample_fastq.Reactive()$seq, subject = init_sample_seq.Reactive())
    })
    
    #### Filter dataframes on phred score -------------------------------------
    # Reactive object containing the phred trimmed data
    filteredData.Reactive = reactive({
      
      # trim sample dataframe on regions of phred filtering
      sample_df = sample_df.Reactive() %>% 
        filter(index >= sample_alignment.Reactive()@subject@range@start) %>%
        filter(index <= sample_alignment.Reactive()@subject@range@start +
                 sample_alignment.Reactive()@subject@range@width - 1) %>%
        mutate(post_filter_index = 1:NROW(index))
      
      # trim control dataframe on regions of phred filtering
      ctrl_df = ctrl_df.Reactive() %>% 
        filter(index >= ctrl_alignment.Reactive()@subject@range@start) %>%
        filter(index <= ctrl_alignment.Reactive()@subject@range@start +
                 ctrl_alignment.Reactive()@subject@range@width - 1) %>%
        mutate(post_filter_index = 1:NROW(index))
      
      # Regenerate primary basecalls
      ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
      sample_seq = sample_df$max_base %>% paste0(., collapse = "")
      
      # Save the pre-cross alignment dataframes
      pre_cross_align_sample_df = sample_df
      pre_cross_align_ctrl_df = ctrl_df
      
      # Return a list containing the information for downstream computation
      return(list("sample_df" = sample_df,
                  "ctrl_df" = ctrl_df,
                  "ctrl_seq" = ctrl_seq,
                  "sample_seq" = sample_seq,
                  "pre_cross_align_sample_df" = pre_cross_align_sample_df)
      )
    })
    
    #### Trim dataframes by aligning sample to control ------------------------
    # Reactive object with trimmed sample and control dataframes

    trimmedData.Reactive = reactive({
      
      # Align sample_seq to ctrl_seq
      trimmed_alignment = align_and_trim(
        filteredData.Reactive()$sample_seq,
        filteredData.Reactive()$ctrl_seq,
        min_continuity = 15)
      
      samp_alignment_seq = trimmed_alignment$alignment@pattern %>% as.character()
      ctrl_alignment_seq = trimmed_alignment$alignment@subject %>% as.character()
      
      # Align the trimmed sequences to the sequences from the data frame
      sample_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$pattern,
                                                   subject = filteredData.Reactive()$sample_seq,
                                                   gapOpening = 1000,
                                                   gapExtension = 1000,
                                                   type = "local")
      
      ctrl_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$subject,
                                                 subject = filteredData.Reactive()$ctrl_seq,
                                                 gapOpening = 1000,
                                                 gapExtension = 1000,
                                                 type = "local")
      
      # Add predicted values from GBM
      # This will likely be deprecated (07.12.2019)
      sample_df = filteredData.Reactive()$sample_df %>%
        dplyr::select(-(A_area:max_base_height), -post_filter_index, index) %>%
        bind_rows(., ., ., .) %>%
        mutate(base = rep(ACGT, each = NROW(index)/4), max_base = base) %>%
        mutate_if(is.character, nucleotide_factor) %>%
        mutate(eta = 1) %>%
        dplyr::select(-base) %>%
        spread(max_base, eta) %>%
        inner_join(filteredData.Reactive()$sample_df, .) %>%
        dplyr::rename(A_eta = A, C_eta = C, G_eta = G, T_eta = `T`) %>%
        mutate(pred_height = {ifelse(max_base == "A", A_eta,
                                     ifelse(max_base == "C", C_eta,
                                            ifelse(max_base == "G", G_eta,
                                                   ifelse(max_base == "T", `T_eta`,NA))))}) %>%
        dplyr::select(A_area:T_perc,
                      max_base, Tot.Area, index, max_base_height, post_filter_index, A_eta:T_eta, pred_height)
      
      # Filter the sample dataframe on the end trimming alignment
      sample_df %<>%
        filter(post_filter_index >= sample_trimmed_alignment@subject@range@start) %<>%
        filter(post_filter_index <= sample_trimmed_alignment@subject@range@start + sample_trimmed_alignment@subject@range@width - 1) %>%
        mutate(post_aligned_index = 1:NROW(index))
      
      # Filter the control dataframe on the end trimming alignment
      ctrl_df = filteredData.Reactive()$ctrl_df %>%
        # error fixed 5.21.2021, use of %<>% not permitted
        filter(post_filter_index >= ctrl_trimmed_alignment@subject@range@start) %>%
        filter(post_filter_index <= ctrl_trimmed_alignment@subject@range@start + ctrl_trimmed_alignment@subject@range@width - 1) %>%
        mutate(post_aligned_index = 1:NROW(index))

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
      
      # join the initial ctrl index to the sample to give a reference
      sample_df = ctrl_df %>%
        dplyr::select(post_aligned_index, index) %>%
        dplyr::rename(ctrl_post_aligned_index = post_aligned_index, ctrl_index = index) %>%
        inner_join(., sample_df)
      
      # Assign the sample and ctrl file names
      sample_df %<>% mutate(sample_file = sample.Reactive(), ctrl_file = ctrl.Reactive())
      ctrl_df %<>% mutate(sample_file = sample.Reactive(), ctrl_file = ctrl.Reactive())
      
      # Return a list containing the information for downstream computation
      return(list("sample_df" = sample_df,
                  "ctrl_df" = ctrl_df,
                  "ctrl_seq" = ctrl_seq,
                  "pre_aligned_sample_df" = pre_aligned_sample_df)
      )
    })
    
    #### Isolate motifs and begin statistical tests ---------------------------
    
    statisticalAnalysis.Reactive = reactive({
      
      # Align the motif of interest to the ctrl_seq
      # In reality this is matching, not alignment
      motif_alignment = matchPattern(pattern = DNAString(motif.Reactive()),
                                     subject = DNAString(trimmedData.Reactive()$ctrl_seq),
                                     fixed = FALSE)
      
      # Count the number of motif matches found
      n_alignments = motif_alignment@ranges %>% length()
      
      # Find the motif positions
      motif_positions = mapply(FUN = seq,
                               from = motif_alignment@ranges@start,
                               to = (motif_alignment@ranges@start + nchar(motif.Reactive()) - 1)) %>%
        as.vector()
      
      # Name the motif positions
      names(motif_positions) = rep(x = c(1:n_alignments), each = nchar(motif.Reactive()))
      
      motif_positions = data.frame(ctrl_post_aligned_index = motif_positions, motif_id = names(motif_positions))
      
      
      
      # Append the sequences from the ctrl df to the sample df
      sample_df = trimmedData.Reactive()$sample_df %>%
        mutate(ctrl_max_base = trimmedData.Reactive()$ctrl_df$max_base,
               ctrl_base_call = trimmedData.Reactive()$ctrl_df$base_call
        ) %>%
        # Add a column for the percent base for the ctrl basecall
        mutate(ctrl_max_base_perc = {ifelse(ctrl_max_base == "A", A_perc,
                                            ifelse(ctrl_max_base == "C", C_perc,
                                                   ifelse(ctrl_max_base == "G", G_perc,
                                                          ifelse(ctrl_max_base == "T", T_perc,0))))})
      
      # Generate null and alternative samples for distribution
      # Perform for both the sample df
      
      # This dataframe consists of all bases in the sample where the motif is not found in the ctrl sequence
      sample_null = sample_df %>%
        filter(!(ctrl_post_aligned_index %in% motif_positions$ctrl_post_aligned_index)) 
      
      # This dataframe consists of all bases in the sample where the motif is found in the ctrl sequence
      sample_alt = motif_positions %>%
        inner_join(., sample_df)
      
      # Find all potential events of significant noise
      filtered_sample_alt = sample_alt %>%
        filter(grepl(wt.Reactive(), ctrl_max_base)) %>% # Use the ctrl wt base for determining data
        dplyr::rename(A = A_area, C = C_area, G = G_area, `T` = T_area) %>%
        gather(base, height, A:`T`) %>%
        filter(grepl(edit.Reactive(), base)) # Filter out hypothetical mutations that are not of interest
      
      # Adjust p-value
      n_comparisons = NROW(filtered_sample_alt)
      # Holm-sidak correction, may need to use the smirnov correction for family-wise error rates
      # Will need to read original paper and cite appropriately
      if(input$adjust_p == TRUE)
      {p_adjust = 1-((1-input$p_value)^(1/n_comparisons))} else
      {p_adjust = input$p_value}
      
      # Generate zG models for each base
      # uses the sample_null to calculate
      zaga_parameters = make_ZAGA_df(sample_null, p_adjust = p_adjust)
      critical_values = zaga_parameters$crit
      
      ### Find significant edits and then apply the GBM adjustments
      # Determine which values are significant
      # Keep significant values and replace all n.s. values with NA
      # Use the significant values to calculated an adjusted height using formula 1
      # Find all significant edits
      # Will need to deprectate the GBM adjustment eventually
      output_sample_alt = gbm_adjust(sample_alt,
                                     wt.Reactive(),
                                     edit.Reactive(),
                                     motif.Reactive(),
                                     sample.Reactive(),
                                     critical_values)
      
      # Calculate the EI_index for each base
      sanger_EI = output_sample_alt %>%
        # group_by(ctrl_max_base) %>%
        # dplyr::summarise(sum_A_perc = sum(A_perc),
        #                  sum_C_perc = sum(C_perc),
        #                  sum_G_perc = sum(G_perc),
        #                  sum_T_perc = sum(T_perc)) %>%
        # dplyr::ungroup() %>%
        # dplyr::mutate(AEI_sanger = (sum_G_perc / ((!! sym(paste0("sum_", wt.Reactive(), "_perc"))) + sum_G_perc)),
        #               CEI_sanger = (sum_T_perc / ((!! sym(paste0("sum_", wt.Reactive(), "_perc"))) + sum_T_perc)),
        #               GEI_sanger = (sum_A_perc / ((!! sym(paste0("sum_", wt.Reactive(), "_perc"))) + sum_A_perc)),
        #               TEI_sanger = (sum_C_perc / ((!! sym(paste0("sum_", wt.Reactive(), "_perc"))) + sum_C_perc))
        # )
        
        dplyr::group_by(ctrl_max_base) %>% # technically shouldn't change anything as it's the same across samples
        dplyr::mutate(AEI_sanger = (sum(G_perc) / (sum(!! sym(paste0(wt.Reactive(), "_perc"))) + sum(G_perc))),
                      CEI_sanger = (sum(T_perc) / (sum(!! sym(paste0(wt.Reactive(), "_perc"))) + sum(T_perc))),
                      GEI_sanger = (sum(A_perc) / (sum(!! sym(paste0(wt.Reactive(), "_perc"))) + sum(A_perc))),
                      TEI_sanger = (sum(C_perc) / (sum(!! sym(paste0(wt.Reactive(), "_perc"))) + sum(C_perc)))
        ) %>%
        dplyr::select(AEI_sanger:TEI_sanger) %>%
        dplyr::distinct() %>%
        ungroup() %>%
        # Keep only the EI_index for each reference base
        gather(EI_base, EI_sanger, AEI_sanger:TEI_sanger) %>%
        mutate(EI_base = gsub("EI_sanger", "", EI_base)) %>%
        dplyr::select(EI_base, EI_sanger)  %>%
        dplyr::rename(Edit =  EI_base, MEI = EI_sanger) %>%
        filter(grepl(wt.Reactive(), Edit)) %>%
        distinct()
      
      return(
        list("sample_df" = sample_df,
             "sample_alt" = sample_alt,
             "output_sample_alt" = output_sample_alt,
             "zaga_parameters" = zaga_parameters,
             "motif_positions" = motif_positions,
             "sanger_EI" =  sanger_EI)
      )
    })
    
    #### Calculate values for plotting ----------------------------------------
    
    editingData.Reactive = reactive({
      calculateEditingData(
        statisticalAnalysis.Reactive()$output_sample_alt,
        edit.Reactive())
    })
    
    #### Render outputs -------------------------------------------------------
    
    # ### Checking ctrl data ----------------------------------------------------
    # # ctrl sanger names
    # output$ctrl_sanger = renderPrint({ctrl_sanger.Reactive()})
    # 
    # # ctrl df
    # output$ctrl_df = renderPrint({head(ctrl_df.Reactive())})
    # 
    # # init_ctrl_seq
    # output$init_ctrl_seq = renderPrint({init_ctrl_seq.Reactive()})
    # 
    # # ctrl_fastq
    # output$ctrl_fastq = renderPrint({as.character(ctrl_fastq.Reactive()$seq)})
    # 
    # # ctrl_alignment
    # output$ctrl_alignment = renderPrint(ctrl_alignment.Reactive())
    # 
    # ### Checking sample data ----------------------------------------------------
    # # sample sanger names
    # output$sample_sanger = renderPrint({sample_sanger.Reactive()})
    # 
    # # sample df
    # output$sample_df = renderPrint({head(sample_df.Reactive())})
    # 
    # # init_sample_seq
    # output$init_sample_seq = renderPrint({init_sample_seq.Reactive()})
    # 
    # # sample_fastq
    # output$sample_fastq = renderPrint({as.character(sample_fastq.Reactive()$seq)})
    # 
    # # sample_alignment
    # output$sample_alignment = renderPrint(sample_alignment.Reactive())
    # 
    # ### Checking computation data ----------------------------------------------------
    # # Phred filtered data
    # output$filteredData = renderPrint({filteredData.Reactive()})
    # 
    # # End trimmed data
    # output$trimmedData = renderPrint({trimmedData.Reactive()})
    # 
    # # Statistical data
    output$rawData = DT::renderDataTable({
      statisticalAnalysis.Reactive()$output_sample_alt %>%
        mutate(A = round(A_perc*100),
               C = round(C_perc*100),
               G = round(G_perc*100),
               `T` = round(T_perc*100)
        ) %>%
        dplyr::select(index, ctrl_index, max_base, A:`T`, A_sig:T_sig, motif, sample_file)
      })
    
    # Render .csv file from dataTable
    output$downloadRawData <- downloadHandler(
      filename = function() {
        paste(sample.Reactive(), "rawData.tsv", sep = "_") %>%
          gsub("[.]ab1|[.]AB1|[.]abif|[.]ABIF", "", .)
      },
      content = function(file) {
        readr::write_tsv(statisticalAnalysis.Reactive()$output_sample_alt %>%
                           mutate(A = round(A_perc*100),
                                  C = round(C_perc*100),
                                  G = round(G_perc*100),
                                  `T` = round(T_perc*100)
                           ) %>%
                           dplyr::select(index, ctrl_index, max_base, A:`T`, A_sig:T_sig, motif, sample_file)
                         , file)
      }
    )
    
    ### Plotting data ----------------------------------------------------------------
    
    ## Plot percent signal of raw sample
    output$plotRawSample = renderPlot({
      plotRawSample(
        sample_df.Reactive(),
        statisticalAnalysis.Reactive()$sample_alt,
        filteredData.Reactive()$pre_cross_align_sample_df
      )
    })
    
    ## Plot percent noise of trimmed sample
    output$plotTrimmedSample = renderPlotly({
      plotTrimmedSample(
        statisticalAnalysis.Reactive()$sample_df,
        filteredData.Reactive()$pre_cross_align_sample_df,
        statisticalAnalysis.Reactive()$output_sample_alt,
        sample_df.Reactive(),
        statisticalAnalysis.Reactive()$sample_alt
      )
    })
    
    
    
    ## Plot summary of editing data
    output$plotEditingData = renderPlot({
      plotEditingData(editingData.Reactive())
    })
    
    ## Render table summarizing editing data
    output$tableEditingData = renderTable({
      tableEditingData(editingData.Reactive())
    })
    
    ## Render table summarizing editing data
    output$tableEditingIndexData = renderTable({
      statisticalAnalysis.Reactive()$sanger_EI
    }, digits = 4) # digits = 4 dictates enough significant figures
    
    ## Render table with ZAGA parameters for bases of interest
    output$tableZAGAParameters = renderTable({
      statisticalAnalysis.Reactive()$zaga_parameters %>%
        filter(grepl(edit.Reactive(), Base))
    })
    
    sampleIndices.Reactive = reactive({
      sampleIndices = {
        if(!input$usePoint & is.null(input$plotTrimmedSample2.Brush$xmin)) {
        return(need(validate(input$plotTrimmedSample2.Brush), "Please highlight points on the plot."))
      } else if(!input$usePoint) {
        
        # Use brush values
        c(input$plotTrimmedSample2.Brush$xmin, input$plotTrimmedSample2.Brush$xmax) %>%
          round()
        
      } else if(input$usePoint & (is.null(input$plotTrimmedSample2.Click))){
        return(need(validate(input$plotTrimmedSample2.Click), "Please click a position of interest on the plot."))
      } else if(input$usePoint & !input$useIndex){
        
        # Use click values for motifs
        motifId.Click = statisticalAnalysis.Reactive()$sample_alt %>%
          slice(which.closest(index, round(input$plotTrimmedSample2.Click$x))) %>% 
          pull(motif_id)
        
        indices = statisticalAnalysis.Reactive()$sample_alt %>%
          filter(motif_id == motifId.Click) %>%
          .$index %>%
          quantile(., c(0,1), names = F)
        
        c(indices[1] - input$pad, indices[2] + input$pad)

      } else if(input$usePoint & input$useIndex){
        
        index = round(input$plotTrimmedSample2.Click$x)
        c(index - input$pad, index + input$pad)
      }
      }
      
      controlIndices = statisticalAnalysis.Reactive()$sample_df %>%
        filter(index %in% sampleIndices) %>%
        pull(ctrl_index)
      
     return(as.numeric(c(sampleIndices, controlIndices)))
      
    })

    output$print = renderText({
      if(is.null(input$plotTrimmedSample2.Brush) & !input$usePoint){
        
        "HIGHLIGHT SELECTION:\nHighlight positions of interest on the plot, or use point selection"
      } else if((is.null(input$plotTrimmedSample2.Click) | !is.null(input$plotTrimmedSample2.Brush$xmax)) & input$usePoint) {
        if(!input$useIndex)
          {"MOTIF SELECTION:\nSelect a single motif of interest on the plot, or use highlight selection"} else
          {"BASE SELECTION:\nSelect a single base of interest on the plot, or use highlight selection"}
      } else {
        sig_indices = statisticalAnalysis.Reactive()$output_sample_alt %>%
          filter(sig) %>%
          .$index
        
        sequence = statisticalAnalysis.Reactive()$sample_df %>%
          filter(ctrl_index >= sampleIndices.Reactive()[3] & ctrl_index <= sampleIndices.Reactive()[4]) %>%
          dplyr::select(index, ctrl_max_base) %>%
          mutate(ctrl_max_base = as.character(ctrl_max_base)) %>%
          mutate(sig_max_base = {ifelse(index %in% sig_indices, paste0("[",ctrl_max_base,"]"), ctrl_max_base)}) %>%
          .$sig_max_base %>%
          paste0(., collapse = "")
        
        sample_positions = paste(sampleIndices.Reactive()[1:2], sep = " - ", collapse = " - ")
        control_positions = paste(sampleIndices.Reactive()[3:4], sep = " - ", collapse = " - ")
        length = gsub("\\[|\\]", "", sequence) %>% nchar
        
        paste0(
          "REGION SELECTED:\n",
          "Sequence:\t", sequence, "\n",
          "Sample Position:\t", sample_positions, "\n",
          "Control Position:\t", control_positions, "\n",
          "Length:\t", length, " nt"
        )
      }
      })
    
    ## Define a sample chromatogram reactive to an action button
    sampleChromatogram.Plot = eventReactive(input$updateChromatograms, {
      sample_sanger.Reactive() %>%
        # geom_chromotagram:
        # inputs a sanger object and the start and end indices of interest
        # returns a ggplot object consisting of the trace and the percent bases of interest
        geom_chromatogram(., sampleIndices.Reactive()[1], sampleIndices.Reactive()[2])

    })
    
    ## Render plot of reactive sample chromatogram
    output$sampleChromatogram = renderPlot({
      sampleChromatogram.Plot()
    })
    
    ## Define a control chromatogram reactive to an action button
    controlChromatogram.Plot = eventReactive(input$updateChromatograms, {
      if(input$use_fasta) {
        return(need(validate(ctrlChromatogram.Plot(), "Provide a control .ab1 file to generate a chromatogram")))
      } else {
        ctrlChromatogram.Plot = ctrl_sanger.Reactive() %>%
          geom_chromatogram(., sampleIndices.Reactive()[3], sampleIndices.Reactive()[4])
        return(ctrlChromatogram.Plot)
      }
    })
    
    ## Render plot of reactive sample chromatogram
    output$ctrlChromatogram = renderPlot({
      controlChromatogram.Plot()
    })
    
    ## Render plot of the trimmed sample for indexing motifs
    output$plotTrimmedSample2 = renderPlot({
      plotTrimmedSample(
        statisticalAnalysis.Reactive()$sample_df,
        filteredData.Reactive()$pre_cross_align_sample_df,
        statisticalAnalysis.Reactive()$output_sample_alt,
        sample_df.Reactive(),
        statisticalAnalysis.Reactive()$sample_alt
      )
    })
  
    ## Render cox plot
    output$CoxPlot = renderPlot(geom_coxplot(editingData.Reactive()))
    
    ## Render app log
    output$log = renderText(
      # statisticalAnalysis.Reactive()$sample_alt
      paste(
        "PARAMETERS USED:",
        paste0("sample file:\t", sample.Reactive()),
        paste0("control file:\t", ctrl.Reactive()),
        paste0("edit motif:\t", motif.Reactive()),
        paste0("edit(s):\t", wt.Reactive(), ">", edit.Reactive()),
        "",
        "OPERATIONS LOG:",
        logError(ctrl_alignment.Reactive(), "Control data loading"),
        logError(sample_alignment.Reactive(), "Sample data loading"),
        logError(filteredData.Reactive(), "Data trimming on phred score"),
        logError(trimmedData.Reactive(), "Sample data alignment to control data"),
        logError(statisticalAnalysis.Reactive(), "Edit detection and quantification"),
        logError(editingData.Reactive(), "Editing table generation"),
        logError(sampleIndices.Reactive(), "Sample and control indices selection"),
        logError(sampleChromatogram.Plot(), "Sample chromatogram generation"),
        logError(controlChromatogram.Plot(), "Control chromatogram generation"),
      sep = "\n"
      )
      )
    
    
    #######
    # BATCH MODE - added by J Chacón
    ######
    # reads the uploaded parameters excel sheet into a reactive
    # with just basenames for seq files
    multieditR_params <- reactive({
      req(input$batch_parameters)
      dat = readxl::read_excel(input$batch_parameters$datapath)
      dat
    })
    
    # shows the parameters table you uploaded
    output$parameter_table <- renderTable({
      req(input$batch_parameters)
      (multieditR_params())
    })
    
    # makes the sequence file names have the full path (on the server)
    # so that multiedtiR can find them
    params = reactive({
      req(c(input$batch_parameters, input$batch_files))
      if (nrow(missing_files()) == 0){      
        params = add_paths(multieditR_params(), input$batch_files)
        params$sample_file = params$sample_path
        params$ctrl_file = params$ctrl_path
        params
      }else{
        NULL
      }
    })
    
    # the function that actually adds the paths to the sample names
    add_paths = function(params, batch_files){
      params$sample_path = NULL
      params$ctrl_path = NULL
      for (i in 1:nrow(params)){
        params$sample_path[i] = batch_files$datapath[batch_files$name == params$sample_file[i]]
        params$ctrl_path[i] = batch_files$datapath[batch_files$name == params$ctrl_file[i]]
      }
      params
    }
    
    # determines from the parameters sheet what sequence files are needed
    sequence_files_required <- reactive({
      req(c(multieditR_params()))
      #print(input$batch_files$name)
      sample_files = multieditR_params()$sample_file
      control_files = multieditR_params()$ctrl_file
      needed_files = unique(c(sample_files, control_files))
      needed_files = data.frame(`Sequence Files Needed:` = needed_files)
      needed_files
    })
    
    # shows which sequence files are needed
    output$needed_sequence_files <- renderTable({
      sequence_files_required()
    })
    
    # after you upload sequence files, determines what you missed (if any)
    missing_files = reactive({
      req(c(sequence_files_required(), input$batch_files))
      uploaded_files = input$batch_files$name
      missing_files = setdiff(sequence_files_required()[,1], uploaded_files)
      data.frame(`Sequence Files Missing:` = missing_files)
    })
    
    # shows the missing files
    output$missing_files = renderTable({
      req(c(sequence_files_required(), input$batch_files))
      missing_files()
    })
    
    # shows button to run multieditR, but only if there aren't missing seq files
    output$run_button <- renderUI({
      if (!is.null(missing_files()) && nrow(missing_files()) == 0){
        actionButton("run_batch_mode", "Run multieditR")
      }
    })
    
    # runs the analysis and stores the result in the reactive output_data
    output_data <- reactive({
      req(c(input$batch_parameters, input$batch_files))
      if (nrow(missing_files()) == 0){
        analysis()
      }
    })
    
    # the actual calling of multiedtiR
    analysis <- reactive({
      req(c(input$batch_parameters, input$batch_files))
      if (nrow(missing_files()) == 0){
        my_params = params()
        all_results = multieditR::detect_edits_batch(my_params)
        all_results
      }
    })    
    
    # futzes with the results table and shows it
    output$combined_results_table <- DT::renderDataTable({
      req(c(input$batch_parameters, input$batch_files))
      if (nrow(missing_files()) == 0){
        all_results = analysis()
        main_table = lapply(all_results, FUN = function(x){x[[1]]}) %>% do.call(rbind, .)
        main_table = main_table[,c("sample_name","target_base","ctrl_max_base",
                                   "A_perc","A_sig",
                                   "C_perc","C_sig",
                                   "G_perc","G_sig",
                                   "T_perc", "T_sig",
                                   "A_pvalue","C_pvalue","G_pvalue","T_pvalue",
                                   "index", "ctrl_index")]
        main_table$A_perc = signif(main_table$A_perc,3)
        main_table$C_perc = signif(main_table$C_perc,3)
        main_table$T_perc = signif(main_table$T_perc,3)
        main_table$G_perc = signif(main_table$G_perc,3)
        main_table$A_pvalue = signif(main_table$A_pvalue, 3)
        main_table$C_pvalue = signif(main_table$C_pvalue, 3)
        main_table$T_pvalue = signif(main_table$T_pvalue, 3)
        main_table$G_pvalue = signif(main_table$G_pvalue, 3)
        
        #main_table = cbind(main_table[,ncol(main_table)], main_table[, 1:(ncol(main_table)-1)])
        print(main_table, n = Inf)
      }else{
        print("Output Will Appear Here")
      }
    },
    server = FALSE,
    extensions = c('Buttons'), 
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel')
    ))
    
    # shows a render-and-download button, after the analysis finishes
    output$render_and_dl_button <- renderUI({
      if (!is.null(analysis())){
        downloadButton("render_and_dl", "Render and Download HTML report\n(this will take a minute)")
      }
    })
    
    # compiles the HTML report and saves it to your downloads
    output$render_and_dl <- downloadHandler(
      filename = "batch_report.html",
      content = function(file) {
        # to show loading message
        showModal(modalDialog("Preparing Report...(please expect 1-3 minutes depending on the number of samples)", footer=NULL))
        on.exit(removeModal())
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed)
        rmd_loc = paste0(system.file(package = "multieditR"),
                         "/batch_report_template.Rmd")
        temp_dir <- tempdir()
        temp_rmd_loc <- file.path(temp_dir, "batch_report_template.Rmd")
        file.copy(rmd_loc, temp_rmd_loc)
        
        
        my_params = params()
        my_data = output_data()
        rmarkdown::render(temp_rmd_loc, 
                          #runtime = "shiny",
                          output_file = file,
                          params = list(params.tbl = my_params, 
                                        results.list = my_data),
                          envir = new.env(parent = globalenv()))
      }
    )

  }
)
