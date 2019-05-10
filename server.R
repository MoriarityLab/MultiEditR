### server.R for multiEditR

###########################################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

### Example data
example_sample = "RP272_cdna_wt.ab1"
example_control = "RP272_cdna_ko.ab1"
example_control_fasta = "RP272_cdna_ko.fasta"
example_motif = "YAR"
example_wt = "A"
example_edit = "G"
example_use_ctrl_seq = FALSE

### Source helper functions
source("helpers.R")

### Server function
shinyServer(
  function(input, output) {

    # Establish reactive input functions
    sampleReactive = reactive({
      # if else things to handle example data loading
      if(input$use_example) {
        return(example_sample)
      } else if(!is.null(input$sample_file)) {
        return(input$sample_file$datapath)
      } else return(validate(
        need(input$sample_file, "Please upload your sanger sequence sample file")
      ))
    })

    ctrlReactive = reactive({
      # if else things to handle example data loading
      if(input$use_example) {
        return(example_control)
      } else if(!is.null(input$ctrl_file)) {
        return(input$ctrl_file$datapath)
      } else return(validate(
        need(input$ctrl_file, "Please upload your control sequence file")
      ))
    })

    motifReactive = reactive({
      # if else things to handle example data loading
      if(input$use_example) {
        return(example_motif)
      } else if(!is.null(input$motif)) {
        return(input$motif)
      } else return(validate(
        need(input$motif, "Please enter a motif of interest")
      ))
    })

    wtReactive = reactive({
      # if else things to handle example data loading
      if(input$use_example) {
        return(example_wt)
      } else if(!is.null(input$wt)) {
        return(input$wt)
      } else return(validate(
        need(input$wt, "Please enter a wt base of interest")
      ))
    })

    editReactive = reactive({
      # if else things to handle example data loading
      if(input$use_example) {
        return(example_edit)
      } else if(!is.null(input$edit)) {
        return(input$edit)
      } else return(validate(
        need(input$edit, "Please enter an edited base of interest")
      ))
    })
    
    # Generate data for app
    # Run the runEditR function to generate all data
    output_data = reactive({
      input$actionButton
      isolate(
      runEditR(sample_file =  sampleReactive() , # have to use a weird work around of calling the first slot in the list containing the sample name...
               ctrl_file = ctrlReactive(), #input$ctrl_file,
               motif = motifReactive(), #input$motif,
               wt = wtReactive(), #input$wt,
               edit = editReactive(),
               use_ctrl_seq = input$use_fasta,
               phred_cutoff = input$phred,
               p_value = input$p_value) # use double brackets to access slot needed else
    )
    })

    # As soon as this reactive function is called the whole program comes to a stand still once deployed.
    # Could then try truncating the function and seeing where it starts breaking.
    output$table1 = renderTable({
      output_data()[[1]] #1
    })

    output$table2 = renderTable({
      output_data()[[5]] #2
    })

    output$plot1 = renderPlot({
      output_data()[[2]] #3
    })

    output$plot2 = renderPlot({
      output_data()[[3]] #4
    })

    output$plot3 = renderPlot({
      output_data()[[4]] #5
    })

    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$sample_name, ".tsv", sep = "")
      },
      content = function(file) {
        write_tsv(output_data()[[6]], file) #6
      }
    )
  }
)
