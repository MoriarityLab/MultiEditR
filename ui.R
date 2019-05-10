### ui.R for multiEditR

###########################################################################################
# Copyright (C) 2018-2019 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

# EditR UI
version <- "1.0.1"

shinyUI(
  pageWithSidebar(
    
    # App title ----
    headerPanel(paste0("MultiEditR v", version)),
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      actionButton(inputId = "actionButton", label = "Submit"),
      p(),
      checkboxInput(inputId = "use_example",
                    label = "Click to see an example."),
      conditionalPanel(condition = "'input.use_example' == FALSE",
                       fileInput(inputId = 'sample_file',
                                 label = 'Choose a .ab1 sample file'),
                       checkboxInput(inputId = "use_fasta",
                                     label = "Click if control file is in fasta format."),
                       fileInput(inputId = 'ctrl_file',
                                 label = 'Choose a .ab1 or .fasta control file'),
                       textInput(inputId = 'sample_name',
                                 label = 'Enter a sample name'),
                       textInput(inputId = 'motif',
                                 label = 'Enter a motif of interest'),
                       textInput(inputId = 'wt',
                                 label = 'Enter the WT base'),
                       textInput(inputId = 'edit',
                                 label = 'Enter the edited base'),
                       numericInput(inputId = 'phred',
                                    label = 'Enter desired phred cutoff',
                                    value = 0.0001),
                       numericInput(inputId = 'p_value',
                                    label = 'Enter P-value cutoff',
                                    value = 0.001)
      
    ),
    p(),
    downloadButton("downloadData", "Download Data")
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    # Output: Formatted text for caption ----
    # insert the graph you want to display
    p("This plot shows the percent signal of each basecall for each position in the sample chromatogram. Highlighted regions are the motif of interest."),
    p(),
    plotOutput(outputId = "plot1"),
    p(),
    p("This plot shows the percent noise across the sample file. Positions were significant edits detected are indiciated with a colored point."),
    p(),
    plotOutput(outputId = "plot2"),
    p(),
    p("This plot shows bases where editing was or was not detected across the entire sample."),
    p(),
    plotOutput(outputId = "plot3"),
    p(),
    p("This table is a summary of editing in the sample"),
    p(),
    tableOutput(outputId = "table1"),
    p(),
    p("Parameters of the zero-adjusted gamma distribution used to determine significance."),
    p(),
    tableOutput(outputId = "table2")
  )
)
)
