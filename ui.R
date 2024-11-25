### ui.R for multiEditR

###########################################################################################
# Copyright (C) 2020-2021 Mitchell Kluesner (klues009@umn.edu)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

# EditR UI
version = "1.1.1"
update = "Feb 28, 2024"

shinyUI(
  navbarPage("multiEditR",
             
             # tabPanel("Single-Sample Mode",
             #          
             #          fluidPage(
             #            theme = shinythemes::shinytheme("cerulean"),
             #            
             #            pageWithSidebar(
             #              
             #              # App title ----
             #              headerPanel(paste0("MultiEditR single-sample mode, v", version)),
             #              
             #              # Sidebar panel for inputs ----
             #              sidebarPanel(
             #                
             #                em(paste0("Please note this application is under development.")),
             #                em(paste0("Updated ", update, ".")),
             #                p(),
             #                tags$a(href="https://www.biorxiv.org/content/biorxiv/early/2019/05/09/633685.full.pdf", "Please cite our bioRχiv preprint!"),
             #                # actionButton(inputId = "actionButton", label = "Update"),
             #                p(),
             #                checkboxInput(inputId = "use_example",
             #                              label = "Click to see an example.", value = FALSE),
             #                fileInput(inputId = 'sample_file',
             #                          label = 'Choose a .ab1 sample file'),
             #                checkboxInput(inputId = "use_fasta",
             #                              label = "Click if control file is in fasta format."),
             #                fileInput(inputId = 'ctrl_file',
             #                          label = 'Choose a .ab1 or .fasta control file'),
             #                # textInput(inputId = 'sample_name',
             #                #           label = 'Enter a sample name'),
             #                textInput(inputId = 'motif',
             #                          label = 'Enter a motif of interest',
             #                          value = ""),
             #                p("Motif is combatible of a sequence of any length and composition, including IUPAC ambiguity codes. (e.g. A, YAR, ACGTYGRGTACA)"),
             #                br(),
             #                textInput(inputId = 'wt',
             #                          label = 'Enter the WT base',
             #                          value = ""),
             #                p("Please enter one base A, C, G, or T under investigation to be edited."),
             #                textInput(inputId = 'edit',
             #                          label = 'Enter the edited base',
             #                          value = ""),
             #                p("Multiple edits of interest (e.g. C > T or G) can be specified by separating bases with | (e.g. T|G)"),
             #                br(),
             #                numericInput(inputId = 'phred',
             #                             label = 'Enter desired phred cutoff',
             #                             value = 0.0001),
             #                numericInput(inputId = 'p_value',
             #                             label = 'Enter P-value cutoff',
             #                             value = 0.0001),
             #                checkboxInput(inputId = "adjust_p",
             #                              label = "Click to apply a multiple comparisons correction."),
             #                p(),
             #                downloadButton("downloadRawData", "Download Raw Data")
             #              ),
             #              
             #              # Main panel for displaying outputs ----
             #              mainPanel(
             #                
             #                tabsetPanel(type = "tabs",
             #                            tabPanel("Instructions",
             #                                     includeMarkdown("instructions.md")
             #                            ),
             #                            
             #                            tabPanel("Analysis",
             #                                     h1("Sample analysis"),
             #                                     br(),
             #                                     h2("Raw sample signal by position"),
             #                                     plotOutput("plotRawSample"),
             #                                     p("This plot shows the signal intensity across all positions in the raw sample file. The greyed regions represent those that have been trimmed by the phred score. The red line represents the average percent signal in the high-quality regions, and the teal regions represent regions containing your motif of interest."),
             #                                     br(),
             #                                     h2("Trimmed sample noise by position"),
             #                                     plotly::plotlyOutput("plotTrimmedSample"),
             #                                     p("This plot shows the noise intensity across all positions in the trimmed sample. The dots represent regions containing the motif of interest. Dots that are colored orange represent noise that measured as non-significant, while teal dots represent noise that measures as significant."),
             #                                     br(),
             #                                     h2("Barplot of editing data"),
             #                                     plotOutput("plotEditingData"),
             #                                     p("This barplot shows the distribution of hypothetical edits that are or are not significant."),
             #                                     br(),
             #                                     h2("Summary table of edits"),
             #                                     tableOutput("tableEditingData"),
             #                                     p("This table shows descriptive statistics of the detected edits."),
             #                                     br(),
             #                                     h2("MultiEditR Editing Index (MEI)"),
             #                                     tableOutput("tableEditingIndexData"),
             #                                     p("This table shows the editing index for each base."),
             #                                     br(),
             #                                     h2("Summary of statistical parameters"),
             #                                     tableOutput("tableZAGAParameters"),
             #                                     p("This table shows how well the zero-adjusted gamma distribution model fits the noise data. A perfect fit is when Filliben's Correlation coefficient is 1. This communicates how reliable the statistical model is.")
             #                                     
             #                            ),
             #                            
             #                            tabPanel("Chromatograms",
             #                                     h1("Chromatograms of sample and control files"),
             #                                     p("Use the bottom interactive plot to choose your specific motif of interest instance. Then change the input motif instance on the left sidepanel."),
             #                                     br(),
             #                                     h2("Edits and their significance"),
             #                                     p("Highlight over points of interest to determine their motif instance number, which can then be entered in the side panel. Points with solid color represent significant edits, while translucent points are not called as significant."),
             #                                     plotOutput("plotTrimmedSample2",
             #                                                click = "plotTrimmedSample2.Click",
             #                                                brush = brushOpts(
             #                                                  id = "plotTrimmedSample2.Brush"
             #                                                )
             #                                     ),
             #                                     verbatimTextOutput('print'),
             #                                     br(),
             #                                     checkboxInput(inputId = "usePoint", label = "Use point selection instead of highlight selection", value = FALSE),
             #                                     conditionalPanel(
             #                                       condition = "input.usePoint==true",
             #                                       checkboxInput(inputId = "useIndex", label = "Use base position selection instead of motif selection", value = FALSE),
             #                                       numericInput(inputId = 'pad',
             #                                                    label = 'Enter desired padding to chromatogram',
             #                                                    value = 5)
             #                                     ),
             #                                     actionButton("updateChromatograms", "Update Chromatograms"),
             #                                     h2("Sample Chromatogram"),
             #                                     p("The motif of interest is centered in the middle of the chromatogram. Additional bases are added based on the desired padding parameter above. Default is 5 bases."),
             #                                     plotOutput("sampleChromatogram"),
             #                                     br(),
             #                                     h2("Control Chromatogram"),
             #                                     p("The motif of interest is centered in the middle of the chromatogram. Additional bases are added based on the desired padding parameter above. Default is 5 bases."),
             #                                     p("If a .fasta file is used as the control this plot will remain blank."),
             #                                     plotOutput("ctrlChromatogram"),
             #                                     br()
             #                            ),
             #                            
             #                            tabPanel("Raw Data", 
             #                                     fluidRow(
             #                                       column(12, DT::dataTableOutput('rawData')
             #                                       )
             #                                     )
             #                            )
             #                            ,
             #                            tabPanel("Logs",
             #                                     h2("Operations Log"),
             #                                     verbatimTextOutput("log"),
             #                                     includeMarkdown("troubleshooting.md")
             #                            )
             #                )
             #              )
             #            )
             #          )
             # ),
             tabPanel("Batch Mode",
                      fluidPage(
                        theme = shinythemes::shinytheme("cerulean"),
                        headerPanel(paste0("MultiEditR batch mode, v1.0")),
                        sidebarLayout(
                          sidebarPanel(
                            em(paste0("Please note this application is under development.")),
                            em(paste0("Updated ", "2024-02-28", ".")),
                            p(),
                            tags$a(href="https://www.biorxiv.org/content/biorxiv/early/2019/05/09/633685.full.pdf", "Please cite our bioRχiv preprint!"),
                            p(),
                            fileInput("batch_parameters", NULL, buttonLabel = "Upload Parameters excel sheet",
                                      accept = c(".xlsx", ".xls"),
                                      label = "Upload a Parameters Excel spreadsheet",
                                      multiple = FALSE),
                            fileInput("batch_files", NULL, buttonLabel = "Upload Sequence Files (.fasta or .ab1)...",
                                      accept = c(".fa",".fasta", ".ab1") ,
                                      label = "Upload All Sequence Files at Once",
                                      multiple = TRUE)
                          ),
                          mainPanel(
                            tabsetPanel(type = "tabs",
                                        tabPanel("Instructions",
                                                 fluidPage(
                                                   tags$iframe(src = "batch_instructions.html",
                                                               allowfullscreen = "true",
                                                               seamless = NA,
                                                               width = 800,
                                                               height = 4800,
                                                               scrolling = "no",
                                                               frameborder = 0)
                                                   #htmltools::includeHTML("www/batch_instructions.html")
                                                   )
                                        ),
                                        tabPanel("Analysis",
                                                 h1("Batch Analysis"),
                                                 em("Upload a parameters table to start seeing output. Click the Instructions tab for details."),
                                                 conditionalPanel("input.parameters",
                                                                  HTML("<b>The parameters table you upload will appear here:</b>")
                                                 ),
                                                 h2("Parameters Table"),
                                                 tableOutput("parameter_table"),
                                                 h2("Sequence files needed based upon Parameters Table"),
                                                 fluidRow(
                                                   column(width = 4, tableOutput("needed_sequence_files")),
                                                   column(width = 4, tableOutput("missing_files"))
                                                 ),
                                                 h2("Click Run to start analysis:"),
                                                 em("Run button appears when all sequence files are uploaded"),
                                                 uiOutput('run_button'),
                                                 conditionalPanel("input.run_batch_mode",
                                                                  em("results and download buttons will appear once analysis is complete..."),
                                                                  HTML("<br>"),
                                                                  h2("Click below to download HTML report"),
                                                                  uiOutput('render_and_dl_button'),
                                                                  HTML("<br>"),
                                                                  h2("Combined Results Table:"),
                                                                  em("note: failed samples do not appear in this table. Check the downloadable report if samples in the parameters sheet do not appear here."),
                                                                  shinycssloaders::withSpinner(
                                                                    DT::dataTableOutput("combined_results_table")
                                                                  ),
                                                 )
                                        ),
                                        tabPanel("FAQ",
                                                 fluidPage(
                                                   htmltools::includeMarkdown(path = "www/FAQ.md")
                                                 )
                                        )
                            )
                          )
                        )
                      )
             )
             
  )
)