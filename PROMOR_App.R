# Upload library ----
library(shiny)
library(shinyjs)
library(shinydashboard)
library(DT)
library(promor)
library(caret)
library(reshape2)
library(stats)
library(pROC)
library(shinyalert)
library(Biobase)
library(BiocGenerics)
library(limma)
library(pROC)
library(pcaMethods)
library(BiocManager)
options(repos = BiocManager::repositories())

ui <- dashboardPage(skin = "blue",
                    dashboardHeader(title = "PROMOR App", titleWidth = 270),
                    ## Sidebar content
                    dashboardSidebar(width = 270,
                      useShinyjs(),
                      sidebarMenu(
                        id = "tabs",
                        menuItem("About", tabName = "home",icon = icon("home")),
                        tags$hr(),
                        menuItem("User Manual", tabName = "manual",icon = icon("book")),
                        tags$hr(),
                        menuItem("Upload Files", expandedName = "upload", icon = icon("file"),
                                 fileInput("file1", "Upload a proteinGroups.txt file produced by MaxQuant software",
                                           accept = c("text/csv", "text/comma-separated-values,text/plain",".csv",
                                                      options(shiny.maxRequestSize=90000*1024^2))),
                                 fileInput("file2", "Upload a expDesign.txt file which describes the experimental design",
                                           accept = c("text/csv", "text/comma-separated-values,text/plain",".csv",
                                                      options(shiny.maxRequestSize=90000*1024^2)))),
                        
                        menuItem("Data Pre-processing:",tabName="upload", icon = icon("cogs"),
                                 checkboxInput("dp1","Filter out empty rows and columns from the data??", TRUE),
                                 checkboxInput("dp2","Filter out reverse proteins??", TRUE),
                                 numericInput("dp3","The minimum number of unique peptides required to identify a protein:",min = 1, max = 20, value = 2),
                                 checkboxInput("dp4","Do you have technical replicates??", FALSE),
                                 checkboxInput("dp5","Repalce zeros with NAs.??", TRUE),
                                 checkboxInput("dp6","Log transformed your data??", TRUE),
                                 numericInput("dp7","Logarithm base",min = 0, max = 10, value = 2),
                                 checkboxInput("dp8","Get Gene Names in the result table.", FALSE),
                                 helpText("*Only select this option, if you have"),
                                 helpText("Gene Names column in your data."),
                                 helpText("*If you select this option, it will only"), 
                                 helpText("provide Gene Names with other results."),
                                 helpText("**Default result table provides"),
                                 helpText("Protein IDs with other results.")
                                 ),
                        
                        menuItem("Data Quality Control:",tabName="upload", icon = icon("flask"),
                                 numericInput("dq1","The proportion of missing data allowed (%):",min = 0.1, max = 2, value = 0.40),
                                 radioButtons("dq2","Filter condition for missing data:", choices = list("either", "each"), inline = TRUE, selected = "either"),
                                 conditionalPanel(condition = "input.dp4 == true",
                                                  radioButtons("dq3","Do you want to remove any sample from data??", choices = list("No", "Yes"), inline = TRUE, selected = "No"),
                                                  conditionalPanel(condition = "input.dq3 == 'Yes'",
                                                  textInput("dq4","Enter Name of the sample to remove:"))),
                                 radioButtons("dq5","Select method for Impute missing data:", choices = list("minProb", "minDet", "RF", "kNN", "SVD"), selected = "minProb"),
                                              conditionalPanel(condition = "input.dq5 == 'minProb'",
                                                numericInput("dq6","tune_sigma:",min = 0.1, max = 99, value = 1),
                                                numericInput("dq7","q:",min = 0.01, max = 10, value = 0.01),
                                                numericInput("dq8","seed:",min = 1, max = 999999, value = 3312),
                                                helpText("Enter any random number in the"),
                                                helpText("above seed option. i.e. 3312"),
                                                helpText("**DO NOT LEAVE BLANK the seed option.")),
                                              conditionalPanel(condition = "input.dq5 == 'minDet'",
                                                numericInput("dq9","q:",min = 0.01, max = 10, value = 0.01),
                                                numericInput("dq10","seed:",min = 1, max = 999999, value = 3312),
                                                helpText("Enter any random number in the"),
                                                helpText("above seed option. i.e. 3312"),
                                                helpText("**DO NOT LEAVE BLANK the seed option.")),
                                              conditionalPanel(condition = "input.dq5 == 'RF'",
                                                numericInput("dq11","maxiter:",min = 1, max = 99999, value = 10),
                                                numericInput("dq12","ntree:",min = 1, max = 99999, value = 20),
                                                numericInput("dq13","seed:",min = 1, max = 999999, value = 3312),
                                                helpText("Enter any random number in the"),
                                                helpText("above seed option. i.e. 3312"),
                                                helpText("**DO NOT LEAVE BLANK the seed option.")),
                                              conditionalPanel(condition = "input.dq5 == 'kNN'",
                                                numericInput("dq14","seed:",min = 1, max = 999999, value = 3312),
                                                helpText("Enter any random number in the"),
                                                helpText("above seed option. i.e. 3312"),
                                                helpText("**DO NOT LEAVE BLANK the seed option.")),
                                              conditionalPanel(condition = "input.dq5 == 'SVD'",
                                                numericInput("dq15","n_pcs:",min = 1, max = 99999, value = 2),
                                                numericInput("dq16","seed:",min = 1, max = 999999, value = 3312),
                                                helpText("Enter any random number in the"),
                                                helpText("above seed option. i.e. 3312"),
                                                helpText("**DO NOT LEAVE BLANK the seed option.")),
                                 radioButtons("dq17","Select the normalization method to use:", choices = list("quantile", "none", "scale", "cyclicloess"), selected = "quantile")),
                      
                      menuItem("Differentially expressed proteins:",tabName="upload", icon = icon("tasks"),
                               numericInput("dex1","Cutoff value for p-values and adjusted p-values:",min = 0.01, max = 999, value = 0.05),
                               numericInput("dex2","Minimum absolute log2-fold change to use as threshold for differential expression (lfc):",min = 1, max = 9999, value = 1),
                               radioButtons("dex3","Select  adjustment method:", choices = list("none", "bonferroni", "holm", "BH", "fdr", "BY"), selected = "BH"),
                               checkboxInput("dex4","Limma Results", FALSE),
                               checkboxInput("dex5","Top Hits Results", FALSE),
                               conditionalPanel(condition = "input.dex5 == true",
                                numericInput("dex6","Enter Top Hits number:",min = 1, max = 99999, value = 10),
                                helpText("Default is 10 but you can change it."))
                               ),
                      
                      menuItem("Protein Modeling:",tabName="upload", icon = icon("th"),
                               radioButtons("pm1","Do you want to build predictive models ??", choices = list("No", "Yes"), selected = "No"),
                               conditionalPanel(condition = "input.pm1 == 'Yes'",
                                h4("Select parameters for pre process:"),
                                radioButtons("pmp1","Criteria to denote significance:", choices = list("adjP", "P"), inline = TRUE, selected = "adjP"),
                                helpText("*adjP for adjusted p-value and *P for p-value."),
                                numericInput("pmp2","Cutoff value for p-values and adjusted p-values in differential expression:",min = 0.01, max = 999, value = 0.05),
                                numericInput("pmp3","Minimum absolute log-fold change to use as threshold for differential expression (lfc):",min = 1, max = 9999, value = 1),
                                numericInput("pmp4","Enter Top Hits number:",min = 1, max = 99999, value = 10),
                                checkboxInput("pmp5","Find highly correlated proteins:", TRUE),
                                numericInput("pmp6","Enter a numeric value specifying the correlation cutoff:",min = 0.01, max = 9999, value = 0.90),
                                checkboxInput("pmp7","Remove highly correlated proteins (predictors or features):", TRUE),
                                checkboxInput("pmp8","Save protein correlation matrix:", TRUE),
                                radioButtons("pmp9","Select  adjustment method:", choices = list("none", "bonferroni", "holm", "BH", "fdr", "BY"), selected = "BH"),
                                h4("Select parameters for split data:"),
                                numericInput("pms1","Enter The size of the training data set:",min = 0.01, max = 99999, value = 0.80),
                                numericInput("pms2","seed:",min = 1, max = 999999, value = 8314),
                                helpText("Enter any random number in the"),
                                helpText("above seed option. i.e. 8314"),
                                helpText("**DO NOT LEAVE BLANK the seed option."),
                                h4("Select parameters for train models:"),
                                selectInput("pmpmt1", "Select resampling method to use:",
                                                         c("None" = "none",
                                                           "boot" = "boot",
                                                           "boot632" = "boot632",
                                                           "optimism_boot" = "optimism_boot",
                                                           "boot_all" = "boot_all",
                                                           "cv" = "cv",
                                                           "repeatedcv" = "repeatedcv",
                                                           "LOOCV" = "LOOCV",
                                                           "LGOCV" = "LGOCV",
                                                           "oob" = "oob",
                                                           "timeslice" = "timeslice",
                                                           "adaptive_cv" = "adaptive_cv",
                                                           "adaptive_boot" = "adaptive_boot",
                                                           "adaptive_LGOCV" = "adaptive_LGOCV"), selected = "repeatedcv"),
                                conditionalPanel(condition = "input.pmpmt1 == 'repeatedcv'",
                                                 numericInput("pmpmt2","The number of complete sets of folds to compute:",min = 1, max = 999999, value = 3)),
                                numericInput("pmpmt3","Number of resampling iterations:",min = 1, max = 999999, value = 10),
                                numericInput("pmpmt4","Random number seed:",min = 1, max = 999999, value = 351),
                                selectInput("pmpmt5", "Select machine learning algorithms to use:",
                                           c("rf" = "rf",
                                             "svmRadial" = "svmRadial",
                                             "glm" = "glm",
                                             "xgbLinear" = "xgbLinear"), multiple = TRUE),
                                helpText("**Select minimum two algorithms above."),
                                helpText("**DO NOT LEAVE BLANK the above option."),
                                h4("Select parameters for test models:"),
                                radioButtons("pmtm1","Type of output:", choices = list("prob", "raw"), inline = TRUE, selected = "prob"),
                                conditionalPanel(condition = "input.pmtm1 == 'raw'",
                                  checkboxInput("pmtm2","Save Confusion matrices:", TRUE))
                                )),
                      actionButton('subm', 'Submit'),
                      tags$hr(),
                      hidden(menuItem("results", tabName = "results"))
                      )
                    ),
                    
                    ## Body content
                    dashboardBody(
                      useShinyjs(),
                      tags$head(includeScript("google-analytics.js")),
                      tabItems(
                        # First tab content
                        tabItem(tabName = "home",
                                fluidRow(
                                  box(
                                    title = "About this App",
                                    h3("PROMOR App Ver 0.1.0"),
                                    tags$ul(tags$li(p("PROMOR App is an interactive web application to analyze and visualize label-free quantification (LFQ) proteomics data preprocessed using MaxQuant software.")),
                                            tags$li(p("PROMOR App also provides an option to build predictive models based on machine learning-based modeling.")),
                                            tags$li(p("PROMOR (PROtein MOdeling using R programming language) App is based on ",a(href="https://caranathunge.github.io/promor/index.html", target="_blank", "promor R package."))),
                                            tags$li(p("PROMOR App source code is availabe at GitHub:", a(href="https://github.com/sgr308/PROMOR_App", target="_blank", "https://github.com/sgr308/PROMOR_App"))),
                                            tags$li(p("Users can use all options and parameters described in promor R package and do analysis using PROMOR App.")),
                                            tags$li(p("Please click on", tags$b("User Manual"), "on the left sidebar for more information about preparing files and how to run PROMOR App."))),
                                    h3("Please Cite following publications/URLs if you have generated any plot or table results in your study using PROMOR App."),
                                    tags$ul(
                                      tags$li(p(tags$b("https://sgrbnf.shinyapps.io/PROMOR_App/"))),
                                      tags$li(p(tags$b("https://github.com/sgr308/PROMOR_App"))),
                                      tags$li(p(tags$b("Ranathunge, C., S. S. Patel, L. Pinky, V. L. Correll, O. J. Semmes, R. K. Armstrong, C. D. Combs, and J. O. Nyalwidhe, 2022: promor: a comprehensive r package for label-free proteomics data analysis and predictive modeling. bioRxiv, doi:10.1101/2022.08.17.503867"))),
                                      ),
                                    h3("Contact Us"),
                                    p("Please contact us if you have any questions or feedback related to PROMOR App or promor R package."),
                                    tags$ul(
                                      tags$li(p(tags$b("Sagar Patel, Ph.D."), "Email: patelss(at)evms.edu, sgr.bnf(at)gmail.com", a(href="https://www.linkedin.com/in/sgr308", target="_blank", "LinkedIn"), a(href="https://twitter.com/sgr308", target="_blank", "Twitter"), a(href="https://github.com/sgr308", target="_blank", "GitHub"), a(href="https://orcid.org/0000-0003-4896-2658", target="_blank", "ORCID"),"(For PROMOR App)")),
                                      tags$li(p(tags$b("Chathurani Ranathunge, Ph.D."), "Email: ranathca(at)evms.edu, caranathunge86(at)gmail.com", a(href="https://www.linkedin.com/in/chathurani-ranathunge", target="_blank", "LinkedIn"), a(href="https://twitter.com/caranathunge", target="_blank", "Twitter"), a(href="https://github.com/caranathunge", target="_blank", "GitHub"), "(For promor R package)"))
                                    ),
                                    h3("News and Updates"),
                                    tags$ul(
                                      tags$li("6-October-2022: PROMOR App released.")
                                    ),
                                    h3("This App is Supported by:"),
                                    br(),
                                    tags$a(href="https://www.evms.edu/", target='_blank', tags$img(src="./Eastern_Virginia_Medical_School_logo.svg", width="400")),
                                    tags$a(href="https://hrbrc.org/", target='_blank', tags$img(src="./HRBRC.jpg", width="400")),
                                    width = 12,
                                    solidHeader = TRUE,
                                    status = "success"
                                  ))
                        ),
                        # Second tab content
                        tabItem(tabName = "results",
                                conditionalPanel(condition = "input.dp4 == true",
                                  fluidRow(
                                  box(
                                    title = "Correlation plots: technical replicates",
                                    fluidRow(
                                      box(numericInput("crp1","Select First Replicate to plot:",min = 1, max = 20, value = 1),
                                      numericInput("crp2","Select Second Replicate to plot:",min = 1, max = 20, value = 3), width = 6),
                                      box(numericInput("crp3","Text Size in the plot:",min = 1, max = 100, value = 10),
                                      numericInput("crp4","Plot resolution (dpi):",min = 20, max = 500, value = 80),width = 6)),
                                    radioButtons("crp5","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                    plotOutput("crp6", height = 600),
                                    radioButtons("crp7","Select the file type for plots:", choices = list("png", "pdf"), inline = TRUE),
                                    downloadButton("crp8","Download Correlation plots"),
                                    width = 6,
                                    solidHeader = TRUE,
                                    status = "info"
                                  )
                                  )),
                                
                                fluidRow(
                                  column(width=6,
                                    tabBox(title = "Data Quality Results", width = 12,
                                           tabPanel(title = "Missing data - Heatmap",
                                                    fluidRow(
                                                      box(numericInput("dqh1","The range of first protein to plot:",min = 1, max = 999999, value = 1),
                                                        numericInput("dqh2","The range of last protein to plot:",min = 1, max = 999999, value = 10),
                                                        helpText("The plot shows protein range from 1:10. If you want to see all proteins in the plot then replace value 10 with the following total number of proteins in the above option:"),
                                                        tags$b((textOutput("prange"))), width = 6),
                                                      box(numericInput("dqh3","The range of first sample to plot:",min = 1, max = 999999, value = 1),
                                                        numericInput("dqh4","The range of last sample to plot:",min = 1, max = 999999, value = 6),
                                                        helpText("The plot shows sample range from 1:6. If you want to see all samples in the plot then replace value 6 with the following total number of samples in the above option:"),
                                                        tags$b((textOutput("srange"))), width = 6)),
                                                    fluidRow(
                                                      box(checkboxInput("dqh5","Reorder samples on the x axis:", FALSE),
                                                        checkboxInput("dqh6","Reorder proteins on the y axis:", FALSE), width = 4),
                                                      box(radioButtons("dqh7","Function to reorder samples along the x axis:", choices = list("mean", "sum"), inline = TRUE, selected = "mean"),
                                                          helpText("This option only works if you select *Reorder samples on the x axis* option on the left side."),
                                                          radioButtons("dqh8","Function to reorder proteins along the y axis:", choices = list("mean", "sum"), inline = TRUE, selected = "mean"),
                                                          helpText("This option only works if you select *reorder proteins along the y axis* on the left side."), width = 4),
                                                      
                                                      box(checkboxInput("dqh9","Label Protein IDs:", FALSE),
                                                        numericInput("dqh10","Text Size in the plot:",min = 1, max = 100, value = 10),
                                                        numericInput("dqh11","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 4)),
                                                    radioButtons("dqh12","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                                    plotOutput("dqh13", height=600),
                                                    radioButtons("dqh14", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                                    downloadButton("dqh15","Download Heat map plot")
                                                    ),
                                           tabPanel(title = "Impute missing data",
                                                    fluidRow(
                                                      box(checkboxInput("dqim1","Global density plot:", TRUE),
                                                      conditionalPanel(condition = "input.dqim1 == false",
                                                        numericInput("dqim2","Number of rows to print sample-wise density plot:",min = 1, max = 9999, value = 3),
                                                        numericInput("dqim3","Number of columns to print sample-wise density plot:",min = 1, max = 9999, value = 2)), width = 6),
                                                      box(numericInput("dqim4","Text Size in the plot:",min = 1, max = 100, value = 10),
                                                        numericInput("dqim5","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 6)),
                                                    radioButtons("dqim6","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                                    plotOutput("dqim7", height=600),
                                                    radioButtons("dqim8", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                                    downloadButton("dqim9","Download Impute missing plot")
                                                    ),
                                           tabPanel(title = "Data Normalization",
                                                    fluidRow(
                                                      box(radioButtons("dqn1", "Type of plot to generate:", choices = list("box", "density"), inline = TRUE), width = 4),
                                                      box(numericInput("dqn2","Text Size in the plot:",min = 1, max = 100, value = 10), width = 4),
                                                      box(numericInput("dqn3","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 4)),
                                                    radioButtons("dqn4","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                                    plotOutput("dqn5", height=600),
                                                    radioButtons("dqn6", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                                    downloadButton("dqn7","Download normalization plot")
                                                    )
                                           )),
                                           
                                  box(
                                    title = "Results Table",
                                    infoBoxOutput("boxres",width = 12),
                                    DT::dataTableOutput("contents"),
                                    downloadButton("downtab","Download all data as .csv file"),
                                    br(),
                                    br(),
                                    conditionalPanel(condition = "input.dex4 == true",
                                      downloadButton("downtab2","Download Limma .csv file")),
                                    br(),
                                    conditionalPanel(condition = "input.dex5 == true",
                                      downloadButton("downtab3","Download Top Hits .csv file")),
                                    width = 6,
                                    solidHeader = TRUE,
                                    status = "danger"
                                  )
                                ),
                                
                                fluidRow(
                                  column(width=6,
                                         tabBox(title = "Results Plots", width = 12,
                                                tabPanel(title = "Volcano Plot",
                                                         fluidRow(
                                                           box(radioButtons("rvp1","Criteria to denote significance:", choices = list("adjP", "P"), inline = TRUE, selected = "adjP"),
                                                               helpText("*adjP for adjusted p-value and *P for p-value."),
                                                               checkboxInput("rvp2","See a dotted line to indicate the lfc threshold in the plot.", TRUE),
                                                               checkboxInput("rvp3","See a dotted line to indicate the p-value cutoff.", TRUE), width = 6),
                                                           box(checkboxInput("rvp4","Label Protein IDs:", TRUE),
                                                               numericInput("rvp5","Text Size in the plot:",min = 1, max = 100, value = 10),
                                                               numericInput("rvp6","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 6)),
                                                         radioButtons("rvp7","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                                         plotOutput("rvp8", height=600),
                                                         radioButtons("rvp9", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                                         downloadButton("rvp10","Download Volcano plot")
                                                ),
                                                tabPanel(title = "Heat Map DE Plot",
                                                         fluidRow(
                                                           box(radioButtons("rhp1","Criteria to denote significance:", choices = list("adjP", "P"), inline = TRUE, selected = "adjP"),
                                                               helpText("*adjP for adjusted p-value and *P for p-value."),
                                                               numericInput("rhp2","Text Size in the plot:",min = 1, max = 100, value = 10),
                                                               numericInput("rhp3","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 12)),
                                                         radioButtons("rhp4","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                                         plotOutput("rhp5", height=600),
                                                         radioButtons("rhp6", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                                         downloadButton("rhp7","Download Heat map DE plot")
                                                )
                                        )),
                                conditionalPanel(condition = "input.pm1 == 'Yes'",
                                  fluidRow(
                                    column(width=6,
                                      tabBox(title = "Protein data modeling Results", width = 12,
                                        tabPanel(title = "Feature variation plots",
                                          fluidRow(
                                            box(radioButtons("pmf1","Plot type:", choices = list("box", "density"), inline = TRUE, selected = "box"),
                                                numericInput("pmf2","Number of rows to print the plots:",min = 1, max = 100, value = 5),
                                                numericInput("pmf3","Number of columns to print the plots:",min = 1, max = 100, value = 2),
                                                helpText("**Number of Rows and columns must be equal or greater than Top hits option,"),
                                                helpText("if Top hits is 10, then you can try any one of the following options:"),
                                                helpText("number of rows:5 and number of columns:2"),
                                                helpText("number of rows:2 and number of columns:5"), width = 6),
                                            box(numericInput("pmf4","Text Size in the plot:",min = 1, max = 100, value = 10),
                                                numericInput("pmf5","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 6)),
                                            radioButtons("pmf6","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                            plotOutput("pmf7", height=600),
                                            radioButtons("pmf8", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                            downloadButton("pmf9","Download Feature plot"),
                                            br(),
                                            br(),
                                            conditionalPanel(condition = "input.pmp8 == true",
                                              downloadButton("downtab4", "Download Protein correlation matrix"))
                                            ),
                                        tabPanel(title = "Performance Plot",
                                          fluidRow(
                                            box(radioButtons("pmpr1","Type of plot to generate:", choices = list("box", "dot"), inline = TRUE, selected = "box"),
                                              numericInput("pmpr2","Text Size in the plot:",min = 1, max = 100, value = 10),
                                              numericInput("pmpr3","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 12)),
                                            radioButtons("pmpr4","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                            plotOutput("pmpr5", height=600),
                                            radioButtons("pmpr6", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                            downloadButton("pmpr7","Download Performance plot")
                                            ),
                                        tabPanel(title = "Variable importance plot",
                                          fluidRow(
                                            box(radioButtons("pmvp1","Plot type:", choices = list("lollipop", "bar"), inline = TRUE, selected = "lollipop"),
                                              numericInput("pmvp2","Number of rows to print the plots:",min = 1, max = 100, value = 4),
                                              numericInput("pmvp3","Number of columns to print the plots:",min = 1, max = 100, value = 2), width = 6),
                                            box(numericInput("pmvp4","Text Size in the plot:",min = 1, max = 100, value = 10),
                                              numericInput("pmvp5","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 6)),
                                              radioButtons("pmvp6","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                            plotOutput("pmvp7", height=600),
                                            radioButtons("pmvp8", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                            downloadButton("pmvp9","Download Variable importance plot")
                                            ),
                                        tabPanel(title = "ROC plot",
                                          conditionalPanel(condition = "input.pmtm1 == 'prob'",
                                          fluidRow(
                                          box(checkboxInput("pmroc1","Multiple plots:", TRUE),
                                              numericInput("pmroc2","Text Size in the plot:",min = 1, max = 100, value = 10),
                                              numericInput("pmroc3","Plot resolution (dpi):",min = 20, max = 500, value = 80), width = 6)),
                                            radioButtons("pmroc4","Select Color Palette for plots:", choices = list("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"), inline = TRUE),
                                            plotOutput("pmroc5", height=600),
                                            radioButtons("pmroc6", "Select the file type", choices = list("png", "pdf"), inline = TRUE),
                                            downloadButton("pmroc7","Download ROC plot")),
                                          conditionalPanel(condition = "input.pmtm1 == 'raw' && input.pmtm2 == true",
                                            downloadButton("downtab5", "Download Confusion matrices"))
                                            ) 
                                      )))),                 
                                    )
                                  ),
                        #Third tab item
                        tabItem(tabName = "manual",
                                fluidRow( 
                                  box(
                                    title = "User Manual of PROMOR App",
                                    h3("PROMOR App Manual"),
                                    tags$ul(
                                      tags$li(p("PROMOR App is based on ",a(href="https://caranathunge.github.io/promor/index.html", target="_blank", "promor R package."))),
                                      tags$li(p("Users can use all options and parameters described in promor R package and do analysis using PROMOR App.")),
                                      tags$li(p("PROMOR App is developed to perform differential expression analysis of label-free quantification (LFQ) proteomics data and build predictive models based on machine learning-based modeling with top protein candidates.")),
                                      tags$li(p("PROMOR App provides a range of quality control and visualization tools at the protein level to analyze label-free proteomics data.")),
                                      tags$li(p("PROMOR App requires two Input files: one is", tags$b("proteinGroups.txt"), "file produced by MaxQuant and second is an", tags$b("expDesign.txt"), "which contains the experimental design of your proteomics data.")),
                                    ),
                                    h3("Input Data"),
                                    tags$ol(
                                      tags$li(p(tags$b("proteinGroups.txt:"), "This is one of the output files generated by MaxQuant program. It is a tab-delimited file that contains information on identified proteins from your peptide data.")),
                                      tags$li(p(tags$b("expDesign.txt:"), "This is a tab-delimited text file that contains the design of your experiment. Note that you will have to create and provide this file when you run PROMOR App with your own data.")),
                                    ),
                                    tags$ul(
                                      tags$li(p("Please click on the following files to prepare your input data as shown in the same format in the files.
                                                After clicking on the following files, you can copy all data and paste in Notepad and save as a new file.
                                                Then you can upload these two files in this PROMOR APP for data analysis as described in YouTube tutorial.
                                                These data are from previously published data set of label-free quantification (LFQ) proteomics data that do not contain technical replicates from",a(href="https://europepmc.org/article/MED/24942700#id609082", target="_blank", "Cox et al. (2014)"))),
                                      tags$li(p(a(href="https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt", target="_blank", "proteinGroups.txt"))),
                                      tags$li(p(a(href="https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt", target="_blank", "expDesign.txt")))
                                      ),
                                    h3("YouTube tutorials for PROMOR App:"),
                                    tags$ol(
                                      tags$li(p(a(href="https://youtu.be/DWQeW74Lluo", target="_blank", "YouTube tutorial of PROMOR App without technical replicates data"))),
                                      tags$li(p(a(href="https://youtu.be/iMXZFCmadc8", target="_blank", "YouTube tutorial of PROMOR App with technical replicates data"))),
                                      tags$li(p(a(href="https://youtu.be/m15gL-0vwC4", target="_blank", "YouTube tutorial of PROMOR App with Modeling data"))),
                                      ),
                                    h3("Please perfom the following steps on the left sidebar to use this App:"),
                                    tags$ol(tags$li(tags$b("Upload Files:"),"Upload a", tags$b("proteinGroups.txt"), "file of Proteomics data produced by MaxQuant and Upload a", tags$b("expDesign.txt"), "file which describes the experiment design."),
                                            tags$li(tags$b("Data Pre-processing:"),"Select your options to pre-process the proteomics data and select if you have replicate data."),
                                            tags$li(tags$b("Data Quality Control:"),"Select options for missing data and select filter condition."),
                                            tags$li(tags$b("Differentially expressed proteins:"),"Select options to find differentially expressed proteins."),
                                            tags$li(tags$b("Protein Modeling:"),"After you select", tags$b("Yes,"), "you can select options for Protein Modeling."),
                                            tags$li(tags$b("Submit:"),"Once you select all options for your data, then click on", tags$b("Select"), "button to get results. Wait for sometime as some algorithm will take time to compute and you will see your results.")),
                                    h3("Results:"),
                                    tags$ul(tags$li(p("Results page will show you plots and tables based on your input data which are:"))),
                                    tags$ol(tags$li(tags$b("Correlation plots:"),"You will only see this plot if you had selected technical replicates in the *Data Pre-processing option."),
                                            tags$li(tags$b("Data Quality Results:"),"You will get results for Data Pre-processing, which are:", tags$b("Missing data - Heatmap,"), tags$b("Impute missing data,"), "and", tags$b("Data Normalization.")),
                                            tags$li(tags$b("Results Table:"),"The table gives information about differentially expressed proteins."),
                                            tags$li(tags$b("Results Plots:"),"You will see", tags$b("Volcano Plot"),"and", tags$b("Heat Map DE Plot."), "of differentially expressed proteins."),
                                            tags$li(tags$b("Protein data modeling Results:"),"You will only see this plots if you had selected", tags$b("Yes,"), " in the *Do you want to build predictive models option."),
                                    ),
                                    tags$ul(tags$li(p("Results in plots can be downloaded in", tags$b(".pdf"), "and", tags$b(".png"), "format. Also, user can change text size, plot color and other options.")),
                                            tags$li(p("Results in Table can be downloaded in", tags$b(".csv"), "format. Also, user can change several options and save the table."))
                                    ),
                                    width = 12,
                                    solidHeader = TRUE,
                                    status = "warning"
                                  ) 
                                )
                        )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$subm, {
    if (is.null(input$file1)) return(NULL)
    if (is.null(input$file2)) return(NULL)
    updateTabItems(session, "tabs", selected = "results" )
    
    ## Shinyalert
    observeEvent(input$subm ,{ 
      if(input$subm==0 ){
        return()
      }
      
      shinyalert("Please wait patiently :)", "Your Data analysis has started, wait until table and plots
                appear on the screen.", type="info",
                 closeOnClickOutside = TRUE,
                 closeOnEsc = TRUE,
                 timer = 8000) # timer in milliseconds (8 sec)
    })
   
    create_df2 <- function(prot_groups = input$file1$datapath,
                          exp_design = input$file2$datapath,
                          filter_na = input$dp1,
                          filter_prot = input$dp2,
                          uniq_pep = input$dp3,
                          tech_reps = input$dp4,
                          zero_na = input$dp5,
                          log_tr = input$dp7,
                          base = input$dp7) {
      # Load the data
      df <- read.csv(prot_groups,
                     sep = "\t",
                     stringsAsFactors = FALSE
      )
      
      # Load the design file, which is a tab-delimited file containing the
      # experimental design.
      design <- read.csv(exp_design,
                         stringsAsFactors = FALSE,
                         sep = "\t"
      )
      
      if (tech_reps == TRUE) {
        # if tech_reps == TRUE, combine all columns to make new sample label
        design$new_label <- paste(design$condition,
                                  design$sample_ID,
                                  design$tech_rep,
                                  sep = "_"
        )
        
        # If tech_reps == FALSE, remove the tech_rep column from the design file and
        # combine remaining columns to make a new sample label
      } else {
        design$tech_rep <- NULL
        design$new_label <- paste(design$condition,
                                  design$sample_ID,
                                  sep = "_"
        )
      }
      # extract number of rows
      orig_rows <- nrow(df)
      # extract number of columns
      orig_col <- ncol(df)
      
      # Filter out empty rows and columns if they exist in the dataframe.
      if (filter_na == TRUE) {
        # Remove proteins (rows) with missing values (NA) across all samples
        df <- df[rowSums(is.na(df)) != ncol(df), ]
        # Remove samples (columns) with missing values (NA) across all proteins
        df <- df[, colSums(is.na(df)) != nrow(df)]
        # calculate number of rows and columns removed
        rem_row <- orig_rows - nrow(df)
        rem_col <- orig_col - ncol(df)
        message(paste0(rem_row, " empty row(s) removed."))
        message(paste0(rem_col, " empty column(s) removed."))
      } else {
        warning("Data frame may contain empty rows and/or columns.")
      }
      # Filter out some proteins based on specific columns. First check if the columns
      # are present in the data frame before removing rows based on the presence of
      # "+" signs.
      # Get number of rows in the df
      orig_rows_1 <- nrow(df)
      
      if (filter_prot == TRUE) {
        if ("Only.identified.by.site" %in% colnames(df)) {
          df <- subset(
            df,
            df$Only.identified.by.site != "+"
          )
          message(paste0(
            orig_rows_1 - nrow(df),
            " protein(s) (rows) only identified by site removed."
          ))
        }
        
        # get number of rows
        orig_rows_2 <- nrow(df)
        
        if ("Reverse" %in% colnames(df)) {
          df <- subset(
            df,
            df$Reverse != "+"
          )
          message(paste0(
            orig_rows_2 - nrow(df),
            " reverse protein(s) (rows) removed."
          ))
        }
        
        # get number of rows
        orig_rows_3 <- nrow(df)
        
        if ("Potential.contaminant" %in% colnames(df)) {
          df <- subset(
            df,
            df$Potential.contaminant != "+"
          )
        }
        
        if ("Contaminant" %in% colnames(df)) {
          df <- subset(
            df,
            df$Contaminant != "+"
          )
          message(paste0(
            orig_rows_3 - nrow(df),
            " protein contaminant(s) (rows) removed."
          ))
        }
        
        # get number of rows
        orig_rows_4 <- nrow(df)
        
        if ("Unique.peptides" %in% colnames(df)) {
          df <- subset(
            df,
            df$Unique.peptides > uniq_pep
          )
          message(paste0(
            orig_rows_4 - nrow(df), " protein(s) identified by ",
            uniq_pep, " or fewer unique peptides removed."
          ))
        }
      } else {
        warning("Proteins have not been filtered")
      }
      
      # Extract majority protein group names
      maj_proteins <- df$Majority.protein.IDs
      
      # Extract gene ID
      gene_names1 <- df$Gene.names
      
      # Subset the data frame to only include the LFQ.intensity columns
      samples <- df[, grepl("LFQ.intensity", names(df))]
      
      # order dataframe columns by column name. Important for next steps involving
      # mapply.
      samples <- samples[, order(names(samples))]
      df <- as.matrix(samples)
      
      # remove LFQ intensity part from the column name
      raw_col <- gsub("LFQ.intensity.", "", colnames(df))
      
      # sort the design table by mq_label so that the order matches to that of
      # raw_col
      design <- design[order(as.character(design$mq_label)), ]
      
      # Compare the mq_label column in the design file with raw_col and replace
      # raw_col with the appropriate new_label.
      raw_col_edited <- mapply(gsub,
                               design$mq_label,
                               design$new_label,
                               raw_col,
                               USE.NAMES = FALSE
      )
      
      # add the newly edited column names to the matrix
      colnames(df) <- raw_col_edited
      
      if (input$dp8 == "TRUE"){
        # add gene names tot the matrix
        rownames(df) <- c(gene_names1)
      }
      else
      # add majority protein names to the matrix
      rownames(df) <- c(maj_proteins)
      
      # Convert zeros to NA
      if (zero_na == TRUE) {
        # Convert matrix to dataframe and convert zeros to NAs
        df <- as.data.frame(df)
        df[df == 0] <- NA
      } else {
        warning("Zeros have not been converted to NAs in the data frame")
      }
      
      # log2 transform the data
      if (log_tr == TRUE) {
        df <- log(df, base)
      } else {
        warning("Intensities have not been log transformed")
      }
      return(df)
      }
    
    raw <- create_df2(input$file1$datapath, input$file2$datapath,
                     filter_na = input$dp1,
                     filter_prot = input$dp2,
                     uniq_pep = input$dp3,
                     tech_reps = input$dp4,
                     zero_na = input$dp5,
                     log_tr = input$dp6,
                     base = input$dp7
                    )
    
    #Correlation plots  
    output$crp6 <- renderPlot({
      # Let's first check the correlation between tech.replicates 1 and 2
      corr_plot(raw, rep_1 = input$crp1, rep_2 = input$crp2, text_size = input$crp3, dpi = input$crp4, palette = input$crp5)})
    #Download Cor plot
    output$crp8 <- downloadHandler(
      filename = function() {
        paste("Correlation_plots",input$crp7, sep = ".")
      },
      content = function(file){
        if(input$crp7 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(corr_plot(raw, rep_1 = input$crp1, rep_2 = input$crp2, text_size = input$crp3, dpi = input$crp4, palette = input$crp5))
        dev.off()
      })
    
    #Data Quality Control
    if (input$dp4 == "TRUE" && input$dq3 == "Yes"){
      raw <- rem_sample(raw, input$dq4)
      raw_ave <- aver_techreps(raw)
      raw_filtered <- filterbygroup_na(raw_ave, set_na = input$dq1, filter_condition = input$dq2)
    }else if (input$dp4 == "TRUE"){
      raw_ave <- aver_techreps(raw)
      raw_filtered <- filterbygroup_na(raw_ave, set_na = input$dq1, filter_condition = input$dq2)
    }else{
      #Filter out proteins with high levels of missing data in each condition
      raw_filtered <- filterbygroup_na(raw, set_na = input$dq1, filter_condition = input$dq2)
    }
      
    #Visualize missing data in a subset of proteins.
    output$prange <- renderText({ nrow(raw_filtered) })
    output$srange <- renderText({ ncol(raw_filtered) })
    output$dqh13 <- renderPlot({
      heatmap_na(raw_filtered, protein_range = input$dqh1:input$dqh2, sample_range = input$dqh3:input$dqh4, reorder_x = input$dqh5, reorder_y = input$dqh6, x_fun = input$dqh7, y_fun = input$dqh8, label_proteins = input$dqh9, text_size = input$dqh10, dpi = input$dqh11, palette = input$dqh12)})
    #Download Heatmap
    output$dqh15 <- downloadHandler(
      filename = function() {
        paste("Missing_data_Heatmap",input$dqh14, sep = ".")
      },
      content = function(file){
        if(input$dqh14 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(heatmap_na(raw_filtered, protein_range = input$dqh1:input$dqh2, sample_range = input$dqh3:input$dqh4, reorder_x = input$dqh5, reorder_y = input$dqh6, x_fun = input$dqh7, y_fun = input$dqh8, label_proteins = input$dqh9, text_size = input$dqh10, dpi = input$dqh11, palette = input$dqh12))
        dev.off()
      })
    
    #Impute missing data
    if (input$dq5 == "minProb"){
      imp_df <- impute_na(raw_filtered, method = input$dq5, tune_sigma = input$dq6, q = input$dq7, seed = input$dq8)
      }else if (input$dq5 == "minDet"){
        imp_df <- impute_na(raw_filtered, method = input$dq5, q = input$dq9, seed = input$dq10)
      }else if (input$dq5 == "RF"){
        imp_df <- impute_na(raw_filtered, method = input$dq5, maxiter = input$dq11, ntree = input$dq12, seed = input$dq13)
      }else if (input$dq5 == "kNN"){
        imp_df <- impute_na(raw_filtered, method = input$dq5, seed = input$dq14)
      }else if (input$dq5 == "SVD"){
        imp_df <- impute_na(raw_filtered, method = input$dq5, n_pcs = input$dq15, seed = input$dq16)
      }
    
    output$dqim7 <- renderPlot({
      impute_plot(original = raw_filtered, imputed = imp_df, global = input$dqim1, n_row = input$dqim2, n_col = input$dqim3, text_size = input$dqim4, dpi = input$dqim5, palette = input$dqim6)})
    #Download Impute missing data
    output$dqim9 <- downloadHandler(
      filename = function() {
        paste("Impute_missing_data",input$dqim8, sep = ".")
      },
      content = function(file){
        if(input$dqim8 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(impute_plot(original = raw_filtered, imputed = imp_df, global = input$dqim1, n_row = input$dqim2, n_col = input$dqim3, text_size = input$dqim4, dpi = input$dqim5, palette = input$dqim6))
        dev.off()
      })
    
    #Normalize intensity data
    norm_df <- normalize_data(imp_df, method = input$dq17)
    output$dqn5 <- renderPlot({
      norm_plot(original = imp_df, normalized = norm_df, type = input$dqn1, text_size = input$dqn2, dpi = input$dqn3, palette = input$dqn4)})
    #Download intensity data
    output$dqn7 <- downloadHandler(
      filename = function() {
        paste("Normalization_data",input$dqn6, sep = ".")
      },
      content = function(file){
        if(input$dqn6 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(norm_plot(original = imp_df, normalized = norm_df, type = input$dqn1, text_size = input$dqn2, dpi = input$dqn3, palette = input$dqn4))
        dev.off()
      })
    
    find_depp <- function(norm_df,
      cutoff = input$dex1,
      lfc = input$dex2,
      adjustmethod = input$dex3){
      
      # Extract group information from colnames
      group <- factor(c(sapply(
        strsplit(colnames(norm_df), "_"),
        getElement, 1
      )))
      
      # create a design based on groups
      design <- model.matrix(~group)
      
      # Fit the model to the protein intensity data based on the experimental design
      fit <- limma::lmFit(norm_df, design)
      fit <- limma::eBayes(fit,
                           robust = T,
                           trend = T
      )
      
      # Make a a list of DE results based on provided criteria
      dec_test <- limma::decideTests(fit,
                                     lfc = lfc,
                                     adjust.method = adjustmethod
      )
    return(fit)
    }
    
    if (input$dp8 == "TRUE"){
    
    #Differentially expressed proteins with Gene name
    output$contents <- DT::renderDataTable({
      cutoff = input$dex1
      lfc = input$dex2
      adjustmethod = input$dex3
      adjustt = input$dex3
      limmaR = input$dex4
      tophits = input$dex5
      n_top = input$dex6
      
      # Extract group information from colnames
      group <- factor(c(sapply(
        strsplit(colnames(norm_df), "_"),
        getElement, 1
      )))
      
      # create a design based on groups
      design <- model.matrix(~group)
      
      # Fit the model to the protein intensity data based on the experimental design
      fit <- limma::lmFit(norm_df, design)
      fit <- limma::eBayes(fit,
                           robust = T,
                           trend = T
      )
      
      # Make a a list of DE results based on provided criteria
      dec_test <- limma::decideTests(fit,
                                     lfc = lfc,
                                     adjust.method = adjustmethod
      )
      
      
      # Write the results of the DE analysis to a text file (tab-separated)
      if (limmaR == TRUE) {
        fileName <- file.path(tempdir(), "limmaF.csv")
        
        # Write the results of the DE analysis to a text file (tab-separated)
        limma::write.fit(fit,
                         file = fileName,
                         adjust = adjustt,
                         results = dec_test,
                         sep="\t"
        )
        # Download DE analysis
        output$downtab2 <- downloadHandler(
          filename = function() {
            paste("Limma_Output",".csv", sep = "")
          },
          content = function(file){
            write.csv(read.csv(fileName, sep="\t"), file, row.names = FALSE)
          })
        
        files <- dir(path=tempdir(), pattern="limmaF*")
        unlink(x=files)
      }
      
      results_de <- limma::topTable(fit,
                                    coef = 2,
                                    adjust.method = adjustmethod,
                                    n = Inf
      )
      
      # add majority protein ids column
      results_de$gene_names <- rownames(results_de)
      rownames(results_de) <- NULL
      
      # rearrange order of columns
      results_de <- results_de[, c(
        "gene_names",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B"
      )]
      
      # extract proteins with absolute logfc > lfc
      results_de <- results_de[abs(results_de$logFC) > lfc, ]
      
      # extract sig. de. proteins and order from smallest to largest p. values
      results_de <- results_de[results_de$adj.P.Val < cutoff, ]
      results_de <- results_de[order(results_de$P.Value, results_de$adj.P.Val), ]
      
      if (nrow(results_de) == 0) {
        stop(
          message(
            paste0(
              "No differentially expressed proteins found at adj.P.value cutoff = ",
              cutoff
            )
          )
        )
      } else {
        message(paste0(
          nrow(results_de),
          " siginificantly differentially expressed proteins found."
        ))
      }
      
      output$boxres <- renderInfoBox({
        infoBox("", paste0(nrow(results_de), " ", "Siginificantly differentially expressed proteins found."), color = "purple", fill = TRUE)
      })
      
      if (tophits == TRUE) {
        if (nrow(results_de) < n_top) {
          fileName2 <- file.path(tempdir(), "TOPH1.csv")
          write.table(results_de[seq_len(nrow(results_de)), ],
                      file = fileName2,
                      quote = FALSE,
                      sep="\t"
          )
          # Download DE analysis
          output$downtab3 <- downloadHandler(
            filename = function() {
              paste("Top_Hits_Results_of_Differential_expression_analysis",".csv", sep = "")
            },
            content = function(file){
              write.csv(read.csv(fileName2, sep="\t"), file, row.names = FALSE)
            })
          
          files1 <- dir(path=tempdir(), pattern="TOPH1*")
          unlink(x=files1)
        } else {
          fileName3 <- file.path(tempdir(), "TOPH2.csv")
          write.table(results_de[1:n_top, ],
                      file = fileName3,
                      quote = FALSE,
                      sep = "\t"
          )
          # Download DE analysis
          output$downtab3 <- downloadHandler(
            filename = function() {
              paste("Top_Hits_Results_of_Differential_expression_analysis",".csv", sep = "")
            },
            content = function(file){
              write.csv(read.csv(fileName3, sep="\t"), file, row.names = FALSE)
            })
          
          files2 <- dir(path=tempdir(), pattern="TOPH2*")
          unlink(x=files2)
        }
      }
      
      #Download all data as .csv file
      output$downtab <- downloadHandler(
        filename = function() {
          paste("Differential_expression_analysis_all_data",".csv", sep = "")
        },
        content = function(file){
          write.csv(results_de, file, row.names = FALSE)
        })
      return(results_de)
    },
    rownames= FALSE,
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1)))))
    
  }
  else
    #Differentially expressed proteins with default Protein IDs results
    output$contents <- DT::renderDataTable({
      cutoff = input$dex1
      lfc = input$dex2
      adjustmethod = input$dex3
      adjustt = input$dex3
      limmaR = input$dex4
      tophits = input$dex5
      n_top = input$dex6
      
      # Extract group information from colnames
      group <- factor(c(sapply(
        strsplit(colnames(norm_df), "_"),
        getElement, 1
      )))
      
      # create a design based on groups
      design <- model.matrix(~group)
      
      # Fit the model to the protein intensity data based on the experimental design
      fit <- limma::lmFit(norm_df, design)
      fit <- limma::eBayes(fit,
                           robust = T,
                           trend = T
      )
      
      # Make a a list of DE results based on provided criteria
      dec_test <- limma::decideTests(fit,
                                     lfc = lfc,
                                     adjust.method = adjustmethod
      )
      
      
      # Write the results of the DE analysis to a text file (tab-separated)
      if (limmaR == TRUE) {
        fileName <- file.path(tempdir(), "limmaF.csv")
        
        # Write the results of the DE analysis to a text file (tab-separated)
        limma::write.fit(fit,
                         file = fileName,
                         adjust = adjustt,
                         results = dec_test,
                         sep="\t"
        )
        # Download DE analysis
        output$downtab2 <- downloadHandler(
          filename = function() {
            paste("Limma_Output",".csv", sep = "")
          },
          content = function(file){
            write.csv(read.csv(fileName, sep="\t"), file, row.names = FALSE)
          })
        
        files <- dir(path=tempdir(), pattern="limmaF*")
        unlink(x=files)
      }
      
      results_de <- limma::topTable(fit,
                                    coef = 2,
                                    adjust.method = adjustmethod,
                                    n = Inf
      )
     
      # add majority protein ids column
      results_de$majority_protein_id <- rownames(results_de)
      rownames(results_de) <- NULL
      
      # rearrange order of columns
      results_de <- results_de[, c(
        "majority_protein_id",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B"
      )]
      
      # extract proteins with absolute logfc > lfc
      results_de <- results_de[abs(results_de$logFC) > lfc, ]
      
      # extract sig. de. proteins and order from smallest to largest p. values
      results_de <- results_de[results_de$adj.P.Val < cutoff, ]
      results_de <- results_de[order(results_de$P.Value, results_de$adj.P.Val), ]
      
      if (nrow(results_de) == 0) {
        stop(
          message(
            paste0(
              "No differentially expressed proteins found at adj.P.value cutoff = ",
              cutoff
            )
          )
        )
      } else {
        message(paste0(
          nrow(results_de),
          " siginificantly differentially expressed proteins found."
        ))
      }
      
      output$boxres <- renderInfoBox({
        infoBox("", paste0(nrow(results_de), " ", "Siginificantly differentially expressed proteins found."), color = "purple", fill = TRUE)
        })
      
      if (tophits == TRUE) {
        if (nrow(results_de) < n_top) {
          fileName2 <- file.path(tempdir(), "TOPH1.csv")
          write.table(results_de[seq_len(nrow(results_de)), ],
                      file = fileName2,
                      quote = FALSE,
                      sep="\t"
          )
          # Download DE analysis
          output$downtab3 <- downloadHandler(
            filename = function() {
              paste("Top_Hits_Results_of_Differential_expression_analysis",".csv", sep = "")
            },
            content = function(file){
              write.csv(read.csv(fileName2, sep="\t"), file, row.names = FALSE)
            })
          
          files1 <- dir(path=tempdir(), pattern="TOPH1*")
          unlink(x=files1)
        } else {
          fileName3 <- file.path(tempdir(), "TOPH2.csv")
          write.table(results_de[1:n_top, ],
                      file = fileName3,
                      quote = FALSE,
                      sep = "\t"
          )
          # Download DE analysis
          output$downtab3 <- downloadHandler(
            filename = function() {
              paste("Top_Hits_Results_of_Differential_expression_analysis",".csv", sep = "")
            },
            content = function(file){
              write.csv(read.csv(fileName3, sep="\t"), file, row.names = FALSE)
            })
          
          files2 <- dir(path=tempdir(), pattern="TOPH2*")
          unlink(x=files2)
        }
      }
      
      #Download all data as .csv file
      output$downtab <- downloadHandler(
        filename = function() {
          paste("Differential_expression_analysis_all_data",".csv", sep = "")
        },
        content = function(file){
          write.csv(results_de, file, row.names = FALSE)
        })
      return(results_de)
    },
    rownames= FALSE,
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1)))))
  
    #DE
    fit_df <- find_depp (norm_df, cutoff = input$dex1, lfc = input$dex2, adjustmethod = input$dex3)
    
    #Volcano Plot
    output$rvp8 <- renderPlot({
      volcano_plot(fit_df, adj_method = input$dex3, cutoff = input$dex1, lfc = input$dex2, n_top = input$dex6, sig = input$rvp1, line_fc = input$rvp2, line_p = input$rvp3, label_top = input$rvp4, text_size = input$rvp5, dpi = input$rvp6, palette = input$rvp7)})
      #Download Volcano Plot
    output$rvp10 <- downloadHandler(
      filename = function() {
        paste("Volcano_plot",input$rvp9, sep = ".")
      },
      content = function(file){
        if(input$rvp9 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(volcano_plot(fit_df, adj_method = input$dex3, cutoff = input$dex1, lfc = input$dex2, n_top = input$dex6, sig = input$rvp1, line_fc = input$rvp2, line_p = input$rvp3, label_top = input$rvp4, text_size = input$rvp5, dpi = input$rvp6, palette = input$rvp7))
        dev.off()
      })
    
    #Heat Map DE
    output$rhp5 <- renderPlot({
      heatmap_de(fit_df, norm_df, cutoff = input$dex1, lfc = input$dex2, n_top = input$dex6, sig = input$rhp1, text_size = input$rhp2, dpi = input$rhp3, palette = input$rhp4)})
    #Download Heat Map DE
    output$rhp7 <- downloadHandler(
      filename = function() {
        paste("Heat_map_DE_plot",input$rhp6, sep = ".")
      },
      content = function(file){
        if(input$rhp6 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(heatmap_de(fit_df, norm_df, cutoff = input$dex1, lfc = input$dex2, n_top = input$dex6, sig = input$rhp1, text_size = input$rhp2, dpi = input$rhp3, palette = input$rhp4))
        dev.off()
      })
    
    #Protein Data Modeling
    pre_process1 <- function(fit_df, norm_df,
    sig = input$pmp1,
    sig_cutoff = input$pmp2,
    fc = input$pmp3,
    n_top = input$pmp4,
    find_highcorr = input$pmp5,
    corr_cutoff = input$pmp6,
    rem_highcorr = input$pmp7,
    save_corrmatrix = input$pmp8){
    
    adjustmd = input$pmp9
    #binding for global variable
    logFC <- P.Value <- adj.P.Val <-  NULL
    
    # Extract the results from the differential expression analysis.
    exp_de <- limma::topTable(fit_df,
                              coef = colnames(fit_df)[2],
                              n = length(fit_df$df.total),
                              adjust.method = adjustmd
    )
    
    # Subset results by logFC and p-value cutoff
    if (sig == "P") {
      top_proteins <- rownames(subset(exp_de,
                                      abs(logFC) > fc & P.Value < sig_cutoff,
                                      drop = FALSE
      ))
      
      # Or default: based on adj.P value
    } else {
      top_proteins <- rownames(subset(exp_de,
                                      abs(logFC) > fc & adj.P.Val < sig_cutoff,
                                      drop = FALSE
      ))
    }
    
    # If the total number of DE proteins < n_top, replace n_top with that number.
    if (length(top_proteins) < n_top) {
      n_top <- length(top_proteins)
      message(paste0(
        "Total number of differentially expressed proteins (",
        n_top, ") ", "is less than n_top."
      ))
    }
    
    # Extract the top n_top hits from the top hit list
    top_proteins <- top_proteins[1:n_top]
    
    # Check if there are sig. proteins before moving on to pre-processing
    if (identical(top_proteins, character(0))) {
      stop(message
           (paste0(
             "No significant proteins found at ",
             sig,
             " < ",
             sig_cutoff,
             "."
           )))
    } else {
      # Extract intensity values for top proteins based on fc and sig_cutoff
      top_intensity <- subset(norm_df,
                              rownames(norm_df) %in% top_proteins,
                              drop = FALSE
      )
    }
    # Extract group or condition information from sample names in the data frame
    group <- factor(c(sapply(
      strsplit(colnames(top_intensity), "_"),
      getElement, 1
    )))
    
    # Transpose the data frame. Columns are now proteins and rows are samples.
    topint_trans <- as.data.frame(t(top_intensity))
    
    # Remove sample names.
    rownames(topint_trans) <- NULL
    
    # Add a new column with the group or condition information.
    # condition column is now the rightmost column in the data frame.
    topint_trans$condition <- group
    
    # For correlation calculations, make a matrix without the condition column
    topint_cor <- topint_trans[, seq_len(ncol(topint_trans)) - 1]
    
    # Create a correlation matrix
    cor_matrix <- cor(topint_cor)
    
    if (save_corrmatrix == TRUE) {
      fileName4 <- file.path(tempdir(), "PMD1.csv")
      write.table(cor_matrix,
                  file = fileName4,
                  row.names = FALSE,
                  col.names = FALSE,
                  sep = "\t",
                  quote = FALSE
      )
      # Download DE analysis
      output$downtab4 <- downloadHandler(
        filename = function() {
          paste("Protein_correlation_matrix",".csv", sep = "")
        },
        content = function(file){
          write.csv(read.csv(fileName4, sep="\t"), file, row.names = FALSE)
        })
      files4 <- dir(path=tempdir(), pattern="PMD1*")
      unlink(x=files4)
    }
    
    if (find_highcorr == TRUE) {
      # Identify protein columns with high pairwise-correlation to remove
      highcor <- findCorrelation(cor_matrix, cutoff = corr_cutoff, names = TRUE)
      if (length(highcor != 0)) {
        message(
          "Following protein(s) show high pairwise-correlation"
        )
      } else {
        message(
          "None of the proteins show high pair-wise correlation."
        )
      }
      message(paste0(highcor, collapse = "\n"))
      
      if (rem_highcorr == TRUE) {
        topint_trans_1 <- topint_trans[, !(colnames(topint_trans) %in% highcor)]
        if (ncol(topint_trans_1) == ncol(topint_trans)) {
          message("No highly correlated proteins to be removed.")
        } else {
          message("Proteins with high pairwise-correlation have been removed.")
        }
      } else {
        warning("Proteins with high pairwise-correlation have NOT been removed.",
                call. = FALSE
        )
        topint_trans_1 <- topint_trans
      }
    } else {
      warning("Your data could have proteins with high pairwise-correlation.",
              call. = FALSE
      )
      topint_trans_1 <- topint_trans
    }
    
    # Convert condition names to R compatible names
    topint_trans_1$condition <- make.names(topint_trans_1$condition)
    
    # Convert condition to a factor (important for varimp calculations)
    topint_trans_1$condition <- factor(topint_trans_1$condition)
    return(topint_trans_1)
    }

    model_df <- pre_process1(fit_df, norm_df, sig = input$pmp1, sig_cutoff = input$pmp2, fc = input$pmp3, n_top = input$pmp4, find_highcorr = input$pmp5, corr_cutoff = input$pmp6, rem_highcorr = input$pmp7, save_corrmatrix = input$pmp8)
    
    #Feature Plot
    output$pmf7 <- renderPlot({
      feature_plot(model_df, type = input$pmf1, n_row = input$pmf2, n_col = input$pmf3, text_size = input$pmf4, dpi = input$pmf5, palette = input$pmf6)})
    #Download Feature Plot
    output$pmf9 <- downloadHandler(
      filename = function() {
        paste("Feature_plot",input$pmf8, sep = ".")
      },
      content = function(file){
        if(input$pmf8 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(feature_plot(model_df, type = input$pmf1, n_row = input$pmf2, n_col = input$pmf3, text_size = input$pmf4, dpi = input$pmf5, palette = input$pmf6))
        dev.off()
      })
    
    #Split Data
    split_df <- split_data(model_df, train_size = input$pms1, seed = input$pms2)
    
    #Model list
    model_list <- train_models(split_df, resample_method = input$pmpmt1, num_repeats = input$pmpmt2, resample_iterations = input$pmpmt3, seed = input$pmpmt4, algorithm_list = input$pmpmt5)
    
    #Performance plot
    output$pmpr5 <- renderPlot({
      performance_plot(model_list, type = input$pmpr1, text_size = input$pmpr2, dpi = input$pmpr3, palette = input$pmpr4)})
    #Download performance Plot
    output$pmpr7 <- downloadHandler(
      filename = function() {
        paste("Performance_plot",input$pmpr6, sep = ".")
      },
      content = function(file){
        if(input$pmpr6 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(performance_plot(model_list, type = input$pmpr1, text_size = input$pmpr2, dpi = input$pmpr3, palette = input$pmpr4))
        dev.off()
      })
    
    #Variable importance plot
    output$pmvp7 <- renderPlot({
      varimp_plot(model_list, type = input$pmvp1, n_row = input$pmvp2, n_col = input$pmvp3, text_size = input$pmvp4, dpi = input$pmvp5, palette = input$pmvp6)})
    #Download Variable importance plot
    output$pmvp9 <- downloadHandler(
      filename = function() {
        paste("Variable_importance_plot",input$pmvp8, sep = ".")
      },
      content = function(file){
        if(input$pmvp8 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(varimp_plot(model_list, type = input$pmvp1, n_row = input$pmvp2, n_col = input$pmvp3, text_size = input$pmvp4, dpi = input$pmvp5, palette = input$pmvp6))
        dev.off()
      })
    
    # test models
    test_modelss <- function(model_list,
                             split_df,
                             type = input$pmtm1,
                             save_confusionmat = input$pmtm2
    ) {
      
      # Extract test data from the split_df object
      test_data <- split_df$test
      
      # Predict test data
      pred_list <- lapply(
        model_list,
        function(x) {
          message(paste0(
            "\n",
            "Testing ",
            x$method,
            "...",
            "\n"
          ))
          predict(x,
                  test_data,
                  type = type
          )
        }
      )
      message(paste0("\n", "Done!"))
      
      # Get confusion matrices and associated statistics
      if (type == "raw") {
        cm_list <- lapply(
          pred_list,
          function(x) {
            confusionMatrix(
              x,
              test_data$condition
            )
          }
        )
        
        if (save_confusionmat == TRUE) {
          
          # Convert c.matrices to long-form data frames
          cm_df <- lapply(
            cm_list,
            function(x) reshape2::melt(as.table(x))
          )
          # Get the list of methods
          method_list <- names(cm_df)
          
          # Add method names to the data frames
          cm_dfm <- lapply(
            seq_along(method_list),
            function(x) {
              cm_df[[x]]["method"] <- method_list[x]
              cm_df[[x]]
            }
          )
          
          # Combine all data frames into one
          cm_dfm_long <- do.call(
            "rbind",
            cm_dfm
          )
          
          # Add column names before saving
          colnames(cm_dfm_long) <- c("Prediction", "Reference", "Value", "Method")
          
          fileName5 <- file.path(tempdir(), "CONFMM.csv")
          
          # Save data in a text file
          write.table(cm_dfm_long,
                      file = fileName5,
                      quote = FALSE,
                      row.names = FALSE,
                      sep = "\t"
          )
          
          # Download confusion matrix
          output$downtab5 <- downloadHandler(
            filename = function() {
              paste("Confusion_matrices",".csv", sep = "")
            },
            content = function(file){
              write.csv(read.csv(fileName5, sep="\t"), file, row.names = FALSE)
            })
          
          files5 <- dir(path=tempdir(), pattern="CONFMM*")
          unlink(x=files5)
        }
      }
      
      return(pred_list)
    }
    
    prob_list <- test_modelss(model_list, split_df, type = input$pmtm1, save_confusionmat = input$pmtm2)
    
    #ROC plot
    output$pmroc5 <- renderPlot({
      roc_plot(prob_list, split_df, multiple_plots = input$pmroc1, text_size = input$pmroc2, dpi = input$pmroc3, palette = input$pmroc4)})
    #Download ROC plot
    output$pmroc7 <- downloadHandler(
      filename = function() {
        paste("ROC_plot",input$pmroc6, sep = ".")
      },
      content = function(file){
        if(input$pmroc6 == "png")
          png(file) # open the png device
        else
          pdf(file) # open the pdf device
        print(roc_plot(prob_list, split_df, multiple_plots = input$pmroc1, text_size = input$pmroc2, dpi = input$pmroc3, palette = input$pmroc4))
        dev.off()
      })

    #final
    })
}

shinyApp(ui, server)