# NMR QC Tool
#disclaimer and link to paper for conditions




library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(dplyr)
library(DT)
library(plotly)
library(ggplot2)
library(colourpicker)
library(rmarkdown)
library(forcats)
library(tibble)
library(lubridate)
library(ggrepel)
library(treemapify)
library(viridis)

#data input
LM_pre<-read.csv("Precent_LM_all_final.csv",stringsAsFactors = FALSE)
LM_post<-read.csv("Postcent_LM_all_final.csv",stringsAsFactors = FALSE)
LM_pre_error<-read.csv("Precent_error_final.csv",stringsAsFactors = FALSE)
LM_post_error<-read.csv("Postcent_error_final.csv",stringsAsFactors = FALSE)

# Scripts -----------------------------------------------------------------
#Calculator, calculates the % based on the LME models and using the inputs of pre/post centrifugation times
calculator<-function(data_pre,data_post,t1,t2){
  #Model pre
  LM_pre<-data_pre
  LM_post<-data_post

  #rearrange data frames and rename
  LM_pre_reduced<-data.frame(LM_pre$name,LM_pre$intercept,LM_pre$slope)
  LM_pre_reduced<-setNames(LM_pre_reduced,c("Metabolite","PreInt","PreSlope"))
  LM_post_reduced<-data.frame(LM_post$name,LM_post$intercept,LM_post$slope)
  LM_post_reduced<-setNames(LM_post_reduced,c("Metabolite","PostInt","PostSlope"))
  LM_combo<-inner_join(LM_pre_reduced,LM_post_reduced,by="Metabolite")

  #calulate deltas
  LM_combo$x2<-t1*(LM_combo$PreSlope)+LM_combo$PreInt
  LM_combo$x3<-t2*(LM_combo$PostSlope)+LM_combo$x2
  LM_combo$delta1<-(LM_combo$x2)-(LM_combo$PreInt)
  LM_combo$delta1_percent<-((LM_combo$delta1))*100/(LM_combo$PreInt)
  LM_combo$delta2<-((LM_combo$x3)-(LM_combo$x2))
  LM_combo$delta2_percent<-((LM_combo$delta2))*100/(LM_combo$x2)
  LM_combo$combined_percent<-LM_combo$delta1_percent+LM_combo$delta2_percent

  #combine final data frame
  final<-data.frame(LM_combo$Metabolite,LM_combo$delta1_percent,LM_combo$delta2_percent,LM_combo$combined_percent)
  final<-setNames(final,c("Metabolite","PreCent","PostCent","Combined"))
}

calculator_multi<-function(df,data_pre,data_post){

  LM_pre<-data_pre
  LM_post<-data_post
  LM_pre_reduced<-data.frame(LM_pre$name,LM_pre$intercept,LM_pre$slope)
  LM_pre_reduced<-setNames(LM_pre_reduced,c("Metabolite","PreInt","PreSlope"))
  LM_post_reduced<-data.frame(LM_post$name,LM_post$intercept,LM_post$slope)
  LM_post_reduced<-setNames(LM_post_reduced,c("Metabolite","PostInt","PostSlope"))
  LM_combo<-inner_join(LM_pre_reduced,LM_post_reduced,by="Metabolite")

  LM_final<-as.data.frame(cbind(df$ID,df$PreCent,df$PostCent,df$TTF))
  LM_final<-setNames(LM_final,c("ID","PreCent","PostCent","TTF"))

  LM_combo<-LM_combo %>% column_to_rownames(var="Metabolite")

  temp_data<-data.frame()
  lm_testy<-data.frame()

  tempy<-list()
  x2<-list()
  x3<-list()
  delta1<-list()
  delta1_percent<-list()
  delta2<-list()
  delta2_percent<-list()
  combined_percent<-list()
  all_change<-list()
  pre_change<-list()
  post_change<-list()

  #double for loop
  LM_final<-LM_final %>% column_to_rownames(var="ID")

  for(n in rownames(LM_final)){
    for(i in rownames(LM_combo)){
      x2[i]<-as.numeric(LM_final[n,]$PreCent)*LM_combo[i,]$PreSlope+LM_combo[i,]$PreInt
      x3[i]<-as.numeric(LM_final[n,]$PostCent)*(LM_combo[i,]$PostSlope)+x2[[i]]
      delta1[i]<-(x2[[i]])-(LM_combo[i,]$PreInt)
      delta1_percent[i]<-((delta1[[i]]))*100/(LM_combo[i,]$PreInt)
      delta2[i]<-((x3[[i]])-(x2[[i]]))
      delta2_percent[i]<-((delta2[[i]]))*100/(x2[[i]])
      combined_percent[i]<-delta1_percent[[i]]+delta2_percent[[i]]
    }
    pre_change[n]<-list(change=delta1_percent)
    post_change[n]<-list(change=delta2_percent)
    all_change[n]<-list(change = combined_percent)
    all_change_dataframe<-as.data.frame(t(do.call(cbind,all_change)))
    all_change_dataframe<-tibble::rownames_to_column(all_change_dataframe,"ID")
  }
  LM_final<-tibble::rownames_to_column(LM_final,"ID")
  final_df<-inner_join(LM_final,all_change_dataframe,by="ID")
}


date_helper<-function(df,sprec,centtime){
  #create time tables from csv using lubridate

  df$Centrifugation<-parse_date_time(df$Centrifugation,order="dmy HM")
  df$Freeze<-parse_date_time(df$Freeze,order="dmy HM")
  df$Draw<-parse_date_time(df$Draw,order="dmy HM")
  df$PreCent<-as.numeric(difftime(df$Centrifugation,df$Draw),units="hours")
  df$PostCent<-as.numeric(difftime(df$Freeze,df$Centrifugation),units="hours")
  df$TTF<-as.numeric(difftime(df$Freeze,df$Draw),units="hours")
  df$PostCent<-df$PostCent-centtime/60
  if(sprec=="Yes"){
    df$SPREC_PRE<-ifelse(df$PreCent < 0.5,"A1",
                       ifelse(df$PreCent < 2, "A",
                              ifelse(df$PreCent < 4, "C",
                                     ifelse(df$PreCent < 8, "E",
                                            ifelse(df$PreCent < 12, "G",
                                                   ifelse(df$PreCent < 24, "I",
                                                          ifelse(df$PreCent < 48,"K","M")
                                                   )
                                            )
                                     )
                              )
                       )
   )
   df$SPREC_POST<-ifelse(df$PostCent < 1,"B",
                        ifelse(df$PostCent <2, "D",
                               ifelse(df$PostCent <8,"F",
                                      ifelse(df$PostCent < 24, "H",
                                             ifelse(df$PostCent < 48, "J","M")
                                             )
                                      )
                              )
                        )
 }
  df<-df
}



# UI ----------------------------------------------------------------------

# Menu --------------------------------------------------------------------


ui <- dashboardPage(
    dashboardHeader(title="PRIMA-Panel"),

    dashboardSidebar(
      sidebarMenu(
                  collapsed=FALSE,
                  menuItem("Home", tabName = "home", icon = icon("home"),selected=T),
                  menuItem("Data",tabName="blood_data_plots",icon=icon(("newspaper")),
                    menuSubItem("Pre-Centrifugation",tabName="blood_data_pre",icon=icon(("right-to-bracket"))), #alternative arrow-right-to-bracket, and arrow-left-to-bracket, right-from-bracket.left-from-bracket
                    menuSubItem("Post-Centrifugation",tabName="blood_data_post",icon=icon(("right-from-bracket")))
                  ),
                  menuItem("QC Panels",tabName="blood_tables",icon=icon(("table")),
                    menuSubItem("Single Sample",tabName="blood_tables",icon=icon("folder")),
                    menuSubItem("Batch",tabName="blood_tables_batch",icon=icon("folder-tree"))
                  ),
                  menuItem("Tools",icon=icon("toolbox"),
                           menuSubItem("Delay Calculator",tabName="blood_tables_helper",icon=icon("clock"))
                           )
      )
    ),

# Home --------------------------------------------------------------------


    dashboardBody(

      tags$head(tags$style(HTML('
            /* logo */
        .skin-blue .main-header .logo {
                              background-color: #ffffff;
                              color:#000000;
                              border-right: 3px solid #254264;

        }
                 /* logo when hovered */
        .skin-blue .main-header .logo:hover{
                              background-color: #ffffff;
                              color:#000000;

                              }
                   /* toggle button when hovered */
        .skin-blue .main-header .navbar .sidebar-toggle {
                              background-color: #ffffff;
                              color:#000000;
                              }
         /* toggle button when hovered */
        .skin-blue .main-header .navbar .sidebar-toggle:hover {
                              background-color: #ffffff;
                              color:#000000;
                              }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #ffffff;
                              color:#000000;

        }

        /* main sidebar */
        .skin-blue .main-sidebar {
                              background-color: #ffffff;
                              color:#000000;
                              border-right: 3px solid #254264;
        }
        /* active selected tab in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #ffffff;
                              color:#000000;
                              font-weight:bold;
                              }
        /* other links in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #ffffff;
                              color: #000000;
        }
        /* other links in the sidebarmenu when hovered */
         .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #ffffff;
                              color:#000000;
                              font-style:italic;
                              font-weight:bold;
         }

        /* box background
         .skin-blue .box.box-solid.box-primary>.box-header {
                              background-color: #ffffff;
                              color:#000000;
         }
        /* body */
        .content-wrapper, .right-side {
                              background-color: #ffffff;
        }
        .box{
        -webkit-box-shadow: none;
        -moz-box-shadow: none;
        box-shadow: none;}

      .box-info{border-top-style: none;}

        '))),





      tabItems(
          tabItem(tabName="home",
                  fluidRow(
                    box(width=12, status = "info",

                        br(),h1(HTML("<b>PRIMA</b>-Panel - <b>Pr</b>e-Analytical <b>I</b>nvestigator for NMR-based <b>M</b>et<b>a</b>bolomics")),br(),br(),
                        fluidRow(
                          box(width=2,status="info",
                              img(src='logo.png',width="63%",style="display: block; margin-left: auto; margin-right: auto;"),
                          ),
                          box(title=HTML("<b>Author</b>"),status="info",width=4,
                            p("Developed by",a("Alexander Funk",href="https://www.uniklinikum-dresden.de/de/das-klinikum/kliniken-polikliniken-institute/klinische-chemie-und-laboratoriumsmedizin/forschung/copy_of_EMS")),
                            p("Institute for Clinical Chemistry and Laboratory Medicine"),
                            p("University Hospital Dresden, Fetscherstr. 74, 01307 Dresden")

                          ),
                          box(title=HTML("<b>Links</b>"),status="info",width=2,
                            p(a("GitHub",href="https://github.com/funkam/QC-Tool")),
                            p(a("shinyapps.io",href="https://funkam.shinyapps.io/QC-Tool/")),
                            p(a("Publication",href=""))
                          ),

                          box(width=4,status="info",
                              #img(src='logo.png',width="20%",style="float:left"),
                              img(src='qmp.png',width="40%",style="float:left"),
                              img(src='uni_logo.jfif',width="50%",style="float:left"),

                          )
                        ),
                  fluidRow(
                    box(title=HTML("<b>The PRIMA Panel</b>"),status="warning",width=12,
                        p("The PRIMA-Panel is tool to investigate the effect processing delays on metabolic parameters in samples of peripheral blood (plasma / serum)."),
                        p("Linear mixed models were used to estimate the change for each metabolic parameter. The data is split into pre- and postcentrinfugation delays"),
                        p("The data is presented in so called stability timepoints. Such a timepoint is defined as the time it takes for a parameter to change by a specific percentage. For example: A value of 0.2 for lactic acid in serum would mean it would take 0.2 hours when the % threshold is set to 20 % change."),
                        p("Further, the effect of these processing delays can be explored entering real pre-analytical data and observing the direct effect of the dealys on the parameterts. Here, interactive HTML reports can be created.")
                    )
                  ),
                  fluidRow(
                    box(title=HTML("<b>Modules</b>"),status="success",width=12,
                    box(title=HTML("<b>Data</b>"),status="success",width=4,
                        p("The Data tab shows different ways of highlighting the different stability time-points in minutes. The time-points are sorted according to their SPREC classification. In addition, there is an alternative way of presenting the date in form of lollipop plots."),
                        p("The data is split according the two different delays (pre- and post-centrifugation")
                        ),
                    box(title=HTML("<b>QC-Panel</b>"),status="success",width=4,
                        p("The 'Single' tab allows the user to set a pre-centrifugation time and a post-centrifugation using a Slider. A table is then generated that higlights a minor and a major change for each metabolite in that given timeframe. The colors, as well as the % threshholds, can be adjusted. The table can be downloaded as a .HTML report or the table directly as .csv (or other formats)."),
                        p("The batch tab allows the same procedure but with an uploaded table with pre- and post-centrifugation times to allow batch processing")
                        ),
                    box(title=HTML("<b>Tools</b>"),status="success",width=4,
                        p("An additional tab for Tools is available. Currently it consists of a tool for calculating the pre- and post-centrifugation times from the differences of date+time stamps. However a specific format is needed. See examples.")
                        )
                  )
                  )
                  )
                    )
          ),

# Pre ---------------------------------------------------------------------


          tabItem(tabName="blood_data_pre",
                  fluidRow(
                    column(width=3,
                      box(width=10, title=HTML("<b>Pre-Centrifugation Input</b>"), status="primary",
                                  numericInput("pb_plots_pre_percent","Enter threshhold (%)",value=20),
                                  numericInput("pb_plots_pre_cutoff","Enter plot cutoff (h)",value=10)
                        )
                    ),
                    column(width=3,
                        box(width=12,title=HTML("<b>Info</b>"),status="primary",
                            p("Enter the threshhold for the stability (i.e. 30% means a 30% change from it original value at 0 hours."),
                            p("The plot cutoff can be used to declutter the lollipop plots"),
                            p("The table displays the different parameters according to the SPREC caterogies.")
                        )
                    )
                    ),
                  fluidRow(
                    tabBox(width=12,
                      tabPanel(title="Table",
                               fluidRow(
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: A1 (<30min) ",
                                     DT::dataTableOutput('pb_table_pre_a1')
                                 ),
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: A (<2h)",
                                     DT::dataTableOutput('pb_table_pre_a')
                                 ),
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: C (2-4h)",
                                     DT::dataTableOutput('pb_table_pre_c')
                                 ),
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: E (4-8h)",
                                     DT::dataTableOutput('pb_table_pre_e')
                                 )
                               ),
                               fluidRow(
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: G (8-12h)",
                                     DT::dataTableOutput('pb_table_pre_g')
                                 ),
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: I (12-24h)",
                                     DT::dataTableOutput('pb_table_pre_i')
                                 ),
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: K (24-48h)",
                                     DT::dataTableOutput('pb_table_pre_k')
                                 ),
                                 box(width=3,solidHeader=TRUE,status="primary",title="SPREC: M (>48h)",
                                     DT::dataTableOutput('pb_table_pre_m')
                                 )
                               )
                      ),
                      tabPanel(title="Plots",
                               div(style="height:1000px",plotlyOutput("pb_lollis_pre"))

                      )
                    )
                  )

                  ),


# Post --------------------------------------------------------------------


          tabItem(tabName="blood_data_post",
                  fluidRow(
                    column(width=3,
                           box(width=10, title=HTML("<b>Post-Centrifugation Input</b>"), status="primary",
                               numericInput("pb_plots_post_percent","Enter threshhold (%)",value=20),
                               numericInput("pb_plots_post_cutoff","Enter plot cutoff (h)",value=10)
                           )
                    ),
                    column(width=3,
                           box(width=12,title=HTML("<b>Info</b>"),status="primary",
                               p("Enter the threshhold for the stability (i.e. 30% means a 30% change from it original value at 0 hours."),
                               p("The plot cutoff can be used to declutter the lollipop plots"),
                               p("The table displays the different parameters according to the SPREC caterogies.")
                           )
                    )
                  ),
                  fluidRow(
                    tabBox(width=12,
                           tabPanel(title="Table",
                                    fluidRow(
                                      box(width=3,solidHeader=TRUE,status="primary",title="SPREC: B (<1h) ",
                                          DT::dataTableOutput('pb_table_post_b')
                                      ),
                                      box(width=3,solidHeader=TRUE,status="primary",title="SPREC: D (1-2h)",
                                          DT::dataTableOutput('pb_table_post_d')
                                      ),
                                      box(width=3,solidHeader=TRUE,status="primary",title="SPREC: F (2-8h)",
                                          DT::dataTableOutput('pb_table_post_f')
                                      ),
                                      box(width=3,solidHeader=TRUE,status="primary",title="SPREC: H (8-24h)",
                                          DT::dataTableOutput('pb_table_post_h')
                                      )
                                    ),
                                    fluidRow(
                                      box(width=3,solidHeader=TRUE,status="primary",title="SPREC: J (24-48h)",
                                          DT::dataTableOutput('pb_table_post_j')
                                      ),
                                      box(width=3,solidHeader=TRUE,status="primary",title="SPREC: M (>48h)",
                                          DT::dataTableOutput('pb_table_post_m')
                                      ),
                                    )
                           ),
                           tabPanel(title="Plots",
                                    div(style="height:1000px",plotlyOutput("pb_lollis_post"))

                           )
                    )
                  )

          ),

# Single Table ------------------------------------------------------------



          tabItem(tabName="blood_tables",
                  fluidRow(
                          box(
                              width=2,title=HTML("<b>Sample Type</b>"),status="primary",
                                column(
                                    width=12,align="center",
                                    radioGroupButtons("sample_type","Select sample type:",choices=c("EDTA"="edta","Lithium-Heparin"="lihep","Serum"="serum"),selected="edta",direction="vertical")
                                )
                          ),
                          box(
                              width=6,title=HTML("<b>Style</b>"),status="success",
                              fluidRow(
                                column(width=3,
                                       box(width=12,
                                           numericInput("minor","Minor color %:",value=10),
                                           colourInput("minor_color","Minor Color:",value="orange"))
                                ),
                                column(width=4,
                                       box(width=12,
                                           numericInput("major","Major color %:",value=20),
                                           colourInput("major_color","Major Color:",value="red")
                                       )
                                ),
                                column(width=3,
                                       box(width=12,
                                           colourInput("neutral_color","Neutral Color:",value="")
                                       )
                                )
                              )
                          ),
                          box(width=2,title=HTML("<b>Tolerable Error</b>"),status="success",
                              p("Enter tolerable error in percent"),
                              numericInput("error","
                                           Error in %", value=15)),
                          box(
                            width=2,title=HTML("<b>Download Output</b>"),status="success",
                            p("Download a report as interactive HTML file.",align="justify"),
                            downloadButton("report", "Generate HTML report")
                          )

                  ),
                  fluidRow(
                          box(width=8,title=HTML("<b>Centrifugation Times</b>"),status="warning",
                              sliderInput("TTZ",
                                          "Pre-Centrifugation Time:",
                                          min = 0,
                                          max = 10,
                                          value = 0,
                                          step=0.1
                              ),
                              sliderInput("TTF",
                                          "Post-Centifugation Time:",
                                          min = 0,
                                          max = 10,
                                          value = 0,
                                          step=0.1)
                          )

                  ),
                  # fluidRow(
                  #   box(width=12,
                  #       DT::dataTableOutput('datatable3')
                  #   )
                  #       ),
                  # fluidRow(
                  #   box(width=12,
                  #       DT::dataTableOutput('datatable4')
                  #   )
                  # ),
                  # fluidRow(
                  #   box(width=12,
                  #       DT::dataTableOutput('datatable2')
                  #   )
                  # ),
                  fluidRow(
                          box(width=12,status="primary",
                            DT::dataTableOutput('datatable')
                          )
                  )
          ),
          # tabItem(tabName="blood_plots",
          #         sliderInput("input_time_ttz", "Adjust highest pre-centrifugation time shown:",
          #                     min=0,
          #                     max=24,
          #                     value=8,
          #                     step=0.1),
          #         sliderInput("input_time_ttf", "Adjust highest pre-centrifugation time shown:",
          #                     min=0,
          #                     max=24,
          #                     value=8,
          #                     step=0.1),
          #         plotlyOutput("pre_lolli")
          #         ),




# Batch Tables ------------------------------------------------------------
          tabItem(tabName="blood_tables_batch",
                  fluidRow(
                    column(width=3,
                           box(width=12,title=HTML("<b>Sample Type</b>"),status="primary",
                               fluidRow(

                                 column(
                                   width=6,align="center",
                                   radioGroupButtons("batch_sample_type","Select sample type:",choices=c("EDTA"="edta","Lithium-Heparin"="lihep","Serum"="serum"),selected="edta",direction="vertical")
                                 ),
                                 column(
                                   width=6,
                                   fileInput("batch_file", "File input", multiple=FALSE),
                                 )
                               ),
                               p("The input file should be a .CSV with the following column headers: ID, PreCent, PostCent. Ideally it was created with the Batch Table Helper"),
                           ),
                           box(
                             width=12,title=HTML("<b>Style</b>"),status="primary",
                             fluidRow(
                               column(width=6,
                                      box(width=12,
                                          numericInput("b_minor","Minor color %:",value=10),
                                          colourInput("b_minor_color","Minor Color:",value="orange")),
                                      box(width=12,
                                          colourInput("b_neutral_color","Neutral Color:",value="")
                                      )
                               ),
                               column(width=6,
                                      box(width=12,
                                          numericInput("b_major","Major color %:",value=20),
                                          colourInput("b_major_color","Major Color:",value="red")
                                      )
                               ),
                             )
                           ),
                           box(width=6,title=HTML("<b>Tolerable Error</b>"),status="primary",
                               p("Enter tolerable error in percent"),
                               numericInput("b_error","
                                           Error in %", value=15)),
                           box(
                             width=6,title=HTML("<b>Download Output</b>"),status="primary",
                             p("Download a report as interactive HTML file.",align="justify"),
                             downloadButton("batch_report", "Generate Report")
                           )

                           ),



                    column(width=9,
                           box(width=4, title=HTML("<b>Pre-Centrifugation Histogramm</b>"), status="warning",
                               plotlyOutput("pre_histo")

                           ),
                           box(width=4, title=HTML("<b>Post-Centrifugation Histogramm</b>"), status="success",
                               plotlyOutput("post_histo")
                           ),
                           box(width=4, title=HTML("<b>2D Hex Histogramm</b>"), status="danger",
                               plotlyOutput("combo_histo")
                           ),
                           box(width=4, title=HTML("<b>Pre-Centrifugation SPREC Distribution</b>"),status="warning",
                               plotOutput("pre_sprec_pie")
                           ),
                           box(width=4, title=HTML("<b>Post-Centrifugation SPREC Distribution</b>"), status="success",
                               plotOutput("post_sprec_pie")
                           )
                           #box(width=4, title="SPREC Distribution",solidHeader=TRUE, status="primary",
                          #     plotlyOutput("pre_sprec_bar")
                           #),
                           #box(width=4, title="SPREC Distribution",solidHeader=TRUE, status="primary",
                          #     plotlyOutput("post_sprec_bar")
                           #)

                           )
                  ),
                  fluidRow(
                    box(width=12,title=HTML("<b>Output</b>"),status="primary",
                        DT::dataTableOutput('batch_datatable')
                    )
                    # box(width=12,title="Output",solidHeader=TRUE,status="primary",
                    #     DT::dataTableOutput('temp')
                    # ),
                    # box(width=12,title="Output",solidHeader=TRUE,status="primary",
                    #     DT::dataTableOutput('temp2')
                    # )
                  )
          ),


# Tools -------------------------------------------------------------------

          tabItem(tabName="blood_tables_helper",
                  box(width=12,status="info",
                  fluidRow(

                  box(
                    width=2,title=HTML("<b>Sample File</b>"),status="primary",
                    fileInput("helper_file", "File input", multiple=FALSE),
                  ),
                  box(status="primary",
                    p("This tool converts time stamp entries in the format of DMY H:M, the delimiter between DMY does not matter. / or : or - will all work"),
                    p("It then calclualted pre-centrifugation, post-centrifugation and total time-to-freeze to be used by the bath processor"),
                    p("Currently it requires the following columnnames: ID, Draw, Centrifugation, Freeze"),
                    p("An example table is included on the GitHub page.")
                  ),
                  box(width=2,title=HTML("<b>Download</b>"),status="primary",
                      downloadButton("download_helper", "Download Template")
                  )
                  ),
                  fluidRow(
                    box(width=2,title=HTML("<b>SPREC</b>"),status="warning",
                    shinyWidgets::radioGroupButtons("sprec","Create SPREC Classification columns?",
                                                    choices = list("Yes" = "Yes", "No" = "No"),
                                                    selected = "Yes")

                    ),
                    box(width=5,title=HTML("<b>Centrifugation</b>"),status="warning",
                        sliderInput("centtime",
                                    "Duration of Centrifugation:",
                                    min = 0,
                                    max = 20,
                                    value = 0,
                                    step=1),
                        p("If the centrifguation time stamp is from the beginning of the centrifugation, the duration of the centrifugation needs to be subtracted for an accurate post-centrifugation SPREC.")
                    )
                  )
                  ),
                    box(title=HTML("<b>Output</b>"),status="primary",
                        width=12,
                    DT::dataTableOutput('helper_datatable')

                  )
          )
      )
    )
)

# Server ------------------------------------------------------------------
server <- function(input, output) {

  # Tables ------------------------------------------------------------------
    pb_plots_pre_percent<-reactive({input$pb_plots_pre_percent})
    pb_plots_pre_cutoff<-reactive({input$pb_plots_pre_cutoff})
    pb_plots_post_percent<-reactive({input$pb_plots_post_percent})
    pb_plots_post_cutoff<-reactive({input$pb_plots_post_cutoff})

  #create reactive table
    pb_table_pre_lollis<-reactive({
      req(input$pb_plots_pre_percent)
      precent_lollis<-LM_pre
      precent_lollis$percent<-precent_lollis$intercept*(pb_plots_pre_percent()/100)
      precent_lollis$limit<-ifelse((precent_lollis$slope < 0), c(precent_lollis$intercept-precent_lollis$percent),c(precent_lollis$intercept+precent_lollis$percent))
      precent_lollis$timepoint<-((precent_lollis$limit-precent_lollis$intercept)/precent_lollis$slope)
      precent_lollis<-precent_lollis %>% mutate(across(where(is.numeric), round, digits=2))
    })

    pb_table_post_lollis<-reactive({
      req(input$pb_plots_post_percent)
      postcent_lollis<-LM_post
      postcent_lollis$percent<-postcent_lollis$intercept*(pb_plots_post_percent()/100)
      postcent_lollis$limit<-ifelse((postcent_lollis$slope < 0), c(postcent_lollis$intercept-postcent_lollis$percent),c(postcent_lollis$intercept+postcent_lollis$percent))
      postcent_lollis$timepoint<-((postcent_lollis$limit-postcent_lollis$intercept)/postcent_lollis$slope)
      postcent_lollis<-postcent_lollis %>% mutate(across(where(is.numeric), round, digits=2))

    })

    pb_table_pre_a1<-reactive({
      precenta1<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint < 0.5)
      precenta1<-subset(precenta1,select =c(name,Type,timepoint))
    })

    pb_table_pre_a<-reactive({
      precenta<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint >= 0.5  & pb_table_pre_lollis()$timepoint <2 )
      precenta<-subset(precenta,select =c(name,Type,timepoint))
    })
    pb_table_pre_c<-reactive({
      precentc<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint >= 2 & pb_table_pre_lollis()$timepoint <4 )
      precentc<-subset(precentc,select =c(name,Type,timepoint))
    })
    pb_table_pre_e<-reactive({
      precente<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint >= 4 & pb_table_pre_lollis()$timepoint <8 )
      precente<-subset(precente,select =c(name,Type,timepoint))
    })

    pb_table_pre_g<-reactive({
      precentg<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint >= 8 & pb_table_pre_lollis()$timepoint <12 )
      precentg<-subset(precentg,select =c(name,Type,timepoint))
    })

    pb_table_pre_i<-reactive({
      precenti<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint >= 12 & pb_table_pre_lollis()$timepoint <24 )
      precenti<-subset(precenti,select =c(name,Type,timepoint))
    })
    pb_table_pre_k<-reactive({
      precentk<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint >= 24 & pb_table_pre_lollis()$timepoint <48 )
      precentk<-subset(precentk,select =c(name,Type,timepoint))
    })

    pb_table_pre_m<-reactive({
      precentm<-subset(pb_table_pre_lollis(),pb_table_pre_lollis()$timepoint >48 )
      precentm<-subset(precentm,select =c(name,Type,timepoint))
    })

    pb_table_post_b<-reactive({
      postcentb<-subset(pb_table_post_lollis(), pb_table_post_lollis()$timepoint < 1)
      postcentb<-subset(postcentb,select =c(name,Type,timepoint))
    })

    pb_table_post_d<-reactive({
      postcentd<-subset(pb_table_post_lollis(), pb_table_post_lollis()$timepoint >= 1 & pb_table_post_lollis()$timepoint <2 )
      postcentd<-subset(postcentd,select =c(name,Type,timepoint))
    })
    pb_table_post_f<-reactive({
      postcentf<-subset(pb_table_post_lollis(), pb_table_post_lollis()$timepoint >= 2 & pb_table_post_lollis()$timepoint <8 )
      postcentf<-subset(postcentf,select =c(name,Type,timepoint))
    })
    pb_table_post_h<-reactive({
      postcenth<-subset(pb_table_post_lollis(), pb_table_post_lollis()$timepoint >= 8 & pb_table_post_lollis()$timepoint <24 )
      postcenth<-subset(postcenth,select =c(name,Type,timepoint))
    })

    pb_table_post_j<-reactive({
      postcentj<-subset(pb_table_post_lollis(), pb_table_post_lollis()$timepoint >= 24 & pb_table_post_lollis()$timepoint <48 )
      postcentj<-subset(postcentj,select =c(name,Type,timepoint))
    })

    pb_table_post_m<-reactive({
      postcentm<-subset(pb_table_post_lollis(), pb_table_post_lollis()$timepoint >= 48)
      postcentm<-subset(postcentm,select =c(name,Type,timepoint))
    })

  #create plot
    pb_plot_pre_lollis<-reactive({
      req(input$pb_plots_pre_cutoff)
      precent_lollis<-pb_table_pre_lollis()
      precent_lollis<-precent_lollis %>% filter(timepoint<pb_plots_pre_cutoff())
      precent_lollis<-precent_lollis %>% mutate(name = fct_reorder(name, timepoint))
      gg_timeline_pre<-ggplot(precent_lollis,aes(x=timepoint,y=Type,color=name))+
        geom_vline(xintercept=0.5,alpha=0.2)+
        geom_vline(xintercept=2,alpha=0.2)+
        geom_vline(xintercept=4,alpha=0.2)+
        geom_vline(xintercept=8,alpha=0.2)+
        # geom_vline(xintercept=12,alpha=0.2)+
        # geom_vline(xintercept=24,alpha=0.2)+
        # geom_vline(xintercept=48,alpha=0.2)+
        geom_hline(yintercept="Serum",alpha=0.6)+
        geom_hline(yintercept="EDTA",alpha=0.6)+
        geom_hline(yintercept="LiHep",alpha=0.6)+
        geom_point(aes(x=timepoint,y=Type),size=6)+
        geom_label(aes(label=name))+
        ylab("")+
        xlab("")+
        theme_classic()+
        theme(legend.position="blank",
              axis.text=element_text(size=10))+
        annotate("text",x=0.13,y=3.3,label="A1",size=6)+
        annotate("text",x=1.2,y=3.3,label="A",size=6)+
        annotate("text",x=3,y=3.3,label="C",size=6)+
        annotate("text",x=6,y=3.3,label="E",size=6)
        # annotate("text",x=10,y=3.52,label="G",size=6)+
        # annotate("text",x=18,y=3.52,label="I",size=6)+
        # annotate("text",x=36,y=3.52,label="K",size=6)
    })

    #create plot
    pb_plot_post_lollis<-reactive({
      req(input$pb_plots_post_cutoff)
      postcent_lollis<-pb_table_post_lollis()
      postcent_lollis<-postcent_lollis %>% filter(timepoint<pb_plots_post_cutoff())
      postcent_lollis<-postcent_lollis %>% mutate(name = fct_reorder(name, timepoint))
      gg_timeline_pre<-ggplot(postcent_lollis,aes(x=timepoint,y=Type,color=name))+
        geom_vline(xintercept=0.5,alpha=0.2)+
        geom_vline(xintercept=2,alpha=0.2)+
        geom_vline(xintercept=4,alpha=0.2)+
        geom_hline(yintercept="Serum",alpha=0.1)+
        geom_hline(yintercept="EDTA",alpha=0.1)+
        geom_hline(yintercept="LiHep",alpha=0.1)+
        geom_point(aes(x=timepoint,y=Type),size=6)+
        geom_label(aes(label=name))+
        ylab("")+
        xlab("")+
        theme_classic()+
        theme(legend.position="blank",
              axis.text=element_text(size=10))+
        annotate("text",x=0.13,y=3.3,label="A1",size=6)+
        annotate("text",x=1.2,y=3.3,label="A",size=6)+
        annotate("text",x=3,y=3.3,label="C",size=6)+
        annotate("text",x=6,y=3.3,label="E",size=6)
    })

  #plot output
     output$pb_lollis_pre<-renderPlotly({
       pb_plot_pre_lollis()
     })

     output$pb_lollis_post<-renderPlotly({
       pb_plot_post_lollis()
     })

     #table output
     output$pb_table_pre_a1<-DT::renderDataTable(
       pb_table_pre_a1(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_pre_a<-DT::renderDataTable(
       pb_table_pre_a(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_pre_c<-DT::renderDataTable(
       pb_table_pre_c(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_pre_e<-DT::renderDataTable(
       pb_table_pre_e(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_pre_g<-DT::renderDataTable(
       pb_table_pre_g(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_pre_i<-DT::renderDataTable(
       pb_table_pre_i(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_pre_k<-DT::renderDataTable(
       pb_table_pre_k(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_pre_m<-DT::renderDataTable(
       pb_table_pre_m(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_post_b<-DT::renderDataTable(
       pb_table_post_b(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_post_d<-DT::renderDataTable(
       pb_table_post_d(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_post_f<-DT::renderDataTable(
       pb_table_post_f(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_post_h<-DT::renderDataTable(
       pb_table_post_h(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_post_j<-DT::renderDataTable(
       pb_table_post_j(),options=list(autoWidth=TRUE),rownames=FALSE
     )
     output$pb_table_post_m<-DT::renderDataTable(
       pb_table_post_m(),options=list(autoWidth=TRUE),rownames=FALSE
     )

     output$pb_lollis_post_table<-DT::renderDataTable(
       pb_table_post_lollis()
     )


     # Reports single ----------------------------------------------------------
  pre<-reactive({input$TTZ})
  post<-reactive({input$TTF})
  error<-reactive({input$error})
  type<-reactive({switch(input$sample_type,"edta"="EDTA","lihep"="LiHep","serum"="Serum")})
  major<-reactive({input$major})
  minor<-reactive({input$minor})
  major_color<-reactive({input$major_color})
  minor_color<-reactive({input$minor_color})
  neutral_color<-reactive({input$neutral_color})


 #setup error data filer
  data_setup_pre_error<-reactive({
    type<-type()
    pre<-pre()
    error<-error()
    if(type=="Serum"){
      df<-LM_pre_error %>% filter(Type=="Serum")
    } else if(type=="EDTA"){
      df<-LM_pre_error %>% filter(Type=="EDTA")
    } else if(type=="LiHep"){
      df<-LM_pre_error %>% filter(Type=="LiHep")
    }
    if(pre<=2){
        df<- df %>% filter(Time==2)
    } else if(pre>2 & pre<=4){
        df<- df %>% filter(Time==4)
    } else if(pre>4 & pre<=6){
        df<- df %>% filter(Time==6)
    } else if (pre>6){
        df<- df %>% filter(Time==8)
    }
    df$percent<-abs(df$percent)
    df<-df %>% filter(percent<error)
    names(df)[names(df) =="ID"]<-"name"
    df<-subset(df,select=name)
      })

  data_setup_post_error<-reactive({
    type<-type()
    post<-post()
    error<-error()
    if(type=="Serum"){
      df<-LM_post_error %>% filter(Type=="Serum")
    } else if(type=="EDTA"){
      df<-LM_post_error %>% filter(Type=="EDTA")
    } else if(type=="LiHep"){
      df<-LM_post_error %>% filter(Type=="LiHep")
    }
    if(post<=2){
      df<- df %>% filter(Time==2)
    } else if(post>2 & post<=4){
      df<- df %>% filter(Time==4)
    } else if(post>4 & post<=6){
      df<- df %>% filter(Time==6)
    } else if (post>6){
      df<- df %>% filter(Time==8)
    }
    df$percent<-abs(df$percent)
    df<-df %>% filter(percent<error)
   names(df)[names(df) =="ID"]<-"name"
   df<-subset(df,select=name)

  })

  #combine_list of pre and post errors
  data_error<-reactive({
    df<-full_join(data_setup_pre_error(),data_setup_post_error(),by="name")
  })

  #filter data
  data_pre<-reactive({
    type<-type()
    data_error<-data_error()
    if(type=="Serum"){
      df<-LM_pre %>% filter(Type=="Serum")
    } else if(type=="EDTA"){
      df<-LM_pre %>% filter(Type=="EDTA")
    } else if(type=="LiHep"){
      df<-LM_pre %>% filter(Type=="LiHep")
    }
    #filter according to errors
    df<-right_join(df,data_error(),by="name")
  })

  #filter data
  data_post<-reactive({
    type<-type()
    if(type=="Serum"){
      df<-LM_post %>% filter(Type=="Serum")
    } else if(type=="EDTA"){
      df<-LM_post %>% filter(Type=="EDTA")
    } else if(type=="LiHep"){
      df<-LM_post %>% filter(Type=="LiHep")
    }
    #filter according to errors
    df<-right_join(df,data_error(),by="name")

  })

  # #output data table
  # output$datatable2<-renderDataTable({
  #   data_error()
  # })
  #
  # output$datatable3<-renderDataTable({
  #   data_setup_pre_error()
  # })
  # output$datatable4<-renderDataTable({
  #   data_setup_post_error()
  # })



#calculate percentages and create colored table
    calc_reactive<-reactive({
      req(data_pre(),data_post())
      final<-calculator(data_pre(),data_post(),pre(),post())
      color_changes<-c(-major(),-minor(),minor(),major())
      colors<-c(major_color(),minor_color(),neutral_color(),minor_color(),major_color())
      final<-datatable(final,
                       extensions = 'Buttons',
                       options = list(dom = 'lfrtpB',
                                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                      lengthMenu = list(c(25, 100, -1), c('10','25','100', 'All')),
                                      pageLength=10
                                  )
                       )  %>%
                           formatStyle("PreCent",  backgroundColor = styleInterval(color_changes,colors),
                                      fontWeight = 'bold') %>%
                           formatStyle("PostCent",  backgroundColor = styleInterval(color_changes,colors),
                                      fontWeight = 'bold') %>%
                           formatStyle("Combined",  backgroundColor = styleInterval(color_changes,colors),
                                      fontWeight = 'bold') %>%
                           formatRound(columns=c("PreCent","PostCent","Combined"),digits=2)
      })

#output data table
output$datatable<-renderDataTable({
  calc_reactive()
 })

#Download R Markdown html
  output$report <- downloadHandler(
    filename = "QC.Report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)

      # Set up parameters to pass to Rmd document
      params <- list(
                      TTZ = input$TTZ,
                      TTF=input$TTF,
                      table=calc_reactive(),
                      major=input$major,
                      minor=input$minor
                      )

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      out<-rmarkdown::render(
                                tempReport,
                                params = params,
                                envir = new.env(parent = globalenv())
      )
      file.rename(out,file)
    }
  )


  # Batch -------------------------------------------------------------------
  batch_file <- reactive({
    inFile1<-input$batch_file
    if (is.null(inFile1)) {
      return("")
    }
    # actually read the file
    read.csv(file = inFile1$datapath)
  })

  b_type<-reactive({switch(input$batch_sample_type,"edta"="EDTA","lihep"="LiHep","serum"="Serum")})
  b_major<-reactive({input$b_major})
  b_minor<-reactive({input$b_minor})
  b_major_color<-reactive({input$b_major_color})
  b_minor_color<-reactive({input$b_minor_color})
  b_neutral_color<-reactive({input$b_neutral_color})
  b_error<-reactive({input$b_error})

  b_data_pre<-reactive({
    b_type<-b_type()
    if(b_type=="Serum"){
      df<-LM_pre %>% filter(Type=="Serum")
    } else if(b_type=="EDTA"){
      df<-LM_pre %>% filter(Type=="EDTA")
    } else if(b_type=="LiHep"){
      df<-LM_pre %>% filter(Type=="LiHep")
    }
  })
  #filter data
  b_data_post<-reactive({
    b_type<-b_type()
    if(b_type=="Serum"){
      df<-LM_post %>% filter(Type=="Serum")
    } else if(b_type=="EDTA"){
      df<-LM_post %>% filter(Type=="EDTA")
    } else if(b_type=="LiHep"){
      df<-LM_post %>% filter(Type=="LiHep")
    }
  })

  batch_combined<-reactive({
    req(input$batch_file)
    b_final<-calculator_multi(batch_file(),b_data_pre(),b_data_post())
  })
  batch_creation<-reactive({
    b_final<-batch_combined()
    color_changes<-c(-b_major(),-b_minor(),b_minor(),b_major())
    colors<-c(b_major_color(),b_minor_color(),b_neutral_color(),b_minor_color(),b_major_color())
    columns<-colnames(b_final[5:ncol(b_final)])
    b_final<-datatable(b_final,columns,
                     extensions = 'Buttons',
                     options = list(dom = 'lfrtpB',
                                    scrollX=TRUE,
                                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                    lengthMenu = list(c(25, 100, -1), c('10','25','100', 'All')),
                                    pageLength=10
                     )
    )  %>%
      formatStyle(columns,  backgroundColor = styleInterval(color_changes,colors),
                  fontWeight = 'bold') %>%
      formatRound(columns,digits=2)
  })


  output$batch_datatable<-renderDataTable({
    batch_creation()
  })

  #filter data
  output$batcher<-renderDataTable({
    b_data_pre()
  })

  pre_histogram<-reactive({
    req(input$batch_file)
    data<-batch_file()
    ggplot(data,aes(x=PreCent,fill=cut(PreCent,50),color=cut(PreCent,50)))+
      geom_histogram(binwidth = 0.5)+
      theme_classic()+
      xlab("Hours")+
      ylab("")+
      theme(legend.position = "none",axis.text =element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      scale_fill_viridis(discrete=TRUE,option="H")+
      scale_color_viridis(discrete=TRUE,option="H")
  })

  output$pre_histo<-renderPlotly(
    pre_histogram()
  )

  post_histogram<-reactive({
    req(input$batch_file)
    data<-batch_file()
    ggplot(data,aes(x=PostCent,fill=cut(PostCent,50),color=cut(PostCent,50)))+
      geom_histogram(binwidth = 0.5)+
      theme_classic()+
      xlab("Hours")+
      ylab("")+
      theme(legend.position = "none",axis.text =element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      scale_fill_viridis(discrete=TRUE,option="H")+
      scale_color_viridis(discrete=TRUE,option="H")
  })


  output$post_histo<-renderPlotly(
    post_histogram()
  )

  combo_histogram<-reactive({
    req(input$batch_file)
    data<-batch_file()
    ggplot(data,aes(x=PreCent,y=PostCent))+
      geom_hex(bins=20)+
      scale_fill_viridis(discrete=FALSE,option="H")+
      theme_classic()+
      xlab("Pre-Cent / Hours")+
      ylab("Post-Cent / Hours")+
      labs(fill="Count")+
      theme(axis.text =element_text(size=16),axis.title=element_text(size=18,face="bold"))
  })


  output$combo_histo<-renderPlotly(
    combo_histogram()
  )

  add_sprec<-reactive({
    req(input$batch_file)
    data<-batch_file()
    data$SPREC_PRE<-ifelse(data$PreCent < 0.5,"A1",
                           ifelse(data$PreCent < 2, "A",
                                  ifelse(data$PreCent < 4, "C",
                                         ifelse(data$PreCent < 8, "E",
                                                ifelse(data$PreCent < 12, "G",
                                                       ifelse(data$PreCent < 24, "I",
                                                              ifelse(data$PreCent < 48,"K","M")
                                                       )
                                                )
                                         )
                                  )
                           )
    )
    data$SPREC_POST<-ifelse(data$PostCent < 1,"B",
                            ifelse(data$PostCent <2, "D",
                                   ifelse(data$PostCent <8,"F",
                                          ifelse(data$PostCent < 24, "H",
                                                 ifelse(data$PostCent < 48, "J","M")
                                          )
                                   )
                            )
    )
    data<-data
  })

  pre_sprec_plot_prepper<-reactive({
    data<-add_sprec()
    df<-data %>% count(SPREC_PRE)
    df$SPREC_PRE<-factor(df$SPREC_PRE,levels=c("A1","A","C","E","G","I","K","M"))
    df<-df
  })

  post_sprec_plot_prepper<-reactive({
    data<-add_sprec()
    df<-data %>% count(SPREC_POST)
    df$SPREC_POST<-factor(df$SPREC_POST,levels=c("B","D","F","H","M"))
    df<-df
  })

  pre_sprec_pie<-reactive({
    df<-pre_sprec_plot_prepper()
    ggplot(df,aes(area=n,fill=SPREC_PRE,label=SPREC_PRE))+
      geom_treemap(color="white")+
      geom_treemap_text(color="white",place="centre",fontface="bold")+
      scale_fill_brewer(palette="Set1")+
      theme(legend.position = "none")

  })

  pre_sprec_bar<-reactive({
    df<-pre_sprec_plot_prepper()
    ggplot(df,aes(x=SPREC_PRE,y=n,fill=SPREC_PRE))+
      geom_bar(stat="identity")+
      theme_classic()+
      theme(legend.position="none",axis.text =element_text(size=16))+
      scale_fill_viridis(discrete=TRUE,option="C")
  })


  post_sprec_pie<-reactive({
    df<-post_sprec_plot_prepper()
    ggplot(df,aes(area=n,fill=SPREC_POST,label=SPREC_POST))+
      geom_treemap(color="white")+
      geom_treemap_text(color="white",place="centre",fontface="bold")+
      scale_fill_brewer(palette="Set1")+
      theme(legend.position = "none")


  })

  post_sprec_bar<-reactive({
    df<-post_sprec_plot_prepper()
    ggplot(df,aes(x=SPREC_POST,y=n,fill=SPREC_POST))+
      geom_bar(stat="identity")+
      theme_classic()+
      xlab("")+
      ylab("")+
    theme(legend.position="none",axis.text =element_text(size=16))+
    scale_fill_viridis(discrete=TRUE,option="C")
  })

  output$pre_sprec_pie<-renderPlot({
    pre_sprec_pie()
  })

  output$pre_sprec_bar<-renderPlotly({
    pre_sprec_bar()
  })

  output$post_sprec_pie<-renderPlot({
    post_sprec_pie()
  })

  output$post_sprec_bar<-renderPlotly({
    post_sprec_bar()
  })

  output$temp<-renderDataTable({
    pre_sprec_plot_prepper()
  })

  output$temp2<-renderDataTable({
    post_sprec_plot_prepper()
  })



    output$batch_report <- downloadHandler(
    filename = "QC_Report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)

      # Set up parameters to pass to Rmd document
      params <- list(
        table=batch_creation(),
        major=input$major,
        minor=input$minor,
        pre_histo=pre_histogram(),
        post_histo=post_histogram(),
        combo_histo=combo_histogram(),
        gg_pre_bar=pre_sprec_bar(),
        gg_post_bar=pre_sprec_bar(),
        gg_pre_pie=pre_sprec_pie(),
        gg_post_pie=post_sprec_pie(),
        pre_sprec=pre_sprec_plot_prepper(),
        post_sprec=post_sprec_plot_prepper(),
        error_input=input$b_error

      )

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      out<-rmarkdown::render(
        tempReport,
        params = params,
        envir = new.env(parent = globalenv())
      )
      file.rename(out,file)
    }
  )




  # Tools -------------------------------------------------------------------
    sprec<-reactive({switch(input$sprec, "Yes"="Yes","No"="No")})
    centtime<-reactive({input$centtime})


    helper_creation<-reactive({
    req(input$helper_file)
    date_helper(helper_file(),sprec(),centtime())
  })

  helper_file <- reactive({
    inFile<-input$helper_file
    if (is.null(inFile)) {
      return("")
    }
    # actually read the file
    read.csv(file = inFile$datapath)
    })

  output$helper_datatable<-renderDataTable({
    helper_creation()
  })
  output$helper<-renderDataTable({
    helper_file()
  })

  output$download_helper <- shiny::downloadHandler(
    filename = function() {
      paste("_dated", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(helper_creation(), file, row.names = FALSE)
    }

  )
}

# Run the application
shinyApp(ui = ui, server = server)
