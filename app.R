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

#data input
LM_pre<-read.csv("Precent_LM_all_final.csv",stringsAsFactors = FALSE)
LM_post<-read.csv("Postcent_LM_all_final.csv",stringsAsFactors = FALSE)

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


date_helper<-function(df){
  #create time tables from csv using lubridate
  df$Centrifugation<-as.POSIXct(df$Centrifugation,format="%d.%m.%Y %H:%M")
  df$Freeze<-as.POSIXct(df$Freeze,format="%d.%m.%Y %H:%M")
  df$Draw<-as.POSIXct(df$Draw,format="%d.%m.%Y %H:%M")
  df$PreCent<-as.numeric(difftime(df$Centrifugation,df$Draw),units="hours")
  df$PostCent<-as.numeric(difftime(df$Freeze,df$Centrifugation),units="hours")
  df$TTF<-as.numeric(difftime(df$Freeze,df$Draw),units="hours")
  df<-df
}



# UI ----------------------------------------------------------------------

# Menu --------------------------------------------------------------------


ui <- dashboardPage(
    dashboardHeader(title="NMR QC"),
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
                           menuSubItem("Batch Table Helper",tabName="blood_tables_helper",icon=icon("clock"))
                           )

      )
    ),

# Home --------------------------------------------------------------------


    dashboardBody(
      tabItems(
          tabItem(tabName="home",
                  fluidRow(
                    box(width=6,
                        h1("NMR Metabolomics Quality Control Panel"),
                        #h3(HTML("<b>NMR</b>")),
                        p("A tool for investigating the quality and stability of metabolites in peripheral blood."),
                        p("Developed by",a("Alexander Funk",href="https://www.uniklinikum-dresden.de/de/das-klinikum/kliniken-polikliniken-institute/klinische-chemie-und-laboratoriumsmedizin/forschung/copy_of_EMS")),
                        p("Institute for Clinical Chemistry and Laboratory Medicine"),
                        p("University Hospital Dresden, Fetscherstr. 74, 01307 Dresden")
                        ),
                     box(width=2,title="Links",solidHeader=TRUE,status="primary",
                             p(a("GitHub",href="https://github.com/funkam/QC-Tool")),
                             p(a("shinyapps.io",href="https://funkam.shinyapps.io/QC-Tool/")),
                             #p(a("Publication",href=""))
                             )
                      ),
                      # column(width=6,
                      #        br(),
                      #        img(src='hexlogo.png',width="92%")
                      # )


                  fluidRow(
                    box(title="Description",solidHeader=TRUE,status="primary",
                      p("Linear mixes models were created to show the change of multiple metabolites over time."),
                      p("The data is then presented as change (in %) to its original value for pre-centrifugation times."),
                      p("For post-centrifugation times the change is calculated as an additinal effect on the value already altered by the pre-centrifugation time."),
                      p("The Panel is split into producing plots, as well as color-coded tables.")
                    )
                  ),
                  fluidRow(
                    box(title="Plots",solidHeader=TRUE,status="primary",
                        p("In the 'Plots' tab the user can enter a percentage limit for metabolite change (default 30%), as well as a end time point for the plot (in hours).",br(),br(),"Lollipop plots are then generated dynamically with the given input, giving an overview of the different metabolites and the calucalted time for the given change",br(),br(),"Datapoints are shown for each of the different coagulation tubes.")
                        ),
                    box(title="Tables",solidHeader=TRUE,status="primary",
                        p("The 'Tables' tab allows the user to set a pre-centrifugation time and a post-centrifugation using a Slider.",br(),br(),"A table is then generated that higlights a minor and a major change for each metabolite in that given timeframe",br(),br(),"The colors, as well as the % threshholds, can be adjusted.",br(),br(),"The table can be downloaded as a .HTML report or the table directly as .csv (or other formats).")
                        )
                  )
          ),

# Pre ---------------------------------------------------------------------


          tabItem(tabName="blood_data_pre",
                  fluidRow(
                    column(width=3,
                      box(width=10, title="Pre-Centrifugation Input",solidHeader=TRUE, status="primary",
                                  numericInput("pb_plots_pre_percent","Enter threshhold (%)",value=20),
                                  numericInput("pb_plots_pre_cutoff","Enter plot cutoff (h)",value=10)
                        )
                    ),
                    column(width=3,
                        box(width=12,title="Info",solidHeader=TRUE,status="danger",
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
                               div(style="height:1000px",plotlyOutput(width="60%", height="1000px","pb_lollis_pre"))

                      )
                    )
                  )

                  ),


# Post --------------------------------------------------------------------


          tabItem(tabName="blood_data_post",
                  fluidRow(
                    column(width=3,
                           box(width=10, title="post-Centrifugation Input",solidHeader=TRUE, status="primary",
                               numericInput("pb_plots_post_percent","Enter threshhold (%)",value=20),
                               numericInput("pb_plots_post_cutoff","Enter plot cutoff (h)",value=10)
                           )
                    ),
                    column(width=3,
                           box(width=12,title="Info",solidHeader=TRUE,status="danger",
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
                                    div(style="height:1000px",plotlyOutput(width="60%", height="1000px","pb_lollis_post"))

                           )
                    )
                  )

          ),

# Single Table ------------------------------------------------------------



          tabItem(tabName="blood_tables",
                  fluidRow(
                          box(
                              width=2,title="Sample Type",solidHeader=TRUE,status="primary",
                                column(
                                    width=12,align="center",
                                    radioGroupButtons("sample_type","Select sample type:",choices=c("EDTA"="edta","Lithium-Heparin"="lihep","Serum"="serum"),selected="edta",direction="vertical")
                                )
                          ),
                          box(
                              width=5,title="Style",solidHeader=TRUE,status="primary",
                              fluidRow(
                                column(width=4,
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
                                column(width=4,
                                       box(width=12,
                                           colourInput("neutral_color","Neutral Color:",value="")
                                       )
                                )
                              )
                          ),
                          box(
                            width=3,title="Download Output",solidHeader=TRUE,status="primary",
                            p("Download a report as interactive HTML file.",align="justify"),
                            downloadButton("report", "Generate HTML report")
                            )
                  ),
                  fluidRow(
                          box(width=10,title="Centrifugation Times",solidHeader=TRUE,status="primary",
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
                  fluidRow(
                          box(width=12,
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
                    box(width=4,title="Sample Type",solidHeader=TRUE,status="primary",
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
                      width=4,title="Style",solidHeader=TRUE,status="primary",
                      fluidRow(
                        column(width=4,
                               box(width=12,
                                   numericInput("b_minor","Minor color %:",value=10),
                                   colourInput("b_minor_color","Minor Color:",value="orange"))
                        ),
                        column(width=4,
                               box(width=12,
                                   numericInput("b_major","Major color %:",value=20),
                                   colourInput("b_major_color","Major Color:",value="red")
                               )
                        ),
                        column(width=4,
                               box(width=12,
                                   colourInput("b_neutral_color","Neutral Color:",value="")
                               )
                        )
                      )
                    ),
                    box(
                      width=3,title="Download Output",solidHeader=TRUE,status="primary",
                      p("Download a report as interactive HTML file.",align="justify"),
                      downloadButton("batch_report", "Generate HTML report")
                    )
                  ),
                  fluidRow(
                    box(width=12,title="Output",solidHeader=TRUE,status="primary",
                        DT::dataTableOutput('batch_datatable')
                    )
                  )
          ),


# Tools -------------------------------------------------------------------

          tabItem(tabName="blood_tables_helper",
                  fluidRow(
                  box(
                    width=2,title="Sample File",solidHeader=TRUE,status="primary",
                    fileInput("helper_file", "File input", multiple=FALSE),
                  ),
                  box(
                    p("This tool converts time stamp entries in the format of d.m.Y H:M to a format recognized by R"),
                    p("It then calclualted pre-centrifugation, post-centrifugation and total time-to-freeze to be used by the bath processor"),
                    p("Currently it requires the following columnnames: ID, Draw, Centrifugation, Freeze"),
                    p("An example table is included on the GitHub page.")
                  ),
                  box(width=2,title="Download",solidHeader=TRUE,status="primary",
                      downloadButton("download_helper", "Download Template")
                  )
                  ),
                  fluidRow(
                    box(title="Output",solidHeader=TRUE,status="primary",
                    DT::dataTableOutput('helper_datatable')
                    )
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
      precent_lollis$timepoint<-(precent_lollis$limit-precent_lollis$intercept)/precent_lollis$slope
      precent_lollis<-precent_lollis %>% mutate(across(where(is.numeric), round, digits=2))
    })

    pb_table_post_lollis<-reactive({
      req(input$pb_plots_post_percent)
      postcent_lollis<-LM_post
      postcent_lollis$percent<-postcent_lollis$intercept*(pb_plots_post_percent()/100)
      postcent_lollis$limit<-ifelse((postcent_lollis$slope < 0), c(postcent_lollis$intercept-postcent_lollis$percent),c(postcent_lollis$intercept+postcent_lollis$percent))
      postcent_lollis$timepoint<-(postcent_lollis$limit-postcent_lollis$intercept)/postcent_lollis$slope
      postcent_lollis<-postcent_lollis %>% mutate(across(where(is.numeric), round, digits=2))

    })


    pb_table_pre_a1<-reactive({
      precenta1<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint < 0.5)
      precenta1<-subset(precenta1,select =c(name,Type,timepoint))
    })

    pb_table_pre_a<-reactive({
      precenta<-subset(pb_table_pre_lollis(), pb_table_pre_lollis()$timepoint >= 0.5 & pb_table_pre_lollis()$timepoint <2 )
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
      pre_lollis_plot<-ggplot(precent_lollis,aes(x=name,y=timepoint,color=Type))+
        geom_segment(aes(x=name,xend=name,y=0,yend=timepoint),color="black")+
        geom_point(aes(fill=Type,size=2,shape=Type)) +
        scale_y_continuous(breaks = round(seq(0, max(precent_lollis$timepoint), by = 2),1)) +
        xlab("")+
        ylab("Stability / hours")+
        theme_minimal()+
        theme(legend.position="blank")+
        coord_flip()
    })
    #create plot
    pb_plot_post_lollis<-reactive({
      req(input$pb_plots_post_cutoff)
      postcent_lollis<-pb_table_post_lollis()
      postcent_lollis<-postcent_lollis %>% filter(timepoint<pb_plots_post_cutoff())
      postcent_lollis<-postcent_lollis %>% mutate(name = fct_reorder(name, timepoint))
      post_lollis_plot<-ggplot(postcent_lollis,aes(x=name,y=timepoint,color=Type))+
        geom_segment(aes(x=name,xend=name,y=0,yend=timepoint),color="black")+
        geom_point(aes(fill=Type,size=2,shape=Type)) +
        scale_y_continuous(breaks = round(seq(0, max(postcent_lollis$timepoint), by = 2),1)) +
        xlab("")+
        ylab("Stability / hours")+
        theme_minimal()+
        theme(legend.position="blank")+
        coord_flip()
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
  type<-reactive({switch(input$sample_type,"edta"="EDTA","lihep"="LiHep","serum"="Serum")})
  major<-reactive({input$major})
  minor<-reactive({input$minor})
  major_color<-reactive({input$major_color})
  minor_color<-reactive({input$minor_color})
  neutral_color<-reactive({input$neutral_color})



#filter data
  data_pre<-reactive({
    type<-type()
    if(type=="Serum"){
    df<-LM_pre %>% filter(Type=="Serum")
    } else if(type=="EDTA"){
      df<-LM_pre %>% filter(Type=="EDTA")
      } else if(type=="LiHep"){
        df<-LM_pre %>% filter(Type=="LiHep")
        }
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
  })

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

  batch_creation<-reactive({
    req(input$batch_file)
    b_final<-calculator_multi(batch_file(),b_data_pre(),b_data_post())
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




  # Tools -------------------------------------------------------------------
  helper_creation<-reactive({
    req(input$helper_file)
    date_helper(helper_file())
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
