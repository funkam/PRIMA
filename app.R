# NMR QC Tool

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

# UI ----------------------------------------------------------------------
ui <- dashboardPage(
    dashboardHeader(title="NMR QC"),
    dashboardSidebar(
      sidebarMenu(
                  collapsed=FALSE,
                  menuItem("Home", tabName = "home", icon = icon("home"),selected=T),
                  menuItem("peripheral Blood",icon=icon("syringe"),
                           menuSubItem("Plots",tabName="blood_plots",icon=icon(("chart-bar"))),
                           menuSubItem("Tables",tabName="blood_tables",icon=icon(("table")))
                           )
      )
    ),

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

          tabItem(tabName="blood_plots",
                  fluidRow(
                    column(width=6,
                      box(width=4, title="Pre-Centrifugation Input",solidHeader=TRUE, status="primary",
                                  numericInput("pb_plots_pre_percent","Enter threshhold (%)",value=30),
                                  numericInput("pb_plots_pre_cutoff","Enter plot cutoff (h)",value=10)
                        )
                    ),
                    column(width=6,
                           box(width=4, title="Post-Centrifugation Input",solidHeader=TRUE, status="primary",
                               numericInput("pb_plots_post_percent","Enter threshhold (%)",value=30),
                               numericInput("pb_plots_post_cutoff","Enter plot cutoff (h)",value=10)
                           )
                    ),
                    ),
                        fluidRow(
                          column(width=6,
                            box(width=12,title="Pre-Centrifugation",solidHeader=TRUE,status="primary",
                                div(style="height:1000px",plotlyOutput(width="100%", height="1000px","pb_lollis_pre"))
                            )
                          ),
                          column(width=6,
                                 box(width=12,title="Post-Centrifugation",solidHeader=TRUE,status="primary",
                                     div(style="height:1000px",plotlyOutput(width="100%", height="1000px","pb_lollis_post"))
                                 )
                          )
                        )
                  ),


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
                                           numericInput("minor","Minor color %:",value=15),
                                           colourInput("minor_color","Minor Color:",value="orange"))
                                ),
                                column(width=4,
                                       box(width=12,
                                           numericInput("major","Major color %:",value=30),
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
          tabItem(tabName="blood_plots",
                  sliderInput("input_time_ttz", "Adjust highest pre-centrifugation time shown:",
                              min=0,
                              max=24,
                              value=8,
                              step=0.1),
                  sliderInput("input_time_ttf", "Adjust highest pre-centrifugation time shown:",
                              min=0,
                              max=24,
                              value=8,
                              step=0.1),
                  plotlyOutput("pre_lolli")
                  )
      )
    )
)

# Server ------------------------------------------------------------------
server <- function(input, output) {

####Blood PLots
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
      precent_lollis<-precent_lollis
    })

    pb_table_post_lollis<-reactive({
      req(input$pb_plots_post_percent)
      postcent_lollis<-LM_post
      postcent_lollis$percent<-postcent_lollis$intercept*(pb_plots_post_percent()/100)
      postcent_lollis$limit<-ifelse((postcent_lollis$slope < 0), c(postcent_lollis$intercept-postcent_lollis$percent),c(postcent_lollis$intercept+postcent_lollis$percent))
      postcent_lollis$timepoint<-(postcent_lollis$limit-postcent_lollis$intercept)/postcent_lollis$slope
      postcent_lollis<-postcent_lollis
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


####Blood Tables
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
}

# Run the application
shinyApp(ui = ui, server = server)
