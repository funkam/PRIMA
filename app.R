# NMR QC Tool
#
#
#

library(shiny)
library(shinydashboard)
library(dplyr)
library(DT)
library(shinyWidgets)
#data input
LM_pre<-read.csv("LME_Pre_All.csv",stringsAsFactors = FALSE)
LM_post<-read.csv("LME_Post_All.csv",stringsAsFactors = FALSE)

# Scripts -----------------------------------------------------------------
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
ui <- fluidPage(


    titlePanel("NMR QC Panel"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
        radioGroupButtons("sample_type","Select sample type",choices=c("EDTA"="edta","Lithium-Heparin"="lihep","Serum"="serum"),selected="edta",individual=TRUE),
        sliderInput("TTZ",
                        "Pre-Centrifugation Time:",
                        min = 0,
                        max = 8,
                        value = 0,
                        step=0.1
                    ),

        sliderInput("TTF",
                    "Post Centifugation Time:",
                    min = 0,
                    max = 8,
                    value = 0,
                    step=0.1)

    ),
    mainPanel(
          DT::dataTableOutput('datatable')


        )
    )
    )


# Server ------------------------------------------------------------------
server <- function(input, output) {
  pre<-reactive({input$TTZ})
  post<-reactive({input$TTF})
  type<-reactive({switch(input$sample_type,"edta"="EDTA","lihep"="LiHep","serum"="Serum")})

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



    calc_reactive<-reactive({
      req(data_pre(),data_post())
    final<-calculator(data_pre(),data_post(),pre(),post())
  })


#output data table and color 15% and 30%
  output$datatable<-renderDataTable({
    datatable(calc_reactive())  %>% formatStyle("PreCent",  backgroundColor = styleInterval(c(15, 30), c('', 'orange', 'red')),
                                fontWeight = 'bold') %>%
                                formatStyle("PostCent",  backgroundColor = styleInterval(c(15, 30), c('', 'orange', 'red')),
                                fontWeight = 'bold') %>%
                                formatStyle("Combined",  backgroundColor = styleInterval(c(15, 30), c('', 'orange', 'red')),
                                fontWeight = 'bold') %>%
                                formatRound(columns=c("PreCent","PostCent","Combined"),digits=2)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
