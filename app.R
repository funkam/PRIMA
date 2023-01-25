#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(dplyr)
library(DT)


# Scripts -----------------------------------------------------------------
calculator<-function(t1,t2){
  #Model pre
  LM_edta_pre<-read.csv("A:/Projekte/2022/Analysis/PreAnalytik/Combo/no24/LM/Mixedmodels/LME_Pre_EDTA.csv",stringsAsFactors = FALSE)
  LM_edta_post<-read.csv("A:/Projekte/2022/Analysis/PreAnalytik/Combo/no24/LM/Mixedmodels/LME_Post_EDTA.csv",stringsAsFactors = FALSE)

  LM_edta_pre_reduced<-data.frame(LM_edta_pre$name,LM_edta_pre$intercept,LM_edta_pre$slope)
  LM_edta_pre_reduced<-setNames(LM_edta_pre_reduced,c("Metabolite","PreInt","PreSlope"))
  LM_edta_post_reduced<-data.frame(LM_edta_post$name,LM_edta_post$intercept,LM_edta_post$slope)
  LM_edta_post_reduced<-setNames(LM_edta_post_reduced,c("Metabolite","PostInt","PostSlope"))
  LM_edta_combo<-inner_join(LM_edta_pre_reduced,LM_edta_post_reduced,by="Metabolite")

  #calulate overall change
  LM_edta_combo$x2<-t1*(LM_edta_combo$PreSlope)+LM_edta_combo$PreInt

  LM_edta_combo$x3<-t2*(LM_edta_combo$PostSlope)+LM_edta_combo$x2
  LM_edta_combo$delta1<-(LM_edta_combo$x2)-(LM_edta_combo$PreInt)
  LM_edta_combo$delta1_percent<-((LM_edta_combo$delta1))*100/(LM_edta_combo$PreInt)

  LM_edta_combo$delta2<-((LM_edta_combo$x3)-(LM_edta_combo$x2))
  LM_edta_combo$delta2_percent<-((LM_edta_combo$delta2))*100/(LM_edta_combo$x2)
  LM_edta_combo$Combined_percent<-LM_edta_combo$delta1_percent+LM_edta_combo$delta2_percent
  final<-data.frame(LM_edta_combo$Metabolite,LM_edta_combo$delta1_percent,LM_edta_combo$delta2_percent,LM_edta_combo$Combined_percent)
  final<-setNames(final,c("Metabolite","Pre-Cent Change","Post-Cent Change","Combined_Change"))
}



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("NMR QC Panel"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
        sliderInput("TTZ",
                        "Pre-Centrifugation Time:",
                        min = 0,
                        max = 8,
                        value = 1,
                        step=0.1
                    ),

        sliderInput("TTF",
                    "Post Centifugation Time:",
                    min = 0,
                    max = 8,
                    value = 1,
                    step=0.1)

    ),
    mainPanel(
          div(style = 'overflow-x: scroll',  DT::dataTableOutput('datatable'))

        )
    )
    )


# Define server logic required to draw a histogram
server <- function(input, output) {
  pre<-reactive({input$TTZ})
  post<-reactive({input$TTF})


  calculatro_reactive<-reactive({
    final<-calculator(pre(),post())
  })

  output$datatable<-renderDataTable({
    testy<-datatable(
    calculatro_reactive(),

  })
}

# Run the application
shinyApp(ui = ui, server = server)
