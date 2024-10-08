calculator<-function(t1,t2){

#Model pre
LM_pre<-read.csv("D:/Code/Projects/Main/QC-Tool/LME_Pre_All.csv",stringsAsFactors = FALSE)
LM_post<-read.csv("D:/Code/Projects/Main/QC-Tool/LME_Post_All.csv",stringsAsFactors = FALSE)
Type<-"EDTA"

if(Type=="Serum"){
  LM_edta_pre<-LM_pre %>% filter(Type=="Serum")
  LM_edta_post<-LM_post %>% filter(Type=="Serum")
} else if(Type=="EDTA"){
  LM_edta_pre<-LM_pre %>% filter(Type=="EDTA")
  LM_edta_post<-LM_post %>% filter(Type=="EDTA")
} else if(Type=="LiHep"){
  LM_edta_post<-LM_post %>% filter(Type=="LiHep")
}

LM_edta_pre_reduced<-data.frame(LM_edta_pre$name,LM_edta_pre$intercept,LM_edta_pre$slope)
LM_edta_pre_reduced<-setNames(LM_edta_pre_reduced,c("Metabolite","PreInt","PreSlope"))
LM_edta_post_reduced<-data.frame(LM_edta_post$name,LM_edta_post$intercept,LM_edta_post$slope)
LM_edta_post_reduced<-setNames(LM_edta_post_reduced,c("Metabolite","PostInt","PostSlope"))
LM_edta_combo<-inner_join(LM_edta_pre_reduced,LM_edta_post_reduced,by="Metabolite")


#input times
t1<-c(0.5)
t2<-c(1)

#calulate overall change
LM_edta_combo$x2<-t1*(LM_edta_combo$PreSlope)+LM_edta_combo$PreInt

LM_edta_combo$x3<-t2*(LM_edta_combo$PostSlope)+LM_edta_combo$x2
LM_edta_combo$delta1<-(LM_edta_combo$x2)-(LM_edta_combo$PreInt)
LM_edta_combo$delta1_percent<-((LM_edta_combo$delta1))*100/(LM_edta_combo$PreInt)

LM_edta_combo$delta2<-((LM_edta_combo$x3)-(LM_edta_combo$x2))
LM_edta_combo$delta2_percent<-((LM_edta_combo$delta2))*100/(LM_edta_combo$x2)
LM_edta_combo$Combined_percent<-LM_edta_combo$delta1_percent+LM_edta_combo$delta2_percent

write.csv(LM_edta_combo,"test.csv")

datatable(LM_edta_combo) %>% formatStyle("delta1_percent",  backgroundColor = styleInterval(c(15, 30), c('green', 'orange', 'red')),
                                         fontWeight = 'bold') %>%
                             formatStyle("delta2_percent",  backgroundColor = styleInterval(c(15, 30), c('green', 'orange', 'red')),
                                         fontWeight = 'bold') %>%
                             formatStyle("Combined_percent",  backgroundColor = styleInterval(c(15, 30), c('green', 'orange', 'red')),
                                        fontWeight = 'bold')
}
