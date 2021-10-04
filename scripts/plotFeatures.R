library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(plyr)

args <- commandArgs(trailingOnly = TRUE) 
suppressPackageStartupMessages(library(data.table))

fn<-args[1] 
out<-args[2]
features <- read.csv(fn,check.names = FALSE)

metacol<-c("chrom","chromStart","chromEnd","name","strand","hyb","up","footprint","down","neighbourhood","stacking")
features <- features %>% select(-metacol)
features$score<-as.factor(features$score)
col_names <- features %>% select(-score)%>%colnames()

group.colors=c(A="#4b51ac",C="#2d9e40",G="#000000",T="#cd1e25")
score.colors=c('0'="#7f7f7f",'1'="#f8766d")

pdf(out,width = 3.5, height = 2.5)
for (c in col_names){
  if(grepl('template$',c)){
    p<-ggplot(features,aes(get(c),fill=get(c),alpha=score))+
      geom_bar(position='dodge',aes(y = (..count..)*200/sum(..count..)))+
      theme_classic(base_size = 10)+
      scale_alpha_discrete(range=c(0.5, 1))+
      scale_fill_manual(values=group.colors)+
      labs(x='Nucleotide',y='Percentage of sites',alpha='Pausing',fill='',title=c)
    print(p)
  }
  else if(grepl('pct$',c)){
    p<-ggplot(features,aes(get(c)*100,color=score,fill=score))+
      #  geom_histogram(aes(y=..density..),binwidth = 1, position="identity",alpha=0)+
      scale_fill_manual(values=score.colors)+
      scale_color_manual(values=score.colors)+
      geom_density(alpha=0.05)+
      theme_classic(base_size = 10)+
      labs(y='Density',x=c)+
      xlim(0,100)
    print(p)
  }
  else{
    p<-ggplot(features,aes(get(c),color=score,fill=score))+
      #  geom_histogram(aes(y=..density..),binwidth = 1, position="identity",alpha=0)+
      geom_density(alpha=0.5)+
      scale_fill_manual(values=score.colors)+
      scale_color_manual(values=score.colors)+
      theme_classic(base_size = 10)+
      labs(y='Density',title=c)
    print(p)
  }
}
dev.off()

