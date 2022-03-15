# install packages
install.packages("dplyr")
install.packages("treemap")
install.packages("sunburstR")
install.packages("readr")
install.packages("tidyverse")
install.packages("data.table")

# libraries
library(tidyverse)
library(dplyr)
library(treemap)
library(sunburstR)
library(purrr)
library(readr)
library(data.table)

#read in the files from the working path
multmerge = function(path){
  filenames=list.files(path=path, full.names=TRUE)
  rbindlist(lapply(filenames, fread))
}

path <- "path/to/tsvs"
data <- multmerge(path)

#clean the data
names(data)<-str_replace_all(names(data), c(" " = "_"))

data<- data %>% 
  select(Contig_id,Element_type,Class, Subclass,Gene_symbol)
data$Class <- sub("^$", "N", data$Class)
data$Subclass <- sub("^$", "N", data$Subclass)
#View(data)

colnames(data) <- c("Contig","Element_type","Class", "Subclass", "Gene")
data$Contig<-gsub("-","_",as.character(data$Contig))
data$Element_type<-gsub("-","_",as.character(data$Element_type))
data$Class<-gsub("-","_",as.character(data$Class))
data$Subclass<-gsub("-","_",as.character(data$Subclass))
data$Gene<-gsub("-","_",as.character(data$Gene))

#calculate percentages
table <- round(100*prop.table(table(data$Gene)))
table <- as.data.frame(table)
table <- table %>% 
  select(Var1,Freq)
colnames(table) <- c("Gene","Frequency_percent")
#View(table)
new_df <- left_join(data,table)
write.csv(new_df,"mulberry_output.csv",row.names = FALSE)

# Reformat data for the sunburstR package
data <- new_df %>%
  filter(Class != "") %>%
  mutate(path = paste(Class, Subclass,Gene, sep="-")) %>%
  dplyr::select(path, Frequency_percent)

# Plot
p <- sunburst(data, legend=FALSE)
print(p)

