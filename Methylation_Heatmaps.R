#Activate Packages
library(data.table)
library(dplyr)
library(RColorBrewer)
library(reshape)
library(reshape2)

## Input files are text files produced by the bismark methylation extractor.
#set the working directory
setwd("/Path/to/input/directory/")    ### Change directory###
getwd()
#Inputting data 
input.files <- list.files(pattern = "CpG", full.names = FALSE, recursive = FALSE)
input.files

#Extracting Methylation calls
for (i in 1:length(input.files)) {
  
  samfile <- read.delim(input.files[i], sep = '', skip = 1, header = F)
  
  colnames(samfile) <- c("read","strand","chr","start","meth")
  
  #make a new data frame where columns represent the CG positions, rows are reads, and "Z" values show methylation
  data_frame = dcast(data = samfile,formula = read~start,value.var = "meth")
  data_frame
  #remove read names
  rownames(data_frame) <- data_frame$read
  data_frame$read <- c()
  
  #convert "Z" into numbers
  data_frame[data_frame=="Z"] <- as.numeric(1)
  data_frame[data_frame=="z"] <- as.numeric(0)
  data_frame
  #prepare data for heatmap (i.e. make numeric, remove rows with NA)
  data_matrix <- as.data.frame(sapply(data_frame, as.numeric))
  data_matrix <- data_matrix[complete.cases(data_matrix),]
  data_matrix
  #Draw heatmap
  hmcol<-rev(brewer.pal(5,"RdBu"))
  pdf.name=paste("/Path/to/output/directory/","Heatmap_",input.files[i],".pdf", sep = '')    ### Change directory###
  
  #calculate mean to 1dp. Change the nrow multiplication to the number of CpGs in the read.
  mean.metylation <- sum(data_matrix)/(ncol(data_matrix)*nrow(data_matrix))*100
  round.mean.meth <- round(mean.metylation, 1)
  
  pdf(pdf.name)
  heatmap(as.matrix(data_matrix), Colv = NA, Rowv = TRUE, main = input.files[i], col=hmcol, scale = "none", xlab="CpG Position",
          ylab=paste("Reads = ",nrow(data_matrix),"    ",round.mean.meth," % Methylation", sep = ""))
  dev.off()
}
print("Complete")

