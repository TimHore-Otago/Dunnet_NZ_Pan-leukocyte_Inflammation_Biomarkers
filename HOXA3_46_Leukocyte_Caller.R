#Activate Packages
library(data.table)
library(dplyr)
library(RColorBrewer)
library(reshape)
library(reshape2)

#set the working directory
setwd("/Path/to/input/directory/")    ### Change directory###
getwd()

## Input files are text files produced by the bismark methylation extractor.
#Inputting data 
input.files <- list.files(pattern = "HOXA3", full.names = FALSE, recursive = FALSE)
input.files

RowNames <- c("Non-leukocyte","Leukocyte","Unknown", "SampleID")
CombinedCallData <- as.data.frame(RowNames)


#Extracting Methylation calls
for (i in 1:length(input.files)) {
  samfile <- read.delim(input.files[i], sep = '', skip = 1, header = F)
  
  colnames(samfile) <- c("read","strand","chr","start","meth")
  
  #make a new data frame where columns represent the CG positions, rows are reads, and "Z" values show methylation
  data_frame = dcast(data = samfile,formula = read~start,value.var = "meth")
  rownames(data_frame) <- data_frame$read
  data_frame$read <- c()  
  data_frame <- data_frame[complete.cases(data_frame),]
  
  #MAP4K1 Specific Changes!
  data_frame$`28` = NULL

  levels=unique(do.call(c,data_frame))
  MethPercent <- as.data.frame(sapply(levels,function(x)rowSums(data_frame==x)))
  
  totalCpG <- rowSums(MethPercent[1,])  
  MethPercent$proportion <- MethPercent$Z/totalCpG
  #Sorting Methylated reads into groups of high/med/low methylation per read; can be adjusted to suit.
  MethylationCat <- function(MethCat) {
    if (MethCat >= 6/10) { return("PBMC") }
    if (MethCat > 3/10) {return("Unknown") }
    else if (MethCat <= 3/10) {return("Lung_Organoid")}
  }
  
  MethPercent$Call <- lapply(MethPercent$proportion, FUN = MethylationCat)
  
  Dist <- as.data.frame(MethPercent$Call)
  Dist <- as.data.table(t(as.matrix(Dist)))
  Dist <- as.data.frame(table(Dist$V1))
  colnames(Dist) <- c("Catagory", "Sum")
  Percentage <- round(Dist$Sum/sum(Dist$Sum)*100, digits = 1)
  Lables <- paste(Dist$Cat, Percentage)
  Lables <- paste(Lables, "%", sep = "")
  
  Output_stats <- as.data.frame(cbind(Dist, Percentage))
  SampleID <- as.character(input.files[i])
  Output_stats[4,] <- SampleID
  
  CombinedCallData <- cbind(CombinedCallData, Output_stats)
}

CombinedCallData <- CombinedCallData[ ,-which(names(CombinedCallData) %in% c("Catagory", "Sum"))]
CombinedCallData <- transpose(CombinedCallData)
CombinedCallData[is.na(CombinedCallData)] <- 0

write.table(CombinedCallData, paste("HOXA3-Combined_calls.txt", sep = ""), 
            sep = "\t", row.names = F, col.names = F, quote = F)
