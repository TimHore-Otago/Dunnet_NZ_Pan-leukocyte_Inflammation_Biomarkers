library(data.table)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(reshape)
library(reshape2)
library(pROC)
library(randomForest)

## Input files are text files produced by the bismark methylation extractor.
# Input condition one data
samfile <- read.delim("~/path/to/input/file1.txt", sep = '', skip = 1, header = F)
colnames(samfile) <- c("read","strand","chr","start","meth")
data_frame = dcast(data = samfile,formula = read~start,value.var = "meth")
rownames(data_frame) <- data_frame$read
data_frame$read <- c()  
data_frame <- data_frame[complete.cases(data_frame),]
data_frame$MethAmount <- rowSums(data_frame == "Z")
meth <- as.data.frame(data_frame$MethAmount)
colnames(meth) <- c("Methylaed CpGs")
meth$condition <- 0

# Input condition two data
samfile2 <- read.delim("~/path/to/input/file2.txt", sep = '', skip = 1, header = F)
colnames(samfile2) <- c("read","strand","chr","start","meth")
data_frame2 = dcast(data = samfile2,formula = read~start,value.var = "meth")
rownames(data_frame2) <- data_frame2$read
data_frame2$read <- c()  
data_frame2 <- data_frame2[complete.cases(data_frame2),]
data_frame2$MethAmount <- rowSums(data_frame2 == "Z")
meth2 <- as.data.frame(data_frame2$MethAmount)
colnames(meth2) <- c("Methylaed CpGs")
meth2$condition <- 1

combined_data <- rbind(meth, meth2)


#order dataframe
df <- combined_data[order(combined_data$`Methylaed CpGs`),]
#variables as numerical
`Methylated CpGs` <- as.numeric(df$`Methylaed CpGs`)
condition <- as.numeric(df$condition)
par(mar = c(6, 5, 4, 2))
#plot variables and produce logistic regression curve
plot(x=jitter(`Methylated CpGs`, amount = 0.1), 
     y = jitter(condition, amount = 0.0),       
     pch = 1,                                   
     cex = 1.5,
     col = adjustcolor("black", alpha.f = 1.0),
     xlab = "Methylated CpGs",                
     ylab = "Probablity of of being leukocyte derived",
     cex.lab = 1.5,
     cex.axis = 1.5
)
glm.fit = glm(condition ~ `Methylated CpGs`, family = binomial)
lines(`Methylated CpGs`, glm.fit$fitted.values)
#plot ROC curve
par(pty="s")
roc(condition, glm.fit$fitted.values, plot=T, legacy.axes=T, lwd=2, print.auc=T, print.auc.x=0.4, col="#984EA3")
roc.info <- roc(condition, glm.fit$fitted.values, plot=T, legacy.axes=T, lwd=2, print.auc=T, print.auc.x=0.4, col="#984EA3")
grid(col = "black")
#Produce threshold data
roc.df <- data.frame(
  TPR=roc.info$sensitivities,
  FPR=1-roc.info$specificities,
  thresholds=roc.info$thresholds)
  roc.df <- roc.df %>% mutate(result = TPR - FPR)
roc.df

#scatterplot version of ROC curve
plot(roc.df$FPR, roc.df$TPR)
lines(roc.df$FPR, roc.df$TPR)
write.table(roc.df, paste("Threshold_data.txt", sep = '/t'), row.names = F, col.names = T, quote = F)
