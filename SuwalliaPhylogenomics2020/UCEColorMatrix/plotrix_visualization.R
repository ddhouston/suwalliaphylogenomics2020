library(ape)
library(plyr)
library(stringr)
library(plotrix)

# Read in all nexus files and concatenate into one data frame
# In GitHub, these files are under SuwalliaLineageTree_IQTree/iqtrees/trees_50p.zip
file_list <- list.files(path="C:/Users/hmcar/Documents/Research/For Derek/Suwallia/nexusfiles")
datasets <- lapply(file_list, read.nexus.data)
dataset <- data.frame(colnames(paste(file_list)), datasets)
names(datasets) <- paste(file_list)

# Build a skeleton to attach info to
skeleton <- data.frame(no=c(1:66), dummycol=c(1:66))
skeleton$dummycol <- sprintf("%02d", as.numeric(skeleton$dummycol))
skeleton$dummycol <- paste("Su", skeleton$dummycol, sep="")

skeleton2 <- skeleton
for (i in 1:length(datasets)) {
 x <- as.data.frame(names(datasets[[i]]))
 x$no <- x$`names(datasets[[i]])`  
 x$no <- str_remove(x$no, "Su")
 x$no <- as.integer(x$no) 
 skeleton2 <- merge(skeleton2, x, by="no", all=TRUE)
 rm(x)
}

#Get rid of the columns we don't need
skeleton3 <- skeleton2[,3:298]
colnames(skeleton3) <- str_remove(file_list, ".nexus") # Name our columns

# Code by presence/absence of data
skeleton3[!is.na(skeleton3) ] <- 1
skeleton3[is.na(skeleton3) ] <- 0
skeleton3  <- sapply(skeleton3, as.numeric )

# Create our final data frame
finaldata <- as.data.frame(skeleton3)
finaldata <- finaldata[apply(finaldata[,-1], 1, function(x) !all(x==0)),]
rownames(finaldata) <- skeleton$dummycol

# Define the colors to match the tree
cellcolors <- unlist(finaldata)
cellcolors <- ifelse(finaldata == 0, "white", 
                     ifelse(finaldata != 0 & rownames(finaldata) %in% c("Su60", "Su61", "Su62", "Su63", "Su64"),
                     "black", 
                     ifelse(finaldata != 0 & rownames(finaldata) %in% c("Su31", "Su37", "Su65"), 
                             "blue", "firebrick")
                            )
                     )
# Make the figure
tiff(filename = "matrix_revised.tiff", width = 12, height = 9, units = "in", res = 900, compression = "lzw")
color2D.matplot(finaldata, cellcolors = cellcolors, xlab = "", ylab = "Individual",
                axes = FALSE, cex.lab=0.8)
axis(side=1,at=0.5:ncol(finaldata),labels=colnames(finaldata), las=2, cex.axis=0.7)
axis(side=2, at=0.5:nrow(finaldata),labels = sort(rownames(finaldata), decreasing = TRUE), las=2, cex.axis=0.7)
dev.off()
