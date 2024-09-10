# load libraries, functions and data
library(lawstat)
source("Rcode/dipcensRfunctions.R")
dat <- read.csv("Data/trimmingExample.csv")

# initialize vectors to store ranges before and after trimming
rngvec <- vector("numeric", ncol(dat))
rngtrimvec <- vector("numeric", ncol(dat))

# for each reaction, calculate the range before and after trimming
# trimming is done using dipcensR's remExtremes function
for(coln in 1:ncol(dat)){
  rngvec[coln] <- max(dat[,coln], na.rm = TRUE) - min(dat[,coln], na.rm = TRUE)
  trim <- remExtremes(dat[,coln], 25000)
  rngtrimvec[coln] <- max(trim, na.rm = TRUE) - min(trim, na.rm = TRUE)
}

# combine ranges into a data frame for plotting / levene's test
rngdf <- data.frame(Trimmed = c(rep("yes", ncol(dat)), rep("no", ncol(dat))),
                    Range = c(rngtrimvec, rngvec))
boxplot(Range ~ Trimmed, data = rngdf)

# BF Levene test
levene.test(rngdf$Range, rngdf$Trimmed)

# calculate relative reduction in mad
1 -mad(rngtrimvec)/mad(rngvec)

