################################
#  LOAD FUNCTIONS & LIBRARIES
################################

source("dipcensR/Rcode/dipcensRfunctions.R")

dat <- read.csv("dipcensR/Data/example1.csv")
thr <- read.csv("dipcensR/Data/example1thr.csv")
thr <- unname(unlist(c(thr)))
thr[1] <- 25000
thr[2] <- 13000
thr[3] <- 2800
thr[4] <- 26000
checkThr <- checkThrSens(dat, thr, plot = TRUE, stepsize=0.01, plotRowSingle=TRUE,
                         adjust=FALSE)
checkThrAdjust <- checkThrSens(dat, thr, plot = TRUE, stepsize=0.01, plotRowSingle=TRUE,
                               adjust=TRUE)
plots <- list()
for(colid in 1:ncol(dat)){
  df.dp <- data.frame(dp = sample(dat[,colid]), x=1:nrow(dat))
  plots[[colid]] <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
    theme_minimal()+
    xlab("Partition ID")+
    ylab("Intensity")+
    geom_hline(yintercept=thr[colid], colour=checkThr$qc[colid], size=1)
}

for(colid in 1:ncol(dat)){
  df.dp <- data.frame(dp = sample(dat[,colid]), x=1:nrow(dat))
  plots[[colid+4]] <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
    theme_minimal()+
    xlab("Partition ID")+
    ylab("Intensity")+
    geom_hline(yintercept=checkThrAdjust$newThr[colid], colour=checkThrAdjust$qc[colid], size=1)
}
plot_grid(plotlist=plots, nrow = 2, ncol=4, labels = letters[1:8])
