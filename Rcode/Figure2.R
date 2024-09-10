################################
#  LOAD FUNCTIONS & LIBRARIES
################################

source("dipcensR/Rcode/dipcensRfunctions.R")

dat <- read.csv("dipcensR/Data/example1.csv")
thr <- read.csv("dipcensR/Data/example1thr.csv")
thr <- unname(unlist(c(thr)))
checkThr <- checkThrSens(dat, thr, plot = TRUE, stepsize=0.01, plotRowSingle=TRUE,
                         adjust=FALSE)

thrBounds <- getThrBounds(dat, thr, 0.2)
plots <- checkThr$plot
lpl <- length(plots)
for(colid in 1:ncol(dat)){
  df.dp <- data.frame(dp = sample(dat[,colid]), x=1:nrow(dat))
  plots[[colid+lpl]] <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
    theme(axis.text=element_text(size=8))+
    theme_minimal()+
    xlab("Partition ID")+
    ylab("Intensity")+
    geom_hline(yintercept=thr[colid], colour=checkThr$qc[colid], linewidth=1)+
    geom_hline(yintercept=thrBounds[[1]][colid],
               colour=checkThr$qc[colid], linewidth=0.75,
               linetype="dotted")+
    geom_hline(yintercept=thrBounds[[2]][colid], colour=checkThr$qc[colid],
               linewidth=0.75, linetype="dotted")
}
setwd("dipcensR/Figures/")
pdf("Figure2.pdf", width=12, height=7)
plot_grid(plotlist=plots, nrow = 2, ncol=4, labels = letters[1:8])
dev.off()
