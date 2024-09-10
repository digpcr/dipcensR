################################
#  LOAD FUNCTIONS & DATA
################################

source("dipcensR/Rcode/dipcensRfunctions.R")
dat <- read.csv("definetherain-master/data/Albumin/Alb 10e5/Big test 6.12.12_F08_Amplitude.csv")[,1]
thr <- kmeans(dat, quantile(dat, c(1/3, 2/3)))

# definetherain threshold: mean +- 3*SD of the negative/positive cluster
thr <- c(thr$centers[1,1] + 3*sd(dat[thr$cluster == 1]),
         thr$centers[2,1] - 3*sd(dat[thr$cluster == 2]))
thr <- mean(thr)

setwd("definetherain-master/data/Albumin/all/")
files <- list.files(pattern="*.csv")
dvec <- vector("numeric", length(files))
qcvec <- vector("character", length(files))
# keep: the positive control
#       orange/red flagged reactions
keep <- c(20, 18, 21, 27, 4, 8, 28)
keeplist <- list()
longdat <- c()
for(file in 1:length(files)){
  dat <- read.csv(files[file])[,1]
  longdat <- c(longdat,c(sample(dat)))
  checkThr <- checkThrSens(as.data.frame(dat), thr, plot = TRUE, stepsize=0.01, plotRowSingle=TRUE,
                           adjust=FALSE)
  dvec[file] <- checkThr$d
  qcvec[file] <- checkThr$qc
  df.dp <- data.frame(dp = sample(dat), x=1:length(dat))
  plots <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
    theme_minimal()+
    xlab("Partition ID")+
    ylab("Intensity")+
    geom_hline(yintercept=thr, colour=checkThr$qc, size=1)+ylim(3500, 25000)+
    theme(axis.text=element_text(size=8),
          plot.margin = margin(2, 0.5, 2, 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if(file %in% keep){
    if(file == 20){
      plots <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
        theme_minimal()+
        xlab("")+
        ylab("Intensity")+
        geom_hline(yintercept=thr, colour=checkThr$qc, size=1)+ylim(3500, 25000)+
        theme(axis.text=element_text(size=8),
              plot.margin = margin(2, 0.5, 2, 0.5),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))      
    } else {
      plots <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
        theme_minimal()+
        xlab("")+
        ylab("")+
        geom_hline(yintercept=thr, colour=checkThr$qc, size=1)+ylim(3500, 25000)+
        theme(axis.text=element_text(size=8),
              plot.margin = margin(2, 0.5, 2, 0.5),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))      
    }
  }
  if(file==27){
    plots <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
      theme_minimal()+
      xlab("Partition ID")+
      ylab("")+
      geom_hline(yintercept=thr, colour=checkThr$qc, size=1)+ylim(3500, 25000)+
      theme(plot.margin = unit(c(0.1,0.1,0,1), "lines"))+
      theme(axis.text=element_text(size=8),
            plot.margin = margin(2, 0.5, 2, 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  if(file==28){
    plots <- ggplot(df.dp, aes(x=x, y=dp))+geom_point(size=0.1)+
      theme_minimal()+
      xlab("")+
      ylab("")+
      geom_hline(yintercept=thr, colour=checkThr$qc, size=1)+ylim(3500, 25000)+
      theme(plot.margin = unit(c(0.1,0.1,0.1,1), "lines"))+
      theme(axis.text=element_text(size=8),
            plot.margin = margin(2, 0.5, 2, 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  if(file %in% keep){
    cat(file)
    keeplist[[which(file==keep)]] <- plots
  }
  pdf(paste0(file, ".pdf"), width=8, height=4)
  print(plots)
  checkThr$plot[[1]]
  dev.off()
  
}
setwd("dipcensR/Figures")
pdf("Figure6.pdf", height=4, width=12)
plot_grid(plotlist = list(keeplist[1][[1]], 
                          keeplist[2][[1]],
                          keeplist[3][[1]],
                          keeplist[4][[1]],
                          keeplist[7][[1]],
                          keeplist[5][[1]],
                          keeplist[6][[1]]),
          ncol= 7, labels=letters[1:7])
dev.off()

conc <- c(rep(c(10e4, 10e3, 10e1, 0), 4),
          rep(c(10e5, 10e2, 10e0), 4))
df <- data.frame(qcvec,conc, dvec)
df[order(df$conc),]
