################################
#  LOAD FUNCTIONS & DATA
################################

source("dipcensR/Rcode/dipcensRfunctions.R")
setwd("dipcensR/Figures/insilicodilution")
dat <- read.csv("definetherain-master/data/Albumin/Alb 10e5/Big test 6.12.12_F08_Amplitude.csv")[,1]
set.seed(1)
thr <- kmeans(dat, quantile(dat, c(1/3, 2/3)))

# definetherain threshold: mean +- 3*SD of the negative/positive cluster
thr <- c(thr$centers[1,1] + 3*sd(dat[thr$cluster == 1]),
         thr$centers[2,1] - 3*sd(dat[thr$cluster == 2]))
thr <- mean(thr)

# perform simulation
posp <- dat[dat>=thr]
negp <- dat[dat<thr]
occs <- c(0.001, 0.01, 0.1, 0.2, 0.3, 0.6, 0.9)
res <- seq(1, 2, 1/3)
nsims <- 100
npart <- 15000
newDat <- c()
resList <- list()

plotlist <- list()
for(resc in 1:length(res)){
  posp2 <- (posp - mean(posp))*res[resc]+mean(posp)
  negp2 <- (negp - mean(negp))*res[resc]+mean(negp)
  resmat <- matrix(0, nrow = nsims, ncol = length(occs))
  for(occ in 1:length(occs)){
    cat(occ, " - ")
    for(sim in 1:nsims){
      newDat <- sample(c(sample(posp2, npart*occs[occ], replace = TRUE),
                         sample(negp2, npart*(1-occs[occ]), replace = TRUE)))
      checkThr <- checkThrSens(as.data.frame(newDat), thr,stepsize=0.01,
                               plot=TRUE, adjust = FALSE)
      resmat[sim,occ] <- checkThr$d
      thrBounds <- getThrBounds(as.data.frame(newDat), thr, 0.2)
      if((sim == 1) && (occ<5)){
        dat.df <- data.frame(Intensity=sample(newDat),x=1:length(newDat))
        plotlist[[occ]] <- checkThr$plot
        plotlist[[occ+4]] <- ggplot(dat.df, aes(x=x, y=Intensity))+
          geom_point(size=0.1)+
          theme(axis.text=element_text(size=8))+
          theme_minimal()+
          xlab("Partition ID")+
          ylab("Intensity")+
          geom_hline(yintercept=thr, colour=checkThr$qc, linewidth=1)+
          geom_hline(yintercept=thrBounds[[1]],
                     colour=checkThr$qc, linewidth=0.75,
                     linetype="dotted")+
          geom_hline(yintercept=thrBounds[[2]], colour=checkThr$qc,
                     linewidth=0.75, linetype="dotted")+
          ylim(5500,25500)

      }
    }
  }
  resList[[resc]] <- resmat
}


setwd("dipcensR/Figures/")
pdf("Figure4.pdf", width=12, height=8)
plot_grid(plotlist = list(plotlist[[1]][[1]], plotlist[[2]][[1]], plotlist[[3]][[1]], plotlist[[4]][[1]], plotlist[[5]], plotlist[[6]], plotlist[[7]], plotlist[[8]]),
          ncol = length(res), nrow=2, byrow = TRUE,
          labels = letters[1:8],
          rel_heights = c(1, 0.7))
dev.off()


# FIGURE 5

plotlist2 <- list()
for(resc in 1:length(res)){
  d.df <- data.frame(d = c(resList[[resc]]), Occupancy=as.factor(rep(occs, each = nsims)))

  plotlist2[[resc]] <- ggplot(d.df, aes(x=Occupancy, y=d))+
    geom_boxplot()+
    theme(axis.text=element_text(size=8))+
    theme_minimal()+
    geom_hline(yintercept = 0.1, col="orange")+
    geom_hline(yintercept = 0.2, col="red3")+
    scale_y_continuous(minor_breaks = seq(0.05, 10, 0.05),
                       breaks = seq(0, 10, 0.1),
                       limits=c(0, max(0.2, max(d.df$d))))+
    coord_flip()+
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.6), unit="cm"))


}

setwd("dipcensR/Figures/")
pdf("Figure5.pdf", width=12, height=5)
plot_grid(plotlist = plotlist2,
          ncol = 1, nrow=4, byrow = TRUE,
          labels = letters[1:4])
dev.off()
