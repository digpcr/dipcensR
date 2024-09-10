source("/dipcensR/Rcode/dipcensRfunctions.R")
dat <- read.csv("dipcensR/Data/example1.csv")
thr <- read.csv("dipcensR/Data/example1thr.csv")
thr <- thr$x

# check sensitivity to parameter changes
windowrange <- seq(0.01, 0.9, 0.01)
flagfactor <- c(seq(0.1, 2, 0.1))
flagres <- matrix("", nrow = length(windowrange), ncol = length(flagfactor))
for(flagfac in 1:length(flagfactor)){
  cat(flagfactor[flagfac], " - ")
  for(win in 1:length(windowrange)){
    flagres[win, flagfac] <- checkThrSens(data.frame(dat[,1]), thr[1], window = windowrange[win], warn.level = 0.1*flagfactor[flagfac], stop.level = 0.2*flagfactor[flagfac])$qc
  }  
}



# check sensitivity to parameter changes
windowrange <- seq(0.01, 0.9, 0.01)
flagfactor <- c(seq(0.1, 2, 0.1))
flagres2 <- matrix("", nrow = length(windowrange), ncol = length(flagfactor))
for(flagfac in 1:length(flagfactor)){
  cat(flagfactor[flagfac], " - ")
  for(win in 1:length(windowrange)){
    flagres2[win, flagfac] <- checkThrSens(data.frame(dat[,2]), thr[2], window = windowrange[win], warn.level = 0.1*flagfactor[flagfac], stop.level = 0.2*flagfactor[flagfac])$qc
  }  
}


# check sensitivity to parameter changes
windowrange <- seq(0.01, 0.9, 0.01)
flagfactor <- c(seq(0.1, 2, 0.1))
flagres3 <- matrix("", nrow = length(windowrange), ncol = length(flagfactor))
for(flagfac in 1:length(flagfactor)){
  cat(flagfactor[flagfac], " - ")
  for(win in 1:length(windowrange)){
    flagres3[win, flagfac] <- checkThrSens(data.frame(dat[,3]), thr[3], window = windowrange[win], warn.level = 0.1*flagfactor[flagfac], stop.level = 0.2*flagfactor[flagfac])$qc
  }  
}


# check sensitivity to parameter changes
windowrange <- seq(0.01, 0.9, 0.01)
flagfactor <- c(seq(0.1, 2, 0.1))
flagres4 <- matrix("", nrow = length(windowrange), ncol = length(flagfactor))
for(flagfac in 1:length(flagfactor)){
  cat(flagfactor[flagfac], " - ")
  for(win in 1:length(windowrange)){
    flagres4[win, flagfac] <- checkThrSens(data.frame(dat[,4]), thr[4], window = windowrange[win], warn.level = 0.1*flagfactor[flagfac], stop.level = 0.2*flagfactor[flagfac])$qc
  }  
}

par(mfrow = c(1, 4))
plot(0, 0, type = "n", bty = "n",
     xlim = c(0, max(flagfactor)),
     ylim = c(min(windowrange), max(windowrange)),
     xlab = "Fold-change in warning/stop level",
     ylab = "Window size")
for(flagfac in 1:length(flagfactor)){
  for(win in 1:length(windowrange)){
    points(flagfactor[flagfac], windowrange[win], pch = 15, col = flagres[win, flagfac], cex = 3.5)
  }  
}
abline(v = 1, lty = 3, col = "black")
abline(h = 0.2, lty = 3, col = "black")
abline(a=0,b=0.4/2)


plot(0, 0, type = "n", bty = "n",
     xlim = c(0, max(flagfactor)),
     ylim = c(min(windowrange), max(windowrange)),
     xlab = "Fold-change in warning/stop level",
     ylab = "Window size")
for(flagfac in 1:length(flagfactor)){
  for(win in 1:length(windowrange)){
    points(flagfactor[flagfac], windowrange[win], pch = 15, col = flagres2[win, flagfac], cex = 3.5)
  }  
}
abline(v = 1, lty = 3, col = "black")
abline(h = 0.2, lty = 3, col = "black")
abline(a=0,b=0.4/2)

plot(0, 0, type = "n", bty = "n",
     xlim = c(0, max(flagfactor)),
     ylim = c(min(windowrange), max(windowrange)),
     xlab = "Fold-change in warning/stop level",
     ylab = "Window size")
for(flagfac in 1:length(flagfactor)){
  for(win in 1:length(windowrange)){
    points(flagfactor[flagfac], windowrange[win], pch = 15, col = flagres3[win, flagfac], cex = 3.5)
  }  
}
abline(v = 1, lty = 3, col = "black")
abline(h = 0.2, lty = 3, col = "black")
abline(a=0,b=0.4/2)

plot(0, 0, type = "n", bty = "n",
     xlim = c(0, max(flagfactor)),
     ylim = c(min(windowrange), max(windowrange)),
     xlab = "Fold-change in warning/stop level",
     ylab = "Window size")
for(flagfac in 1:length(flagfactor)){
  for(win in 1:length(windowrange)){
    points(flagfactor[flagfac], windowrange[win], pch = 15, col = flagres4[win, flagfac], cex = 3.5)
  }  
}
abline(v = 1, lty = 3, col = "black")
abline(h = 0.2, lty = 3, col = "black")
abline(a=0,b=0.4/2)
