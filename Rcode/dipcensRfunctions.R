sampleFluo <- function(cluster){
  sampleprob <- rep(1, length(cluster))
  fluo <- unlist(unname(cluster[sample(1:length(cluster), 1, prob = sampleprob)]))
  return(fluo)
}

remExtremes <- function(dat, thr){
  dat <- sort(dat)
  numAbove <- mean(dat>thr)*length(dat)
  numBelow <- mean(dat<thr)*length(dat)

  if(numAbove >= 10){
    # remove 10% maximum vals
    dat[c((length(dat)-floor(numAbove*0.1)):length(dat))] <- NA
  }

  if(numBelow >= 10){
    # remove 10% maximum vals
    dat[c(1:floor(numBelow*0.1))] <- NA
  }
  return(dat)
}
getRange <- function(dat, thr){
  # trim extremes (limited outliers can interfere with flagging)
  for(col in 1:ncol(dat)){
    dat[,col] <- remExtremes(dat[,col],thr[col])
  }

  range <- sapply(dat,max, na.rm=TRUE)-sapply(dat,min,na.rm=TRUE)
  return(range)
}
perturbOccup <- function(perturb, dat, thr){
  # trim extremes (limited outliers can interfere with flagging)
  for(col in 1:ncol(dat)){
    dat[,col] <- remExtremes(dat[,col],thr[col])
  }

  range <- sapply(dat,max, na.rm=TRUE)-sapply(dat,min,na.rm=TRUE)

  thr <- thr+unname(perturb*range)
  pertOccup <- getOccup(dat, thr,subset=FALSE)

  #range may be too large if threshold very low or high, e.g. 3 clusters
  # set to high/very low value to exclude region
  for(col in 1:ncol(dat)){
    if(thr[col]<=min(dat[,col],na.rm=TRUE)){
      pertOccup[col] <- 1e10
    }
    if(thr[col]>=max(dat[,col],na.rm=TRUE)){
      pertOccup[col] <- 1e-10
    }
  }
  return(pertOccup)
}
dipcensR <- function(dat, thr, plot = TRUE, window = 0.2, warn.level = 0.10, stop.level = 0.20,
                     adjust = FALSE, adjustLimit = 0.5, stepsize=0.01, flagNeg = FALSE){
  # check that data is data frame format
  if(is(dat)[1] != "data.frame"){
    errorCondition("Error: data must be a data frame\n")
  }
  # get perturbation values
  perbLevel <- seq(-0.5,0.5,stepsize)
  # calculate sensitivity curve
  sensCurv <- sapply(perbLevel, perturbOccup, dat, thr)
  if(is(sensCurv)[1] != "matrix"){
    # 1-color data, transform
    sensCurv <- t(data.frame(sensCurv))
  }
  # get the expected occupancy
  expec <- sensCurv[,(ncol(sensCurv)+1)/2]

  # initialize variables
  flagCol <- vector("character" ,ncol(dat))
  findInt <- window/stepsize
  rangeCol <- vector("numeric", ncol(dat))
  adjustThr <- vector("numeric", ncol(dat))
  rangeShift <- vector("numeric", ncol(dat))
  newThr <- thr

  # for each color, calculate sensitivity metric and flag color
  for(colid in 1:ncol(dat)){

    # get the first deriv., i.e., change in sensitivity within given window
    diffVals<-diff(sensCurv[colid,], findInt)

    ###########
    # if the threshold needs to be adjusted
    ##########
    if(adjust){
      # disregard values outside the adjust limits
      #   set them to infinity
      diffVals[union(which((sensCurv[colid,]/expec[colid])>(1+adjustLimit)),
                     which((sensCurv[colid,]/expec[colid])<1/(1+adjustLimit)))] <- Inf
      # what is the minimizing window location?
      getMin <- which.min(abs(diffVals))
      minLoc <- floor(median(which(abs(diffVals)==abs(diffVals)[getMin])))

      # get the d metric
      rangeCheck <- abs(diffVals[minLoc])
      rangeShift[colid] <- (sensCurv[colid,]/expec[colid])[minLoc+findInt/2]-1
      adjustThr[colid] <- perbLevel[minLoc+findInt/2]

      #########
      # then re-evaluate range (i.e., run the dipcensR algo again)
      #########
      # adjust the threhsold to the minimizing d metric one
      newThr[colid] <- thr[colid]+adjustThr[colid]*getRange(dat,thr)[colid]
      # recalculate sensitivity
      sensCurv2 <- sapply(perbLevel, perturbOccup, dat, newThr)
      if(is(sensCurv2)[1] != "matrix"){
        # if 1-color data, transform to data frame
        sensCurv2 <- t(data.frame(sensCurv2))
      }
      # get the expected occupancy based on the new threshold
      expec2 <- sensCurv2[,(ncol(sensCurv2)+1)/2]

      # get the minimizing range (no longer needed?)
      diffVals<-diff(sensCurv2[colid,], findInt)
      diffVals[union(which((sensCurv2[colid,]/expec2[colid])>(1+adjustLimit)),
                     which((sensCurv2[colid,]/expec2[colid])<1/(1+adjustLimit)))] <- Inf

      # get the d metric for the new threshold location
      rangeCheck <- sensCurv2[colid,as.integer(length(perbLevel)/2)-findInt/2]/expec2[colid] -
        sensCurv2[colid,as.integer(length(perbLevel)/2)+findInt/2]/expec2[colid]
      rangeCol[colid] <- rangeCheck

    ############
    # if the threshold is only checked for sensitivity towards changes
    ############
    } else {
      # just calculate the d metric in the window around the supplied
      #   threshold
      rangeCheck <- sensCurv[colid,as.integer(length(perbLevel)/2)-findInt/2]/expec[colid] -
        sensCurv[colid,as.integer(length(perbLevel)/2)+findInt/2]/expec[colid]
      rangeCol[colid] <- rangeCheck
    }

    # assign a color based on the estimated d metric
    #   and the warn/stop limits
    if(rangeCheck<warn.level){
      flagCol[colid] <- "green2"
    } else if (rangeCheck<stop.level){
      flagCol[colid] <- "orange"
    } else {
      flagCol[colid] <- "red3"
    }
    # check whether there are negative partitions
    if(sum(dat[,colid]>thr[colid])==length(dat[,colid])){
      # no negative partitions, flag as needing review (saturation)
      flagCol[colid] <- "red3"
    }
    # check whether there are positive partitions, and if they
    #   need to be flagged
    if(sum(dat[,colid]<thr[colid])==length(dat[,colid]) && flagNeg){
      # no positive partitions, and flagNeg is set to TRUE
      flagCol[colid] <- "red3"
    }
  }


  # get the complete vector with adjusted thresholds for plotting
  if(adjust){
    newThr <- thr+adjustThr*getRange(dat,thr)
  }

  # make some plots
  plotlist <- list()
  if(plot==TRUE){
    require("ggplot2")

    for(colid in 1:ncol(dat)){
      df.plot <- data.frame(perb = perbLevel*100,
                            sens = (sensCurv[colid,]/expec[colid])*100-100)
      if(adjust){
        plotlist[[colid]] <- ggplot(df.plot, aes(x=perb, y=sens))+
          theme_minimal()+
          scale_y_continuous(minor_breaks = seq(-60, 60, 5))+
          xlab("Perturbation (%, trimmed range)")+
          ylab("Partition occupancy change (%)")+
          geom_vline(xintercept = c(-window/2*100+adjustThr[colid]*100,
                                    window/2*100+adjustThr[colid]*100),
                     linetype = "dotted", colour="darkgrey")+
          geom_vline(xintercept=adjustThr[colid]*100, colour="darkgrey")+
          geom_hline(yintercept=c(-warn.level/2*100+rangeShift[colid]*100,
                                  warn.level/2*100+rangeShift[colid]*100),
                     linetype = "dashed", colour="orange2")+
          geom_hline(yintercept=rangeShift[colid]*100,
                     linewidth = 0.3, colour="darkgrey")+
          geom_hline(yintercept=c(-stop.level/2*100+rangeShift[colid]*100,
                                  stop.level/2*100+rangeShift[colid]*100),
                     linetype = "dashed", colour="orangered1")+
          # draw intersection of perturbation window and sensitivity curve
          #. tolerance issues, check that difference is less than 1e-5
          annotate("point",
                   x=-window/2*100+adjustThr[colid]*100,
                   y=df.plot$sens[which(abs(df.plot$perb - (-window/2*100+adjustThr[colid]*100)) < 1e-5)],
                   colour=flagCol[colid],
                   size = 1.8)+
          annotate("point",
                   x=window/2*100+adjustThr[colid]*100,
                   y=df.plot$sens[which(abs(df.plot$perb - (window/2*100+adjustThr[colid]*100)) < 1e-5)],
                   colour=flagCol[colid],
                   size = 1.8)+
          annotate("rect",
                   xmin=-window/2*100+adjustThr[colid]*100,
                   xmax=window/2*100+adjustThr[colid]*100,
                   ymin=df.plot$sens[which(abs(df.plot$perb - (window/2*100+adjustThr[colid]*100)) < 1e-5)],
                   ymax=df.plot$sens[which(abs(df.plot$perb - (-window/2*100+adjustThr[colid]*100)) < 1e-5)],
                   fill=flagCol[colid],
                   alpha = 0.3)+coord_cartesian(ylim=c(-60, 60))+
          annotate("text",
                   x = 49.5,
                   y = 57.5,
                   label=paste0("italic(d) == ",
                                formatC(rangeCol[colid], digits=2)),
                   parse=TRUE,
                   size = 3, hjust = "right")+
          annotate("text",
                   x = 49.5,
                   y = warn.level*100/2+2.5+rangeShift[colid]*100,
                   label=formatC(warn.level),
                   size = 2.5, hjust = "right",
                   colour="orange2")+
          annotate("text",
                   x = 49.5,
                   y = stop.level*100/2+2.5+rangeShift[colid]*100,
                   label=formatC(stop.level),
                   size = 2.5, hjust = "right",
                   colour="orangered1")+
          geom_line()
      } else {
        plotlist[[colid]] <- ggplot(df.plot, aes(x=perb, y=sens))+
          theme_minimal()+
          scale_y_continuous(minor_breaks = seq(-60, 60, 5))+
          xlab("Perturbation (%, trimmed range)")+
          ylab("Partition occupancy change (%)")+
          # draw perturbation window, and center of window
          geom_vline(xintercept = c(-window/2*100, window/2*100),
                     linetype = "dotted", colour="darkgrey")+
          geom_vline(xintercept=0, colour="darkgrey")+
          # draw warn/error levels, and center
          geom_hline(yintercept=c(-warn.level/2*100,warn.level/2*100),
                     linetype = "dashed", colour="orange2")+
          geom_hline(yintercept=0,
                     linewidth = 0.3, colour="darkgrey")+
          geom_hline(yintercept=c(-stop.level/2*100,stop.level/2*100),
                     linetype = "dashed", colour="orangered1")+
          # draw intersection of perturbation window and sensitivity curve
          #. tolerance issues, check that difference is less than 1e-5
          annotate("point",
                   x=-window/2*100,
                   y=df.plot$sens[which(abs(df.plot$perb - (-window/2*100)) < 1e-5)],
                   colour=flagCol[colid],
                   size = 1.8)+
          annotate("point",
                   x=window/2*100,
                   y=df.plot$sens[which(abs(df.plot$perb - (window/2*100)) < 1e-5)],
                   colour=flagCol[colid],
                   size = 1.8)+
          annotate("rect",
                   xmin=-window/2*100,
                   xmax=window/2*100,
                   ymin=df.plot$sens[which(abs(df.plot$perb - (window/2*100)) < 1e-5)],
                   ymax=df.plot$sens[which(abs(df.plot$perb - (-window/2*100)) < 1e-5)],
                   fill=flagCol[colid],
                   alpha = 0.3)+coord_cartesian(ylim=c(-60, 60))+
          annotate("text",
                   x = 49.5,
                   y = 57.5,
                   label=paste0("italic(d) == ",
                                formatC(rangeCol[colid], digits=2)),
                   parse=TRUE,
                   size = 3, hjust = "right")+
          annotate("text",
                   x = 49.5,
                   y = warn.level*100/2+2.5,
                   label=formatC(warn.level),
                   size = 2.5, hjust = "right",
                   colour="orange2")+
          annotate("text",
                   x = 49.5,
                   y = stop.level*100/2+2.5,
                   label=formatC(stop.level),
                   size = 2.5, hjust = "right",
                   colour="orangered1")+
          geom_line()
      }


    }
    for(colid in 1:ncol(dat)){

      if(adjust){
        df.plot <- data.frame(perb = perbLevel*100,
                              sens = (sensCurv2[colid,]/expec2[colid])*100-100)
        plotlist[[colid+ncol(dat)]] <- ggplot(df.plot, aes(x=perb, y=sens))+
          theme_minimal()+
          scale_y_continuous(minor_breaks = seq(-60, 60, 5))+
          xlab("Perturbation (%, trimmed range)")+
          ylab("Partition occupancy change (%)")+
          geom_vline(xintercept = c(-window/2*100,
                                    window/2*100),
                     linetype = "dotted", colour="darkgrey")+
          geom_hline(yintercept=c(-warn.level/2*100,
                                  warn.level/2*100),
                     linetype = "dashed", colour="orange2")+
          geom_hline(yintercept=c(-stop.level/2*100,
                                  stop.level/2*100),
                     linetype = "dashed", colour="orangered1")+
          # draw intersection of perturbation window and sensitivity curve
          #. tolerance issues, check that difference is less than 1e-5
          annotate("point",
                   x=-window/2*100,
                   y=df.plot$sens[which(abs(df.plot$perb - (-window/2*100)) < 1e-5)],
                   colour=flagCol[colid],
                   size = 1.8)+
          annotate("point",
                   x=window/2*100,
                   y=df.plot$sens[which(abs(df.plot$perb - (window/2*100)) < 1e-5)],
                   colour=flagCol[colid],
                   size = 1.8)+
          annotate("rect",
                   xmin=-window/2*100,
                   xmax=window/2*100,
                   ymin=df.plot$sens[which(abs(df.plot$perb - (window/2*100)) < 1e-5)],
                   ymax=df.plot$sens[which(abs(df.plot$perb - (-window/2*100)) < 1e-5)],
                   fill=flagCol[colid],
                   alpha = 0.3)+
          annotate("text",
                   x = 49.5,
                   y = 57.5,
                   label=paste0("italic(d) == ",
                                formatC(rangeCol[colid], digits=2)),
                   parse=TRUE,
                   size = 3, hjust = "right")+
          annotate("text",
                   x = 49.5,
                   y = warn.level*100/2+2.5,
                   label=formatC(warn.level),
                   size = 2.5, hjust = "right",
                   colour="orange2")+
          annotate("text",
                   x = 49.5,
                   y = stop.level*100/2+2.5,
                   label=formatC(stop.level),
                   size = 2.5, hjust = "right",
                   colour="orangered1")+
          coord_cartesian(ylim=c(-60, 60))+
          geom_line()
      }
      }

  }

  # return
  #   the flag color(s),
  #   the plot(s) (if requested)
  #   the d metric(s)
  #   the new threshold(s) (if requested)
  if(adjust){
    return(list(qc = flagCol,
                plot = plotlist,
                d = rangeCol,
                newThr=newThr))
  } else {
    return(list(qc = flagCol,
                plot = plotlist,
                d = rangeCol))
  }
}

getOccup <- function(testAmplitudes, thr, subset = FALSE){
  occup <- vector("numeric", length(thr))
  for(color in 1:length(thr)){
    #cat(color)
    complAmp <- testAmplitudes[complete.cases(testAmplitudes[,color]),]
    if(subset){
      complAmpSub <- partSubset(complAmp, color, thr)[,color]
      occup[color] <- mean(complAmpSub>thr[color], na.rm = TRUE)

    } else {
      occup[color] <- mean(testAmplitudes[,color]>thr[color], na.rm = TRUE)
    }

  }
  return(occup)
}
