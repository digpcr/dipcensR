# dipcensR
Evaluation and optimization of partition classification robustness

In this repository you can find
- the dipcensR R functions (/Rcode/dipcensRfunctions.R)
- example data used in the paper (Data folder)
- R code used to generate results in the paper (/Rcode/FigureX.R)

Note that experimental dilution series data from Jones et al. (doi 10.1016/j.jviromet.2014.02.020) is used, available from https://github.com/jacobhurst/definetherain

# Example

Download the dipcensRfunctions.R code and example data to your machine.

A simple example analysis is as follosws:

	# load libraries
	# the ggplot library is only required when dipcensR curve plots are generated
	library(ggplot2)

    # load the dipcensR functions
    setwd("/change/to/your/working/directory")
	source("Rcode/dipcensRfunctions.R")

	# load the example data

	dat <- read.csv("Data/example1.csv")
	thr <- read.csv("Data/example1thr.csv")
	thr <- unname(unlist(c(thr)))

	# check robustness 
	checkThr <- dipcensR(dat, thr)

dipcensR requires:
- partition intensities (numeric), in a data frame format (one column per threshold to evaluate)
- thresholds (numeric), in vector format

Some input parameters can optionally be altered, it is recommended to start an initial evaluation using default parameters.
- plot (default TRUE): whether to generate dipcensR curve plots
- window (default 0.2): perturbation window to consider for assigning a flag
- warn.level (default 0.1): relative occupancy change to exceed to assign a warning flag (orange) (and as long as the stop.level is not exceeded; 0.10 meaning 10% change within the window (argument))
- stop.level (default 0.2): relative occupancy change to exceed to assign a stop flag (red; 0.20 meaning 20% change within the window (argument))
- adjust (default FALSE): whether to adjust the thresholds to their location of maximum robustness, within adjustment limits (adjustLimit argument)
- adjustLimit (default 0.5): relative partition occupancy change allowed, in which to find the robustness maximum
- stepsize (default 0.01): step size of the threshold perturbation (0.01 meaning steps of 1% perturbation)
- flagNeg (default FALSE): whether to flag partitions as stop (red) when there are no positive partitions, even if threshold perturbation results in a green/orange flag

The `checkThr` list contains the following items:
- "qc": the resulting flags (green/orange/red), one for each column of `dat`
- "plots": if requested, the dipcensR curve plots, one for each column of `dat`
- "d": the dipcensR metrics, one for each column of `dat`
- "newThr": if requested, robustness-maximized thresholds, one for each column of `dat`

Make some plots:

	# visualize sensitivity of example dataset
	# use the cowplot library for generating a multi-row, labeled plot
	library(cowplot)

	# get the dipcensR plots
	plots <- checkThr$plot
	# get the number of plots
	lpl <- length(plots)

	# generate partition ID vs. intensity plots, and plot the threshold line
	#   use the flag color for the threshold line
	for(colid in 1:ncol(dat)){
	  # get one column of data
	  df.dp <- data.frame(dp = sample(dat[,colid]), x=1:nrow(dat))

	  # generate a plot of intensities
	  plots[[colid+lpl]] <- ggplot(df.dp, aes(x=x, y=dp))+
	    geom_point(size=0.1)+ # draw intensities as points
	    theme_minimal()+      # use a minimal theme
	    xlab("Partition ID")+ # x axis is for the partitions
	    ylab("Intensity")+    # y axis is for the intensities
	    geom_hline(yintercept=thr[colid], colour=checkThr$qc[colid], size=1)
	                          # add a threshold line, color is the flag color
	}

	# use the cow plot library to generate a 2 by 4 plot
	plot_grid(plotlist=plots, nrow = 2, ncol=4, labels = letters[1:8])

![Example dipcensR plots](/Figures/example.png "Example dipcensR plots")

These can be interpreted as explained in the Figure below. Calculations are explained in detail in the dipcensR manuscript.

![dipcensR plot interpretation](/Figures/interpretation.png "dipcensR interpretation")
