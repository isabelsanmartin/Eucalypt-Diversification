# Modification to change the Y axis so that it does not start in zero (important for speciation and extinction rates if they never drop to zero)
# Works well but it does not show the complete plot with confidence intervals.
# A better solution is to modify the tess.plot.output.mod.R, changing these two lines if the y axis lower bound is 0.4 instead of zero):
# # ylim <- c(0.4, max(quantilesSpeciation, quantilesExtinction))
# ylim <- c(0.4, max(quantilesThisOutput))


#Usage
# Example call: speciation and extinction rates starting at 0.4, where out is "my.tess.output"
#source(tess.plot.output.ymin.R)
#tess.plot.output.mod.ymin(
#  output = out,
#  fig.types = c("speciation rates", "extinction rates"),
#  ymin = 0.4
# )

# Or: let other plots still start at 0
#tess.plot.output.mod.ymin(
#  output = my.tess.output,
#  fig.types = c("net diversification rates", "relative extinction rates"),
#  ymin = 0
#)

tess.plot.output.mod.ymin <- function(output, fig.types = NULL, ylim = NULL, ymin = 0) {
  
  # Loop over figure types
  for (type in fig.types) {
    
    # Select correct data for the plot
    if (type %in% names(output)) {
      dataThisOutput <- output[[type]]
    } else {
      stop(paste("Type", type, "not found in output"))
    }
    
    quantilesThisOutput <- apply(dataThisOutput, 2, quantile, prob = c(0.025, 0.975))
    
    # Set ylim dynamically unless provided
    if (is.null(ylim)) {
      if (type %in% c("speciation rates", "extinction rates")) {
        quantilesSpeciation <- apply(output[["speciation rates"]], 
                                     2, quantile, prob = c(0.025, 0.975))
        quantilesExtinction <- apply(output[["extinction rates"]], 
                                     2, quantile, prob = c(0.025, 0.975))
        ylim <- c(ymin, max(quantilesSpeciation, quantilesExtinction))
      } else {
        ylim <- c(ymin, max(quantilesThisOutput))
      }
    }
    
    # Plot
    plot(NA, xlim = c(0, ncol(dataThisOutput)), ylim = ylim,
         xlab = "Time", ylab = type, main = type)
    
    # Example: add median line
    medians <- apply(dataThisOutput, 2, median)
    lines(medians, col = "blue")
  }
}
