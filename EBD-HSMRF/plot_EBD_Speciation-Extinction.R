################################################################################
#
# Plot correlation factors of episodic birth-death process
# with environmental variables
#
#
# authors: Sebastian Hoehna
#
################################################################################

library(RevGadgets)
library(ggplot2)


# specify the input file
file <- paste0("output/Eucalypts_ML_noOut_trueEuc_EBD.log")

trace_env <- readTrace(path = file, burnin = 0.25)

# plot the prior vs the posterior
plot <- plotTrace(trace_env, vars=c("speciation_global_scale", "extinction_global_scale"))[[1]]  +
     # modify legend location using ggplot2
     theme(legend.position.inside = c(0.80,0.80))

ggsave(paste0("Eucalypts_ML_noOut_trueEuc_EBD.png"), plot, height=5, width=5)

