# Install RevGadget if you haven't done so already
#library(devtools)
#install_github("revbayes/RevGadgets")

#setwd("/Users/isabelsanmartin/Documents/Eucalyptus_Pugnaire/ANALYSES/EBD_NonEnviron_Eucalypts_2016_dated_r8s_ultrametric_binary_noOut_trueEuc_LARGEVARIANCE")

library(RevGadgets)
library(ggplot2)

# specify the output files

# IMPORTANT!! In "speciation_times.log" and "extinction_times.log" files (including the run_1, run_2, etc. files), replace the word "extinction_" and "speciation_" in the column headers by "interval_". This is because the original function processDivRates in RevGadgets searches for this word, but the RevBayes new script uses the parameter name: speciation_times, extinction_times

### All #---------------------------------------------------------------------------------------------------

speciation_time_file <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_times.log"
speciation_rate_file <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_rates.log"
extinction_time_file <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_times.log"
extinction_rate_file <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_rates.log"

# read in and process rates
rates <- processDivRates(speciation_time_log = speciation_time_file,
                          speciation_rate_log = speciation_rate_file,
                          extinction_time_log = extinction_time_file,
                          extinction_rate_log = extinction_rate_file,
                          burnin = 0.25,
                          summary = "median")

# plot rates through time
p <- plotDivRates(rates = rates) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")
p
ggsave("Eucalypts_ML_noOut_trueEuc_EBD_LargeVariance.png", p)

### Modify the X axis geological scale to make it more detailed in ggplot2 #####
# 'scale_x_continuous() customizes the x-axis scale
# 'trans = "reverse" argument reverses the scale of the x-axis scale
# ´custom_breaks() function maintains fine control over tick mark placement. Generates breaks
# Adjust number of breaks or ticks as needed. Here, "n = 20" introduces 20 evenly breaks over x length.

# Assuming p is your ggplot object

# Define a custom function for breaks
custom_breaks <- function(limits) {
  pretty(limits, n = 20)
}

# Reverse the x-axis scale
p <- p + 
  scale_x_continuous(
    breaks = custom_breaks,
    trans = "reverse"
  )

p
ggsave("Eucalypts_ML_noOut_trueEuc_EBD_LargeVariance.png", p)


### Run 1 #---------------------------------------------------------------------------------------------------

speciation_time_file1 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_times_run_1.log"
speciation_rate_file1 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_rates_run_1.log"
extinction_time_file1 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_times_run_1.log"
extinction_rate_file1 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_rates_run_1.log"

# read in and process rates
rates1 <- processDivRates(speciation_time_log = speciation_time_file1,
                         speciation_rate_log = speciation_rate_file1,
                         extinction_time_log = extinction_time_file1,
                         extinction_rate_log = extinction_rate_file1,
                         burnin = 0.25,
                         summary = "median")

# plot rates through time
p1 <- plotDivRates(rates = rates1) +
        xlab("Millions of years ago") +
        ylab("Rate per million years")
p1
ggsave("Eucalypts_ML_noOut_trueEuc_EBD_run_1_LargeVariance.png", p1)

### Run 2 #---------------------------------------------------------------------------------------------------

speciation_time_file2 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_times_run_2.log"
speciation_rate_file2 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_rates_run_2.log"
extinction_time_file2 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_times_run_2.log"
extinction_rate_file2 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_rates_run_2.log"

# read in and process rates
rates2 <- processDivRates(speciation_time_log = speciation_time_file2,
                          speciation_rate_log = speciation_rate_file2,
                          extinction_time_log = extinction_time_file2,
                          extinction_rate_log = extinction_rate_file2,
                          burnin = 0.25,
                          summary = "median")

# plot rates through time
p2 <- plotDivRates(rates = rates2) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")
p2
ggsave("Eucalypts_ML_noOut_trueEuc_EBD_run_2_LargeVariance.png", p2)

### Run 3 #---------------------------------------------------------------------------------------------------

speciation_time_file3 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_times_run_3.log"
speciation_rate_file3 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_rates_run_3.log"
extinction_time_file3 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_times_run_3.log"
extinction_rate_file3 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_rates_run_3.log"

# read in and process rates
rates3 <- processDivRates(speciation_time_log = speciation_time_file3,
                          speciation_rate_log = speciation_rate_file3,
                          extinction_time_log = extinction_time_file3,
                          extinction_rate_log = extinction_rate_file3,
                          burnin = 0.25,
                          summary = "median")

# plot rates through time
p3 <- plotDivRates(rates = rates3) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")
p3
ggsave("Eucalypts_ML_noOut_trueEuc_EBD_run_3_LargeVariance.png", p3)

### Run 4 #---------------------------------------------------------------------------------------------------

speciation_time_file4 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_times_run_4.log"
speciation_rate_file4 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_speciation_rates_run_4.log"
extinction_time_file4 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_times_run_4.log"
extinction_rate_file4 <- "output/Eucalypts_ML_noOut_trueEuc_EBD_extinction_rates_run_4.log"

# read in and process rates
rates4 <- processDivRates(speciation_time_log = speciation_time_file4,
                          speciation_rate_log = speciation_rate_file4,
                          extinction_time_log = extinction_time_file4,
                          extinction_rate_log = extinction_rate_file4,
                          burnin = 0.25,
                          summary = "median")

# plot rates through time
p4 <- plotDivRates(rates = rates4) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")
p4
ggsave("Eucalypts_ML_noOut_trueEuc_EBD_run_4_LargeVariance.png", p4)
