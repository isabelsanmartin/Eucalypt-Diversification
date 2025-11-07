###
# Reversible Jump MCMC algorithm:
# applied to estimate the diversification rates over time.
#
# author: Sebastian Hoehna
#
#library(TESS) loads packages ape and coda and desolve as dependencies
setwd("~/Documents/Eucalyptus_Pugnaire/ANALYSES/COMET")
getwd()

library(TESS)


#### OBS!!! ORIGINAL TREES ARE NOT ULTRAMETRIC AND NOT BINARY ########################
##### # If is.ultrametric(tree) AND is.binary tree (TRUE), move to step (B)

#### (A) TREE IS NOT ULTRAMETRIC AND NON BINARY ######################################

#Error message in tess.analysis(tree = tree, numExpectedRateChanges = EXPECTED_NUM_EVENTS,  : 
# The likelihood function is only defined for ultrametric trees! ###

tree <- read.tree("data/Eucalypts_2016_dated_r8s.tre")

is.binary(tree)
FALSE
is.ultrametric(tree)
FALSE


### (A) COERCE TREE TO BE BINARY  ########################################
library(phytools)
library(ape)

# Uses ape function multi2di, which solves polytomies by assigning nodes with more thant two descendants, branches of length zero

tree <- multi2di(tree, tol = 1e-08)
is.binary(tree)
TRUE


### COERCE TREE TO BE ULTRAMETRIC ########################################

# Uses phytools force.ultrametric function                            *
# NOTE *    force.ultrametric does not include a formal method to    *
#      *    ultrametricize a tree & should only be used to coerce    *
#      *   a phylogeny that fails is.ultrametric due to rounding --  *
#      *    not as a substitute for formal rate-smoothing methods. 



# NOT WORKS <-  default method

# Implements the method = "nnls" from the phangorn package to compute the set of edge lengths that result in a minimized sum-of-squares distance 
# between the patristic distance of the output and input trees. Works but makes strange artifacts in the tree, creating a basal polytomy
# tree2 <- force.ultrametric(tree, method = "nnls") 

## Implements the method = "extend" by simply extending all the external edges of the tree to match the external edge with the greatest total height ##
#  Works! The tree has the same topology as original #
tree <- force.ultrametric(tree, method = "extend")

is.ultrametric(tree)
TRUE

is.binary(tree)
TRUE

write.tree(tree, file = "data/Eucalypts_2016_dated_r8s_ultrametric_binary_noOut_trueEuc.tre")



### SECTION B ## 

#This command assigns the character string "Eucalypts_Mr_Bayes_dated_r8s" to the variable treefile. It creates a variable named treefile and stores the specified #string in it. This variable could be used later in the code to refer to the tree file.

treefile <- "Eucalypts_2016_dated_r8s_ultrametric_binary_noOut_trueEuc"

# This command uses the sprintf function to create a character string with a formatted placeholder. The %s placeholder is replaced by the value of the treefile variable. Specifically, it creates a string starting with the sprintf placeholder,"CoMET_noMEE", or "CoMET_MEE", or "CoMET_Div" (depending on the analysis), followed by the value stored in treefile variable above.

analysis_name_noMEE <- sprintf("CoMET_noMEE_%s",treefile)
analysis_name_MEE <- sprintf("CoMET_MEE_%s",treefile)
analysis_name_Div <- sprintf("CoMET_Div_%s",treefile)

CDT <- "survival"
EXPECTED_NUM_EVENTS <- 2
MCMC_ITERATIONS <- 1000000
ALLOW_MASS_EXTINCTION <- !TRUE

# This section below is not needed if we specify the "analysis_name" as above (no-MEE, MEE, Div)
#if ( ALLOW_MASS_EXTINCTION == TRUE ) {
#   analysis_name <- sprintf("CoMET_%s_ME",treefile)
#} else {
#   analysis_name <- sprintf("CoMET_noMEE_%s",treefile)
#}

# Now, we read the tree stored in the treefile name using the sprintf function. 

tree <- read.tree(file=sprintf("data/%s.tre",treefile) )

# Alternatively, we can just read the tree (easier!!!)
tree <- read.tree("data/Eucalypts_2016_dated_r8s_ultrametric_binary_noOut_trueEuc.tre")



################################################# Assign sampling value at present #############################################

#If there is complete taxon sampling, rho <- 1.0, we do not need to provide this argument

# If there is incomplete taxon sampling, and taxon sampling has been "uniform" across clades (default strategy), we need to estimate the paramete rho, which is the ratio between the number of species sampled in the phylogeny, and the Total diversity in the case-study lineage (the number of extant species classified as belonging to the lineage in the current taxonomy). In the case of eucalytus, there are c. 800 species in the true eucalyptus (Thornhill et al. 2019). The number of tips in the "Eucalypts_2016_dated_r8s_ultrametric_binary_noOut_trueEuc.tre" phylogeny is 711. Therefore, rho = 0.88875.

total <- 800
rho <- (tree$Nnode+1)/total

# In our trees, the sampling is as following. Notice how non-binary trees have equal number of tips but fewer nodes than binary trees (polytomies)
# rho
# [1] 0.71125 (MrBayes tree ultrametric but not binary: 674 tips and 568 internal nodes. NON BINARY. ULTRAMETRIC)
# [1] 0.8425 (MrBayes tree ultrametric binary: 674 tips and 673 internal nodes. BINARY AND ULTRAMETRIC)
# [1] 0.91375 (ML 1 tree with outgroups - Sebastian bifurcating: 731 tips and 730 internal nodes. BINARY AND ULTRAMETRIC)
# [1] 0.895 (ML 1 tree without outgroups but including mesoeucalypts: 716 tips and 715 internal nodes. BINARY AND ULTRAMETRIC)
# [1] 0.88875 (ML 1 tree without outgroups but without mesoeucalypts: 711 tips and 710 internal nodes. BINARY AND ULTRAMETRIC) <- Tree used in EBD analyses


# Assign priorForms; a distribution prior is used for speciation and extinction rates, defined by a mean=X.Z and stdev=X.X, which will be estimated from the data in a preliminary (burnin) phase. These fixed values wil be used as hyper priors for the lognormal distribution used in the final run. 

# Number of rate shifts in extinction and speciation and number MEEs are modeled by three independent Compound Poisson Process with a distribution hyperprior with lambda=2 (0.5 probability given to 0 rate shifts)

#If you do not want to fix the prior distribution but let TESS estimate it with empiricalHyperPriors=TRUE, then leave open the brackets as in TESS manual.
# priorForms <- c()

# Or tell TESS to explore all three alternative distributions

priorForms <- c("lognormal","normal","gamma")

# In this case, we use lognormal priors as in the original tutorial because this is the most common distribution used for scale parameters (conservative). COMET will estimate these empirical hyperpriors in a burnin phase.

# Prior forms does not seem to have a large impact on analysis. We set up for all rate parameters to adopt lognormal prior distributions, the most common.

priorForms <- c("lognormal")

#################################################### RUN ANALYSIS ##########################################################
#To run an analysis with no MEEs, use the code [1] below. This will produce an output called "CoMET_noMEE" containing the files with parameter values. It will also produce at the end a PDF figure named "CoMET_noMEE.pdf" containing vignettes for 4 figures. This function will create a directory with the name of "analysis_name" given above in the current working directory getwd()

##############################
##### [1] WITH NO MEEs #######
##############################

tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name_noMEE)  

# To modify the empirical Hyperprior values (e.g., to increase the variance), you can change the arguments to speficy the user-defined hyperprior distributions: speciationRatePriorMean, speciationRatePriorStDev, extinctionRatePriorMean, extinctionRatePriorStDev, providing appropriate values, and also change the argument empiricalHyperPriorForm to FALSE to disallow the empirical Bayesian estimation.
# Example to correct very low ESS values possible due to large variance. DOES NOT WORK VERY WELL. WORSENS DIAGNOSTICS - ESS VALUES.

tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriors = FALSE, speciationRatePriorMean = 0.07895113, speciationRatePriorStDev = 0.5, extinctionRatePriorMean = 0.04307893, extinctionRatePriorStDev = 0.5, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name_noMEE)  

##############################
##### [2] WITH MEEs #######
##############################
#To run an analysis with MEEs, use the code [2] below. This will produce an output folder called "CoMET_ME_pruned" containing the files with parameter posterior MCMC values. It will also produce at the end a PDF figure named "CoMET_ME_pruned.pdf" containing vignettes for six figures.


ALLOW_MASS_EXTINCTION <- TRUE

tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name_MEE)
          
################################# CHANGE SAMPLING STRATEGY FROM UNIFORM TO DIVERSIFIED ####################################################################
# INCOMPLETE TAXON SAMPLING strategy. Default in tess.analysis function is "samplingStrategy=uniform" (species are sampled uniformly at random), which means that each species has the same probability (rho) of being sampled in the phylogeny. 
# To apply alternative: "samplingStrategy=diversified" (species are sampled to maximize the diversity sampled in the phylogeny: i.e., only the oldest 25% of divergence events are included in the reconstructed phylogeny with sampling probability rho, and all later divergence events are excluded), we can use the function "tess.analysis.diversified.R" created by modifying the argument "samplingStrategy=diversified" in the "tess.likelihood" function part of the code. First, source it:

source("tess.analysis.diversified.R")

tess.analysis.diversified(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, initialSpeciationRate=2.0, initialExtinctionRate=1.0, empiricalHyperPriorInflation = 10.0, empiricalHyperPriorForm = priorForms, samplingProbability=rho, estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, dir = analysis_name_Div)   

######################################################################################################################## 
###########                           PROCESS ANALYSIS OUTPUT                                                          
#########################################################################################################################


############################################# WITH NO MEEs ########################################################


out <- tess.process.output(analysis_name_noMEE, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

# If there are two analyses run: run1 and run2, then we can use the following commands to construct a combined sample to plot

output1 <- tess.process.output(analysis_name1, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)
output2 <- tess.process.output(analysis_name2, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

outputs <- list(output1, output2)

#And use the multichain diagnostics below for these analyses

test.plot.multichain.diagnostics(outputs)


############################################# WITH MEEs ########################################################

outMEE <- tess.process.output(analysis_name_MEE, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, numExpectedMassExtinctions=EXPECTED_NUM_EVENTS)

#Arguments: tree = NULL, numExpectedRateChanges = numExpectedMassExtinctions = EXPECTED_NUM_EVENTS = 2, burnin = 0.25, numIntervals = 100, criticalBayesFactors = c(2,6,10)


### With NO-MEE DIVERSIFIED STRATEGY

outDiv <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

###########################################################################################################
#####                              CHECK MCMC DIAGNOSTICS WITH CODA
###########################################################################################################

#If only one chain and with EES diagnostic
# Remember to give only one parameter or divide into vignettes with "par" function. Otherwise, plots will be overlaid. Similarly, choose between "EES" and "geweke" to avoid overlaying of values. There are other arguments, which we do not consider here: (col=NULL,xaxt="n",yaxt="s")


##################################################### Plot the result diagnostics ######################################

##### WITH NON-MEE 

NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","net-diversification rates","relative-extinction rates", "speciation shift times","extinction shift times")

pdf(sprintf("%s-diagnosticsCODA.pdf",analysis_name_noMEE))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)


tess.plot.singlechain.diagnostics(out,parameters=c("speciation rates", "extinction rates", "net-diversification rates", "relative-extinction rates", "speciation shift times", "extinction shift times"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()

##### WITH MEE

NUM_FIGS <- 8
FIG_TYPES <- c("speciation rates", "extinction rates", "speciation shift times", "extinction shift times", "net-diversification rates", "relative-extinction rates", "mass extinction times", "empty")

pdf(sprintf("%s-diagnosticsCODA.pdf", analysis_name_MEE))
layout.mat <- matrix(1:NUM_FIGS,nrow=4,ncol=NUM_FIGS / 4)
layout(layout.mat)

tess.plot.singlechain.diagnostics(outMEE, parameters=c("speciation rates", "extinction rates", "speciation shift times", "extinction shift times", "net-diversification rates", "relative-extinction rates", "mass extinction times"),diagnostics=c("ESS"),ess.crit=c(100,200),correction="bonferroni",xlab="million years ago",pch=19)
dev.off()


############# If multiple chains, we use different diagnostics ############################################################################

tess.plot.multiplechain.diagnostics(outputs,parameters=c("speciation rates", "extinction rates", "speciation shift times", "extinction shift times", "net-diversification rates", "relative-extinction rates", "mass extinction times"),diagnostics=c("geweke"),geweke.crit=0.05,correction="bonferroni",xlab="million years ago",pch=19)

############# If multiple chains and with Gelman-Rubin diagnostics ############################################################################

# It compares the variance in parameter estimates within a single chain with the variance between chains to determine if the parameter values can be considered to be sampled from the same probability distribution <-  Convergence accepted

# First, we read the outputs of two chains

analysis_name1 <- CoMET_noMEE_Eucalypts_2016_dated_r8s_ultrametric_binary_noOut_trueEuc_CONVERGE_RUN1
analysis_name2 <- CoMET_noMEE_Eucalypts_2016_dated_r8s_ultrametric_binary_noOut_trueEuc_CONVERGE_RUN2

output1 <- tess.process.output(analysis_name1, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

output2 <- tess.process.output(analysis_name2, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

# Then creates a composite output

outputs <- list(output1, output2)

# Applies the diagnostics

tess.plot.multichain.diagnostics(outputs, parameters=c("speciation rates", "speciation shift times", "extinction rates", "extinction shift times"," net-diversification rates", "relative-extinction rates", "mass extinction times"), diagnostics="Gelman-Rubin", gelman.crit=1.05, xlab="million years ago", col=NULL, xaxt="n", yaxt="s", pch=19))


############################################################################################################################################# 
###########                           PLOT RESULTS IN VIGNETTES                                                                     #########
#############################################################################################################################################

### Results will be plotted as vignettes within a PDF: nrows = 2, ncolumns= FIG_TYPES/2
### Figures will be plotted as vignettes: number of figures depends on Fig_Types. 
### If no MEEs are allowed as in the code [1] above, only 4 figures will be plotted
### If MEE (ALLOW_MASS_EXTINCTION = TRUE) as in the code [2] above, 6 figures will be plotted as vignettes.
### Rate shift times and mass extinction times are plotted together with Bayes Factors (2lnBF) significance levels 
### at 2,6,10 on the right y axis (this can be assigned a name with argument "yaxt") 
### Alternatively, 8 vignettes including everything can be plotted

################################################################ OPTIONAL!!!! ###########################################################

##################### To modify the x-axis (geological scale) to make it more detailed ###################################################

### Check "Modify-GeologicalTimeScale-TESS.R" for options
# We used the first alternative. Modify the tess.plot.output function, using script "tess.plot.output.mod.R"
# and changing the label argument with  labels <- pretty(c(0, treeAge), n = 25, min.n = 25)

##################### TO PLOT ONLY ONE FIGURE OR TWO FIGURES, MODIFY THE LAYOUT COMMAND, "NUM_FIGS" AND "FIG_TYPES" OBJECTS #################

NUM_FIGS <- 1
FIG_TYPES <- c("net-diversification rates")
pdf(sprintf("%s.pdf", analysis_name_noMEE))
layout(mat <- matrix(1:NUM_FIGS,nrow=1,ncol=NUM_FIGS))
tess.plot.output.mod(out, fig.types=FIG_TYPES, las=1)
dev.off()

NUM_FIGS <- 2
FIG_TYPES <- c("net-diversification rates","relative-extinction rates")
pdf(sprintf("%s.pdf", analysis_name_noMEE))
layout(mat <- matrix(1:NUM_FIGS,nrow=1,ncol=NUM_FIGS))
tess.plot.output.mod(out, fig.types=FIG_TYPES, las=1)
dev.off()

### Modify the X axis geological scale to make it more detailed: closer evenly spaced ticks  #####
## Adjust number of breaks or ticks as needed. Here, "n = 20" introduces 20 evenly breaks over x length.

# 1. Modify tess.plot.output
# Modify these lines:
## Comment out ## (line 35) labels <- pretty(c(0, treeAge))
## Replace by ##  (line 36)  labels <- pretty(c(0, treeAge), n = 30, min.n = 20)
## Save it as "tess.plot.output.mod.R"
# Source it to load it in the environment

source("tess.plot.output.mod.R")

# Now you can use tree.plot.output.mod as before


## To modify the Y axis to adjust it to the minimum value for speciation and extinction rates ######
# Modification to change the Y axis so that it does not start in zero (important for speciation and extinction rates if they never drop to zero)
# 1. Use the ChatGPT script "tess.plot.output.mod.ymin.R". [Works well but it does not show the complete plot with confidence intervals]
# 2. Modify the script "tess.plot.output.mod.R", changing these two lines:
# ylim <- c(0, max(quantilesSpeciation, quantilesExtinction))
# ylim <- c(0, max(quantilesThisOutput))

# TO (if the y axis lower bound is 0.4 instead of zero)
# ylim <- c(0.4, max(quantilesSpeciation, quantilesExtinction))
# ylim <- c(0.4, max(quantilesThisOutput))

### Plot the temperature curve overlaid on top of the speciation rate and extinction rate figure ###
# First, plot the main figure: speciation rates with the modified X axis and Y axis also modified as above (ylim <- c(0.4....)). Additionally, we also made the mean curve thicker using lwd = 2. We changed this line in the "tess.plot.output.mod.R" script, specifically adding lwd = 2, after type = "l":
# plot(x = plotAt, y = c(meanThisOutput[1], meanThisOutput), 
#                type = "l", lwd = 2, ylim = ylim, xaxt = xaxt, col = col[type], 
#                ylab = "rate", main = type, xlab = xlab, ...)

NUM_FIGS <- 1
FIG_TYPES <- c("speciation rates")
pdf(sprintf("%s_speciation.pdf", analysis_name_noMEE))
layout(mat <- matrix(1:NUM_FIGS,nrow=1,ncol=1))
source("tess.plot.output.mod.R")
tess.plot.output.mod(out, fig.types=FIG_TYPES, las=1)

# Second, plot the temperature. We first read the temperature vector "temp" contaning 1 million-year measures of global paleotemperatures from Scotese et al. (2021)
# WITH 59 intervals: MAX_VAR_AGE=58 (root node is 57.162 Mya). We remove the duplicated interval.

temp <- c(14.5, 13.982262, 14.1140187, 14.7085975, 15.7764077, 14.8181598, 14.8640161, 16.0979162, 16.3342222, 16.102692, 16.2516659, 16.409967, 16.5643955, 16.7215511, 17.084627, 18.468766, 18.1279569, 17.6510765, 17.3640895, 17.4204596, 17.3268091, 17.3560177, 18.8665765, 17.0659801, 18.5800289, 18.9638969, 18.6088709, 18.0644789, 16.4594411, 17.3484709, 17.9351214, 16.8980567, 17.8961749, 17.2611247, 18.7754193, 19.3008148, 18.8329007, 19.6450903, 20.1665978, 20.8042431, 21.4266484, 21.0447375, 23.0863524, 22.3505405, 22.2174331, 22.2424571, 22.9363945, 23.5788511, 24.2127006, 24.8031232, 25.1409294, 25.5099037, 25.1804596, 25.025148, 24.2789908, 24.3539958, 25.217426, 22.9353006)

par(new = TRUE)
# We need to reverse the curve, to start from the last value (left: 58 Myr = 23ºC ) and end with the first value (right: 0 Myr = 14.5 ºC)
maxTime <- 58
tempTime <- seq(from = 0, to = maxTime, length.out = length(temp))
plot(tempTime, temp, type = "l", lwd = 2, col = "black", axes = FALSE, xlab = "", ylab = "", ylim = range(temp), xlim = c(maxTime, 0))
axis(side = 4, col = "black", col.axis = "black")
mtext("Temperature (°C)", side = 4, line = 3, col = "black")
dev.off()

##### To plot as larger size PDF to get margins correct ###################
# If you get an error of "figure margins too large", 
# First remove any previous settings
graphics.off() # closes all open plotting devices
par(mfrow = c(1,1)) # reset layout

# Explicitly set margins large enough
par(mar = c(5, 5, 2, 5))  # bottom, left, top, right

# Not sure why but plotting as a png does not work.
#png("my_plot.png", width = 8, height = 6, res = 300)  # size in inches

# Second, plot the speciation rate curve. No need to split the plot device
NUM_FIGS <- 1
FIG_TYPES <- c("extinction rates")

tess.plot.output.mod(out, fig.types=FIG_TYPES, las=1)

# Third, plot the temperature curve with the correct information
par(new = TRUE)
maxTime <- 58
tempTime <- seq(from = 0, to = maxTime, length.out = length(temp))
plot(tempTime, temp, type = "l", lwd = 2, col = "black", axes = FALSE, xlab = "", ylab = "", ylim = range(temp), xlim = c(maxTime, 0))
axis(side = 4, col = "black", col.axis = "black")
mtext("Temperature (°C)", side = 4, line = 3, col = "black")

#dev.off()

######################################################### WITH NO-MEEs ###################################################

NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","net-diversification rates","relative-extinction rates", "speciation shift times","extinction shift times")

pdf(sprintf("%s.pdf",analysis_name_noMEE))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(out,fig.types=FIG_TYPES,las=1)
dev.off()

### Modify the X axis geological scale to make it more detailed: closer evenly spaced ticks  #####
## Adjust number of breaks or ticks as needed. Here, "n = 20" introduces 20 evenly breaks over x length.

# 1. Modify tess.plot.output
# Modify these lines:
## Comment out ## (line 35) labels <- pretty(c(0, treeAge))
## Replace by ##  (line 36)  labels <- pretty(c(0, treeAge), n = 30, min.n = 20)
## Save it as "tess.plot.output.mod.R"
# Source it to load it in the environment

source("tess.plot.output.mod.R")

# Now you can use tree.plot.output.mod as before

NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","net-diversification rates","relative-extinction rates", "speciation shift times","extinction shift times")

pdf(sprintf("%s.pdf", analysis_name_noMEE))
layout(mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2))
tess.plot.output.mod(out, fig.types=FIG_TYPES, las=1)
dev.off()

###### PLOT A LIMITED SET OF FIGURES #########
NUM_FIGS <- 2
FIG_TYPES <- c("speciation rates","extinction rates")

# TWO FIGURES IN ONE ROW, TWO COLUMNS
pdf(sprintf("%s.pdf", analysis_name_noMEE))
layout(mat <- matrix(1:NUM_FIGS,nrow=1,ncol=2))
tess.plot.output.mod(out, fig.types=FIG_TYPES, las=1)
dev.off()

# TWO FIGURES IN ONE COLUMN, TWO ROWS
pdf(sprintf("%s.pdf", analysis_name_noMEE))
layout(mat <- matrix(1:NUM_FIGS,nrow=2,ncol=1))
tess.plot.output.mod(out, fig.types=FIG_TYPES, las=1)
dev.off()



##### WITH MEEs

NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","net-diversification rates","relative-extinction rates", "mass extinction Bayes factors","mass extinction times")

pdf(sprintf("%s.pdf",analysis_name_MEE))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(outME,fig.types=FIG_TYPES,las=1)
dev.off()

### With NO-MEE BUT DIVERSIFIED STRATEGY

outDiv <- tess.process.output(analysis_name_Div, tree2, numExpectedRateChanges=EXPECTED_NUM_EVENTS)

NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","net-diversification rates","relative-extinction rates", "speciation shift times","extinction shift times")

pdf(sprintf("%s.pdf",analysis_name_Div))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(outDiv,fig.types=FIG_TYPES,las=1)
dev.off()


#######################################################################################################################################################


#### ALTERNATIVE: TO DO EVERYTHING AT ONCE FOLLOWING THE TUTORIAL

NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","speciation shift times","extinction rates","extinction shift times","mass extinction Bayes factors","mass extinction times")
if ( ALLOW_MASS_EXTINCTION == FALSE ) {
   NUM_FIGS <- 4
   FIG_TYPES <- c("speciation rates","extinction rates", "net-diversification rates", "relative-extinction rates")
}

pdf(sprintf("%s.pdf",analysis_name))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(out,fig.types=FIG_TYPES,las=1)
dev.off()
########################################################################################################################################################

