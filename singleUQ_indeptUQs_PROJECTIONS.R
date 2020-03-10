### Bayesian state-space model for estimating long-term trend in nest count time series data
### Original model estimation code adapted from Boyd et al. (2016) paper; code provided by C. Boyd
### Future projections and evaluation of probabilities code by S. Martin 
### with input from T. Eguchi, B. Langseth, A. Yau, J. Baker, T. Jones, R. Ahrens, Z. Siders, N. Ducharme-Barth

#setwd(save.dir)
require(ggplot2)

##==============================================================================
# SELECT MODELING OPTION BELOW depending on what model we're using (**FINAL** options only here for this script): 
##==============================================================================
# ---------------------------------------
# 1. Set up Loggerheads Model = singleUQ (using 3 time series of nest counts from Yakushima)
# ---------------------------------------
if(scenario=="Cc_Yakushima_sUQ"){
    # name of file directory to save model results into
  if (file.exists(save.dir)){                               # if file directory exists, set wd to it
    setwd(file.path(save.dir))
  } else {                                                  # if file directory doesn't exist, create it, then set wd to it
    dir.create(file.path(save.dir))
    setwd(file.path(save.dir))
  }
  rdatafile <- paste(scenario, ".RData", sep="")        # can use load(file=rdatafile) if you've already run code and want to re-load results
  file.tag <- paste(scenario,"_", Sys.Date(), sep="")   # to add to output file names below
  thedata <- thedata.loggers[2:31, c(1,6,5,7)]  # remove first row due to all NAs; put beach w/ first year of data in first data column
  data.cols <- 2:4                              # columns of Annual Females data to analyze
  data.rows <- 1:30                             # rows=Years of data to analyze; must start with a value rather than NA for at least 1 time series
  write.csv(thedata[data.rows , c(1,data.cols)], file=paste(file.tag,"_0_data_used.csv",sep=""),quote=FALSE, row.names=FALSE) # export actual data feeding into model
  modl <- "singleUQ"
  data.type <- "Annual_Females"                 # this is from when we were testing variants of Annual vs. Total Femals vs. Nests
  pop.name <- c("Inakahama", "Maehama", "Yotsuse")   # code below uses these population names, e.g., for plots
  pop.name.combo <- "N. Pacific Loggerheads - Yakushima 3 beaches"   # code below uses for plots
  remig <- remigLH
  clutch.freq <- clutch.freqLH
}

# ---------------------------------------
# 2. Set up Leatherbacks Model = singleUQ (with JM and W separate time series BUT *same* trend/process)
# 2A. Uses ** MEDIAN ** nest count estimates from Tomo's imputation model output
# 2B. Uses ** LOW ** (lower 95% CI) nest count estimates from Tomo's imputation model output
# 2C. Uses ** HIGH ** (upper 95% CI) nest count estimates from Tomo's imputation model output
# ---------------------------------------
if(scenario=="Dc_JM&W_MEDIAN_sUQ"|scenario=="Dc_JM&W_LOW_sUQ"|scenario=="Dc_JM&W_HIGH_sUQ") {
  if (file.exists(save.dir)){
    setwd(file.path(save.dir))
  } else {
    dir.create(file.path(save.dir))
    setwd(file.path(save.dir))
  }
  rdatafile <- paste(scenario, ".RData", sep="")       # can use load(file=rdatafile) if you've already run code and want to re-load results
  file.tag <- paste(scenario,"_", Sys.Date(), sep="")  # to add to output file names below
  thedata <- thedata.leathers
  data.cols <- 2:3               # columns of Annual Females data to analyze
  data.rows <- 1:17              # need to start with a value rather than NA for at least 1 time series
  write.csv(thedata[data.rows , c(1,data.cols)], file=paste(file.tag,"_0_data_used.csv",sep=""),quote=FALSE, row.names=FALSE) # export actual data feeding into model
  modl <- "singleUQ"
  data.type <- "Annual_Females"   # this is from when we were testing variants of Annual vs. Total Femals vs. Nests
  pop.name <- c("Jamursba Medi Leatherback Turtles", "Wermon Leatherback Turtles")
  pop.name.combo <- "Western Pacific Leatherback Turtles"
  remig <- remigLB            # remigration interval in years (for run sum period)
  clutch.freq <- clutch.freqLB
}
##==============================================================================




##==============================================================================
# PLOT DATA
##==============================================================================

# SINGLE TIME SERIES: 
# ---------------------------------------
if(length(pop.name)==1) {           
  png(filename=paste(file.tag,"_0", "_data_rawNests.png", sep=""), width=650, height=575, units="px")
  par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
  plot(thedata[,1], clutch.freq*thedata[ ,data.cols], pch=19, 
       ylim=range(clutch.freq*thedata[ ,data.cols], na.rm=T), xlab="Season", ylab="Nest Counts", main=pop.name, 
       cex.lab=1.5, cex.axis=1.5, cex.main=2, cex=1.75)
  lines(thedata[,1], clutch.freq*thedata[ ,data.cols], lwd=2, cex=1.5)
  dev.off()
  
  png(filename=paste(file.tag,"_1", "_data_lnFemales.png", sep=""), width=650, height=575, units="px")
  par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
  plot(thedata[,1], log(thedata[ ,data.cols]), pch=19, ylim=range(log(thedata[ ,data.cols]), na.rm=T), 
       xlab="Season", ylab="Ln(Annual Females)", main=pop.name, 
       cex.lab=1.5, cex.axis=1.5, cex.main=2, cex=1.75)
  lines(thedata[,1], log(thedata[ ,data.cols]), lwd=2, cex=1.5)
  dev.off()
}

# MULTIPLE TIME SERIES:
# ---------------------------------------
if(length(pop.name)>1) {
  for (i in 1:length(pop.name)) {
    png(filename=paste(file.tag,"_0.", i, "_data_rawNests.png", sep=""), width=650, height=575, units="px")
    par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
    plot(thedata[,1], clutch.freq*thedata[,i+1], pch=19, 
         ylim=range(clutch.freq*thedata[,i+1], na.rm=T), xlab="Season", ylab="Nest Counts", main=pop.name[i], 
         cex.lab=1.5, cex.axis=1.5, cex.main=2, cex=1.75)
    lines(thedata[,1], clutch.freq*thedata[,i+1], lwd=2, cex=1.5)
    dev.off()
  }
  
  for (i in 1:length(pop.name)) {
     png(filename=paste(file.tag,"_1.", i, "_data_lnFemales.png", sep=""), width=650, height=575, units="px")
     par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
     plot(thedata[,1], log(thedata[,i+1]), pch=19, ylim=range(log(thedata[,i+1]), na.rm=T), 
          xlab="Season", ylab="Ln(Annual Females)", main=pop.name[i], 
         cex.lab=1.5, cex.axis=1.5, cex.main=2, cex=1.75)
     lines(thedata[,1], log(thedata[,i+1]), lwd=2, cex=1.5)
     dev.off()
  }
}
##==============================================================================





##==============================================================================
# TRANSFORM DATA: into natural log (Ln) space for the model
##==============================================================================
log.data <- log(thedata[data.rows , data.cols])     # Ln(TOTAL NUMBER OF FEMALES)
data.mat <- t(log.data)    # tranpose to have rows = time series and cols = years of data
data.mat
##==============================================================================





##==============================================================================
## Model setup definitions and priors
##==============================================================================
# Model options:
# 1: single population process (singleUQ)
# 2: independent trend, independent variance, covariance
# 3: independent population processes

# Data
dat <- data.mat     # matrix with rows = different time series and cols = years of data within those time series
n.yrs <- ncol(dat)
n.timeseries <- nrow(dat)
Y <- rbind(dat, NA)   # add a row (ie time series of data) to trick jagsUI into NOT converting single time series matrix to vector
Y

# Set priors
a_mean <- 0
a_sd <- 4

u_mean <- 0
u_sd <- 0.05

q_alpha <- 8
q_beta <- 2

r_alpha <- 8
r_beta <- 2


# Note: The initial state, x0, is treated as a model parameter. 
# the prior matters for Wermon - setting mean to first JM data point skews Wermon to JM trend due to poor W data
# making it specific to each time series for independentUQ

# for singleUQ, set mean of X0 based on first time series (first data point of it) 
if(modl=="singleUQ") {
  x0_mean <- dat[1,1]     # take the first data point of first time series (JM for leathers)
  x0_sd <- 10             # we went with wide sd by testing for both JM and W for leathers; 10 works to make it super wide so as to not influence Wermon
}


# Set MCMC parameters
n.samples <- 10000   # 1000; bumped up to 5000 per CB after geweke.diag for deviance in chain 1 looked high at 2.1830; gives 10,000 across the 2 chains
mcmc.chains <- 2     # 2 is min; AY uses 3
mcmc.thin <- 50      # 50-100 is more than safe; 500 seems excessive per AY
mcmc.burn <- 5000    # 1000; also bumped up to 5000 (same as n.samples) per CB after geweke.diag for deviance in chain 1 still looked high at 2.2
samples2Save <- (mcmc.burn + n.samples) * mcmc.thin
##==============================================================================




##==============================================================================
# Model 1:  'singleUQ' ... multiple time series -> single population process (trend)
# This is the final one we are using (7/31/19), not the Independent UQ model below, but left that in for now
##==============================================================================
if(modl=="singleUQ"){
  
  # Model-specific parameters
  whichPop <- rep(1,n.timeseries) # multiple time series -> single population process
  n.states <- max(whichPop)
  
  Z <- matrix(0,n.timeseries+1,n.states+1)   # matrix with rows as n.timeseries and cols as n.states (pops)
  Z[n.timeseries+1, ] <- NA                  # add a row of NAs to keep jagsUI from converting single time series matrix into vector
  Z[ , n.states+1] <- NA                     # add a col of NAs to keep jagsUI from converting single state matrix into vector
  for(i in 1:length(whichPop)) Z[i,whichPop[i]] <- 1
  Z
  
  
  jags.data <- list("Y","n.yrs","n.timeseries","Z","a_mean","a_sd","u_mean","u_sd","q_alpha","q_beta","r_alpha","r_beta","x0_mean","x0_sd")
  jags.params <- c("A", "U", "Q", "R", "X0", "X")
  model.loc <- paste(main.folder, "singleUQ.txt", sep="")
  
  # Run model
  set.seed <- 132
  jags.model <- jags(jags.data, 
                     inits = NULL, 
                     parameters.to.save= jags.params, 
                     model.file=model.loc, 
                     n.chains = mcmc.chains, 
                     n.burnin = mcmc.burn*mcmc.thin, 
                     n.thin = mcmc.thin, 
                     n.iter = samples2Save, 
                     DIC = T, 
                     parallel=T, 
                     seed=set.seed)
  print(jags.model)
}
##==============================================================================


##==============================================================================
# DIAGNOSTICS: Trace plots - CB didn't like the ones from jagsUI so these were her own... 
##==============================================================================
labpar<-c("U", "Q", "R", "X0", "deviance")
xx <- 1:n.samples

outs <- array(NA, c(n.samples, 2, length(labpar)))  # order is: rows, cols, matrices: 5 matrices w/ nrows=n.samples, ncols=5)
# row, col, matrix
outs[,1,1] <- jags.model$sims.list$U[1:n.samples]
outs[,2,1] <- jags.model$sims.list$U[(n.samples+1):(n.samples*2)]
outs[,1,2] <- jags.model$sims.list$Q[1:n.samples]
outs[,2,2] <- jags.model$sims.list$Q[(n.samples+1):(n.samples*2)]
outs[,1,3] <- jags.model$sims.list$R[1:n.samples]
outs[,2,3] <- jags.model$sims.list$R[(n.samples+1):(n.samples*2)]
outs[,1,4] <- jags.model$sims.list$X0[1:n.samples]
outs[,2,4] <- jags.model$sims.list$X0[(n.samples+1):(n.samples*2)]
outs[,1,5] <- jags.model$sims.list$deviance[1:n.samples]
outs[,2,5] <- jags.model$sims.list$deviance[(n.samples+1):(n.samples*2)]

png(filename=paste(file.tag,"_3_Trace_plots.png", sep=""), width=650, height=575, units="px")
par(mfrow=c(2, 3))
for(j in 1:length(labpar)) {
  plot(xx, outs[,1,j], xlab="cycle number", ylab=labpar[j],type='b',pch=16, ylim=range(outs[,,j]))
  lines(xx, outs[,2,j], type='b',pch=16, col="gray50")
}
dev.off()
##==============================================================================



##==============================================================================
# MORE DIAGNOSTICS - export to files
##==============================================================================
McmcList <- vector("list",mcmc.chains)
for(i in 1:length(McmcList)) McmcList[[i]] <- as.mcmc(outs[,i,])

# Write model diagnostics to txt file
sink(paste(file.tag,"_4_model_diagnostics_", modl, ".txt", sep="")) # open file
print("effectiveSize(McmcList[[1]])")
print(effectiveSize(McmcList[[1]]))
cat("\n")

print("effectiveSize(McmcList[[2]])")
print(effectiveSize(McmcList[[2]]))
cat("\n")

print("geweke.diag(McmcList[[1]])")
print(geweke.diag(McmcList[[1]])) # the test statistic is the standard z score for the equality of two means - value greater than 1.96 (p<0.05) two-tailed
#cat("\n")

print("geweke.diag(McmcList[[2]])")
print(geweke.diag(McmcList[[2]]))
#cat("\n")

print("gelman.diag(McmcList)")
print(gelman.diag(McmcList)) # values of the statistic should be small for each parameter (i.e. 1.00-1.05).
sink()   # close file
##==============================================================================






##==============================================================================
# WRITE the model file"jags.model" to a txt file
##==============================================================================
file.jags <- paste(file.tag,"_2_Jags_model_", modl, ".txt", sep="") 
sink(file.jags) # open file
print(jags.model)
sink()   # close file
##==============================================================================









##==============================================================================
# PLOT MODEL FIT GRAPH: with data points, median model predicted counts, 95% CI band shading
##==============================================================================
# Define a color for confidence interval bounds on plot
col2rgb("gray", alpha=TRUE) # get specs on gray color to use in polygon and tweak transparency
mygray <- rgb(red=190, green=190, blue=190, alpha=200, maxColorValue=255)

# PLOT model fit line (estimated trend line) on top of the data  with 95% CI shading
jags.model$q50$X   # median
length(jags.model$q50$X)  # should be number of years of data 

# SINGLE TIME SERIES:
# ---------------------------------------
if(length(pop.name)==1){
  png(filename=paste(file.tag,"_6_model_fit_med95.png", sep=""), width=650, height=575, units="px")
  par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
  #ylim.specs <- c(range(log.data, na.rm=T)[1]*.95, range(log.data, na.rm=T)[2]*1.05)  # add 5% to min & max values for y axis
  y.low <- min(min(jags.model$q2.5$X), min(log.data, na.rm=T))        # min of data or CI
  y.hi <- 1.05*max(max(jags.model$q97.5$X), max(log.data, na.rm=T))  # max of data or CI
  ylim.specs <- c(y.low, y.hi)
  
  plot(thedata[data.rows,1], log.data, pch=19, ylim=ylim.specs, xlab="Year", ylab="Ln(Annual Females)", main=pop.name, 
       cex.lab=1.5, cex.axis=1.5, cex.main=2, cex=1.75)
  polygon(x= c(thedata[data.rows,1][1], thedata[data.rows,1], rev(thedata[data.rows,1])),                      # add 95% CI shading from JAGS model output
          y= c(jags.model$q2.5$X[1], jags.model$q97.5$X, rev(jags.model$q2.5$X)), 
          col=mygray, lty=0)
  lines(thedata[data.rows,1], jags.model$q50$X, col="blue", lwd=2)         # plot median fit line (better than mean for Bayes)
  points(thedata[data.rows,1], log.data, lwd=2, type="p", pch=16, cex=1.75) # add points back over CI shading
  dev.off()
}

# MULTIPLE TIME SERIES:
# ---------------------------------------
if(length(pop.name)>1){
  y.vec.low <- jags.model$q2.5$X   # create vector with lower 95% values for X 
  y.vec.hi <- jags.model$q97.5$X   # create vector with upper 95% values for X
  
  for (i in 2:length(pop.name)){   # for each time series, the X values are really X+A; put all in one vec to find min & max for plot limits
    y.vec.low <- c(y.vec.low,  jags.model$q2.5$X + jags.model$q2.5$A[i])
    y.vec.hi <- c(y.vec.hi, jags.model$q97.5$X + jags.model$q97.5$A[i])
  }
  
  y.low <- min(min(y.vec.low), min(log.data, na.rm=T))        # plot limits; min of data or CI
  y.hi <- 1.05*max(max(y.vec.hi), max(log.data, na.rm=T))     # plot limits; max of data or CI
  ylim.specs <- c(y.low, y.hi)
  
  for (i in 1:length(pop.name)){
    png(filename=paste(file.tag,"_6.", i, "_model_fit_med95.png", sep=""), width=650, height=575, units="px")
    par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
    plot(thedata[data.rows,1], log.data[data.rows,i], pch=19, ylim=ylim.specs, xlab="Year", ylab="Ln(Annual Females)", main=pop.name[i], 
         cex.lab=1.5, cex.axis=1.5, cex.main=2, cex=1.75)
   
  # Add in 'A' scaling factor from the model estimates 
    polygon(x= c(thedata[data.rows,1][1], thedata[data.rows,1], rev(thedata[data.rows,1])),                      # add 95% CI shading from JAGS model output
            y= c(jags.model$q2.5$X[1]+jags.model$q2.5$A[i], jags.model$q97.5$X+jags.model$q97.5$A[i], rev(jags.model$q2.5$X+jags.model$q2.5$A[i])), 
            col=mygray, lty=0)
    lines(thedata[,1], jags.model$q50$X+jags.model$q50$A[i], col="blue", lwd=2)         # plot median fit line (better than mean for Bayes)
    
    points(thedata[,1], log.data[,i], lwd=2, type="p", pch=16, cex=1.75) # add points back over CI shading
    dev.off()
  }
}
# Save output
saveRDS(jags.model, paste(file.tag, ".rds", sep=""))   # SLM: this is an R workspace file to re-load work thus far if needed
##==============================================================================
fy <- length(thedata[data.rows,1])       # final year of observed data; 22 for loggerheads; 17 for leatherbacks (JM)
yrf=100                                  # 100 = years into the future
nsim=10000                               # number of sim runs; length(jags.model$sims.list$X[ ,fy])  

Umed=jags.model$q50$U                    # for reporting results, median is best
Umean=jags.model$mean$U                  # but for simulations, use mean, sd in distribution
Usd=jags.model$sd$U
Uvar=Usd^2
Uci=c(jags.model$q2.5$U, jags.model$q97.5$U)

lambda.mean=exp(Umean)
lambda.med=exp(Umed)
lambda.var=exp(Uvar)
lambda.ci=exp(Uci)

Qmed=jags.model$q50$Q
Qmean=jags.model$mean$Q
Qsd=jags.model$sd$Q

X.len <- length(jags.model$sims.list$X[,fy-0])
X.thin <- seq(from=1, to=X.len, by = X.len/nsim)  # make length of the X estimate vectors same as nsim, (thin by every other value from MCMC samples)

# **SINGLE** time series
# ---------------------------------------
if(length(pop.name)==1){
  X.fym0 <- exp(jags.model$sims.list$X[,fy-0][X.thin])  # Number of Annual Females in final data year
  X.fym1 <- exp(jags.model$sims.list$X[,fy-1][X.thin])  # final data year minus one
  X.fym2 <- exp(jags.model$sims.list$X[,fy-2][X.thin])  
  X.fym3 <- exp(jags.model$sims.list$X[,fy-3][X.thin])
}

# **MULTIPLE TIME SERIES** (e.g., for leatherbacks, need to account for X representing JM since A=0 and X+A representing W since A!=0)
# e.g., Annual Females estimate for JM & W combined = exp(X) + exp(X+A), where first part is JM and second is for W.
# ---------------------------------------
if(length(pop.name)>1){
  X.fym0 <- exp(jags.model$sims.list$X[,fy-0][X.thin])  # Number of Annual Females in final data year
  X.fym1 <- exp(jags.model$sims.list$X[,fy-1][X.thin])  # final data year minus one
  X.fym2 <- exp(jags.model$sims.list$X[,fy-2][X.thin])
  X.fym3 <- exp(jags.model$sims.list$X[,fy-3][X.thin])

  for (i in 2:length(pop.name)){
    X.fym0 <- X.fym0 + exp(jags.model$sims.list$X[,fy-0][X.thin] + jags.model$sims.list$A[,i][X.thin])  # rescaling to correct magnitude
    X.fym1 <- X.fym1 + exp(jags.model$sims.list$X[,fy-1][X.thin] + jags.model$sims.list$A[,i][X.thin])
    X.fym2 <- X.fym2 + exp(jags.model$sims.list$X[,fy-2][X.thin] + jags.model$sims.list$A[,i][X.thin])
    X.fym3 <- X.fym3 + exp(jags.model$sims.list$X[,fy-3][X.thin] + jags.model$sims.list$A[,i][X.thin])
  }
}

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## plot the fit

X.total <- apply(cbind(jags.model$sims.list$X0,jags.model$sims.list$X), 2, function(v) rowSums(apply(jags.model$sims.list$A[,], 2, function(x) exp(v+x))))

X0.total <- rowSums(apply(jags.model$sims.list$A[,], 2, function(x) exp(jags.model$sims.list$X0+x)))

X.total.med <- apply(log(X.total),2,median)

# X.fit <- apply(sapply(jags.model$sims.list$U, function(x) {x * seq(0,nrow(thedata))}), 1, function(x) x + log(X0.total))

X.q <- apply(log(X.total), 2, quantile, probs=c(0.025, 0.5, 0.975))
#X0
den.X0 <- density(log(X0.total),adj=2)
den.X0$y2 <- (den.X0$y/max(den.X0$y))
q.X0 <- quantile(log(X0.total), probs=c(0.025,0.975))
xid<- sapply(q.X0, function(x) {which.min(abs(x-den.X0$x))})
#N0
den.N0 <- density(log(X.total[,ncol(X.total)]),adj=2)
den.N0$y2 <- (den.N0$y/max(den.N0$y))
q.N0 <- quantile(log(X.total[,ncol(X.total)]), probs=c(0.025,0.975))
xid.n0<- sapply(q.N0, function(x) {which.min(abs(x-den.N0$x))})

yrange <- range(c(log(rowSums(thedata[,2:ncol(thedata)], na.rm=T)), den.X0$x))

png(filename=paste(file.tag,"_model_fit_U.png", sep=""), width=7, height=4.62, units="in", res=300)
  layout(matrix(1:2,1,2),width=c(1,0.23))
  par(mar=c(4,4,1,1))

  plot(thedata[,1], log(rowSums(thedata[,2:ncol(thedata)], na.rm=T)), pch=16, type="n", las=1, ylab="log(Annual Nesters)", xlab="Season", ylim=range(pretty(yrange)), xlim=c(min(thedata[,1])-1.1, max(thedata[,1]+0.25)), xaxs='i')
  polygon(c(c(min(thedata[,1])-1,thedata[,1]), rev(c(min(thedata[,1])-1,thedata[,1]))),
          c(X.q[1,], rev(X.q[3,])),
          col='grey85',
          border=FALSE)
  
  lines(c(min(thedata[,1])-1,thedata[,1]), X.q[2,], lwd=3, col='gray50')

  polygon(x = c(rep((thedata[1,1]-1), length(xid[1]:xid[2])),rev(den.X0$y2[xid[1]:xid[2]]+(thedata[1,1]-1))), 
          y = c(den.X0$x[xid[1]:xid[2]],rev(den.X0$x[xid[1]:xid[2]])), 
          lwd=2, border=FALSE, 
          col=col2rgbA("dodgerblue3", 0.3))

  polygon(x = c(rep(thedata[nrow(thedata),1], length(xid.n0[1]:xid.n0[2])),rev(-den.N0$y2[xid.n0[1]:xid.n0[2]]+thedata[nrow(thedata),1])), 
          y = c(den.N0$x[xid.n0[1]:xid.n0[2]],rev(den.N0$x[xid.n0[1]:xid.n0[2]])), 
          lwd=2, border=FALSE, 
          col=col2rgbA("darkorchid3", 0.3))

  lines(den.X0$y2[xid[1]:xid[2]]+(thedata[1,1]-1), den.X0$x[xid[1]:xid[2]], lwd=2, col='dodgerblue3')
  lines(-den.N0$y2[xid.n0[1]:xid.n0[2]]+thedata[nrow(thedata),1], den.N0$x[xid.n0[1]:xid.n0[2]], lwd=2, col='darkorchid3')

  points(thedata[,1], log(rowSums(thedata[,2:ncol(thedata)], na.rm=T)), pch=16)
  points(c(min(thedata[,1])-1,thedata[,1]), X.total.med, pch=16, col=c('dodgerblue3',rep('red',length(X.total)-2), "darkorchid3"))

  par(mar=c(0,0,0,0))
  plot.new()
  legend("center", legend=c(expression(sum(N[list(obs,j)],j,"")), expression(sum(T[j]+a[j],j,"")), "Median r", "95% r", expression(paste(T[0]," (95%CI)")), expression(paste(N[final]," (95% CI)"))), pch=c(16,16,NA,15,NA,NA), lwd=c(NA,NA,3,NA,2,2), pt.cex=c(1,1,NA,3,NA,NA), col=c("black","red","gray50","gray85","dodgerblue3","darkorchid3"), bty="n", y.intersp = 0.9, xpd=NA)
  pop.name.combo2 <- gsub("Pacific ","Pac.\n",pop.name.combo)
  if(grepl("Leatherback",pop.name.combo2)){
    pop.name.combo2 <- gsub("Western","W.",pop.name.combo2)
    pop.name.combo2 <- gsub(" Turtles","\n(JM & W)",pop.name.combo2)
  }else{
    pop.name.combo2 <- gsub(" 3 beaches","\n(3 beaches)",pop.name.combo2)
    pop.name.combo2 <- gsub(" - ","\n",pop.name.combo2)
  }
   
  fig.lab(pop.name.combo2, xscale=0.5, yscale=0.85, cex=1)
dev.off()

# POSTERIORS for Zach to input into the projections and take model component now (7/24/19)
posts.out <- cbind("U"=jags.model$sims.list$U, "Q"=jags.model$sims.list$Q, "N_fym0"=X.fym0, "N_fym1"=X.fym1, "N_fym2"=X.fym2, "N_fym3"=X.fym3)
head(posts.out)
write.csv(posts.out, paste("Posteriors_U_and_AbundFinalYrs_", file.tag, ".csv", sep=""), row.names=FALSE)



