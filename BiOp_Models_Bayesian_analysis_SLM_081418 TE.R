### N Pacific loggerhead - Bayesian analysis

rm(list=ls())  # clear working space

## Libraries
library(jagsUI)
library(coda)

##==============================================================================
# Data

# setwd("C:/Users/charlotte.boyd/google drive/nrc/retrospective analysis")
#setwd("C:/Users/Summer.Martin/Dropbox/SLM documents/SLM-8-Turtle_Population_Models_and_Trends/Models for BiOp SSLL 2018")
thedata <- read.table("data/Nests_Japan_loggerheads_STAJ_data.txt", header=T, sep="\t")
thedata

# Add column for Total Nest Count data using Run-sum (3 year based on remigration interval)
# can skip this part if want to use Annual Nest Counts instead of Total Nest Counts for data
thedata$TotalNests <- NA
for (y in 3:length(thedata$Nests)) thedata$TotalNests[y] <- sum(thedata$Nests[(y-2):y]) # 3 year remigration interval, so use 3 year run sum
thedata
thedata <- thedata[-c(1,2), -2]  # remove first two rows (due to 3 yr run sum) & remove 2nd column (Annual Nests)
thedata

# convert Annual Nest Count data into Annual Nesting Females using clutch frequency
# doesn't matter if you convert to nesters up front here or not since just dividing all data by clutch freq of 3
data.type <- "runsum_totalFemales"
log.data <- log(1/3 * thedata[ , 2])     # Ln(TOTAL NUMBER OF FEMALES): assumes mean clutch frequency of 3 per year
ylab.txt <- "Ln(Total Females)"  # or could be "Ln(Annual Nesting Females)" if not using run sums
main.txt <- "Total Females"      # or could be "Annual Nesting Females" if not using run sums
par(mfrow=c(1,1))
plot(thedata$Year, log.data, pch=19, ylim=range(log.data, na.rm=T), xlab="Year", ylab=ylab.txt, main=main.txt)
lines(thedata$Year, log.data, lwd=2)

npac.data <- matrix(log.data, nrow=1)   # time series data in 1 matrix row

##==============================================================================
## Model summary: single trajectory

# Data
dat <- npac.data
#n.yrs <- ncol(dat)
n.timeseries <- nrow(dat)
Y <- as.vector(dat)
n.yrs <- length(Y)

# Add non-observed future data:
Y.future <- rep(NA, 100)
#Y <- c(Y, Y.future)
n.yrs.proj <- length(Y.future)

# Set priors
a_mean <- 0
a_sd <- 4

#u_mean <- 0
#u_sd <- 0.5

q_alpha <- 0.01
q_beta <- 0.01

r_alpha <- 0.01
r_beta <- 0.01

x0_mean <- dat[1] # Note: The initial state, x0, is treated as a model parameter.
x0_sd <- 1

#-------------------------------------------------------------------------------
# Run model

# Model-specific parameters
whichPop <- rep(1,n.timeseries) # multiple time series -> single population process
n.states <- max(whichPop)

Z <- matrix(0,n.timeseries,n.states)
for(i in 1:length(whichPop)) Z[i,whichPop[i]] <- 1
Z

# Set MCMC parameters
n.samples <- 10000   # 1000; bumped up to 5000 per CB after geweke.diag for deviance in chain 1 looked high at 2.1830; gives 10,000 across the 2 chains
mcmc.chains <- 2     # 2 is min; AY uses 3
mcmc.thin <- 50      # 50-100 is more than safe; 500 seems excessive per AY
mcmc.burn <- 1000    # 1000; also bumped up to 5000 (same as n.samples) per CB after geweke.diag for deviance in chain 1 still looked high at 2.2
samples2Save <- (mcmc.burn + n.samples) * mcmc.thin  # why multiply by mcmc.thin? TE

#                  "u_mean",
#"u_sd",

jags.data <- list("Y",
                  "n.yrs",
                  "n.yrs.proj",
                  "n.timeseries",
                  "Z",
                  "a_mean",
                  "a_sd",
                  "q_alpha",
                  "q_beta",
                  "r_alpha",
                  "r_beta",
                  "x0_mean",
                  "x0_sd")

jags.params <- c("A", "U", "Q", "R", "X0", "X", "u_mean", "u_sd", "mean.u")
model.loc <- "singleUQ_singletimeseries1_TE.txt"   # paste("standard model files/singleUQ_singletimeseries1.txt",sep="")

# Run model
set.seed(246)
jags.model <- jags(jags.data,
                   inits = NULL,
                   parameters.to.save= jags.params,
                   model.file=model.loc,
                   n.chains = mcmc.chains,
                   n.burnin = mcmc.burn*mcmc.thin,
                   n.thin = mcmc.thin,
                   n.iter = samples2Save,
                   DIC = T, parallel=T)
print(jags.model)

# plot the model fit and future projection (median + 95% credible interval)
plot(1:length(jags.model$q50$X), jags.model$q50$X, ylim=c(0,21), pch=16, cex=.7)
lines(1:length(jags.model$q50$X), jags.model$q50$X, col="blue")
lines(1:length(jags.model$q2.5$X), jags.model$q2.5$X, col="red")
lines(1:length(jags.model$q97.5$X), jags.model$q97.5$X, col="red")

hist(jags.model$sims.list$u_mean)
hist(jags.model$sims.list$U)
#-------------------------------------------------------------------------------
# Trace plots - CB didn't like the ones from jagsUI so does her own...

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

par(mfrow=c(2, 3))
for(j in 1:length(labpar)) {
  plot(xx, outs[,1,j], xlab="cycle number", ylab=labpar[j],type='b',pch=16, ylim=range(outs[,,j]))
  lines(xx, outs[,2,j], type='b',pch=16, col="gray50")
}

#-------------------------------------------------------------------------------
# Diagnostics

McmcList <- vector("list",mcmc.chains)
for(i in 1:length(McmcList)) McmcList[[i]] <- as.mcmc(outs[,i,])

effectiveSize(McmcList[[1]])
# var1     var2     var3     var4     var5
# 10263.78 10375.48 10000.00  9793.11 10000.00

effectiveSize(McmcList[[2]])
# var1      var2      var3      var4      var5
# 9433.693 10000.000 10000.000  9693.893 10000.000

geweke.diag(McmcList[[1]]) # the test statistic is the standard z score for the equality of two means - value greater than 1.96 (p<0.05) two-tailed
# var1    var2    var3    var4    var5
# -0.5666  0.2067 -2.3968  1.7167 -2.0420

geweke.diag(McmcList[[2]])
# var1    var2    var3    var4    var5
# 0.3448  0.1055 -0.4066 -0.4308 -0.9120

gelman.diag(McmcList) # values of the statistic should be small for each parameter (i.e. 1.00-1.05).
# Potential scale reduction factors:
#
#   Point est. Upper C.I.
# [1,]          1          1
# [2,]          1          1
# [3,]          1          1
# [4,]          1          1
# [5,]          1          1
#
# Multivariate psrf
#
# 1


#-------------------------------------------------------------------------------
# Plot U and Q

par(mfrow=c(1,2))
hist(jags.model$sims.list$U, col="grey", breaks=100, xlab="", main="U")
abline(v=0, lwd=2, lty=2)
hist(jags.model$sims.list$Q, col="grey", breaks=100, xlab="", main="Q")

#-------------------------------------------------------------------------------
# Save output

saveRDS(jags.model, "XXX.rds")   # SLM: not sure what this file type is...

# RDS files are conveinient because you can "load" them to different names,
# unlike load and save



# ----------------------------------------------------------------------------------------------------------------------
# FOR PROJECTION: Forecast using simulation; calculate quantiles across simulation runs to determine how many runs fall above certain quantiles
fy <- 22                         # final year of observed data
yrf=100                          # 100 = years into the future
nsims=10000                      # 50 = number of simulation runs

rmean=jags.model$mean$U
rsd=jags.model$sd$U
rvar=rsd^2
rci=c(jags.model$q2.5$U, jags.model$q97.5$U)

lambda.mean=exp(rmean)
lambda.ci=exp(rci)

nlast=jags.model$mean$X[fy]      #  = mean Total Females estimate for final year of observed data (starting point for simulation)
sdnlast=jags.model$sd$X[fy]    #  = sd Total Females for that last year of data

pN.rstat=matrix(nrow=yrf+1,ncol=nsims)   # +1 for yr=0; for static 'r' projection; predicted N (nests) = matrix with YEARS as rows and SIMULATIONS as columns
pN.rdyn=matrix(nrow=yrf+1,ncol=nsims)    # +1 for yr=0; for dynamic 'r' projection; predicted N (nests) = matrix with YEARS as rows and SIMULATIONS as columns
dim(pN.rstat)
# ----------------------------------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------
# FUTURE PROJECTION SIMULATION for PVA -- use *either* simulated or real data from above
# Random draw from dist of r for year 1 of sim *AND* each subsequent year 2:yrf
# reminder: here pN is matrix with YEARS as rows and SIMULATIONS as columns and NUMBER OF TOTAL FEMALES as content
for(sim in 1:nsims)
{
  N0 <- exp(rnorm(1,nlast,sdnlast))   # for both rdynamic & rstatic approaches; start year 0 of sim with random draw from N
  if(N0>=0) {
    pN.rdyn[1,sim]= N0
    pN.rstat[1,sim]= N0
  }
  if(N0<0) {
    pN.rdyn[1,sim] <- 0    # if pop starting point is below 0, use 0 rather than go negative (not biologically feasible)
    pN.rstat[1,sim] <- 0
  }

  rstat.sim=rnorm(1,rmean,rsd)            # for r dynamic, a single constant 'r' to carry through all future years of ONE simulation run

  for(yr in 1:yrf)                   # for each year of sim...
  {
    i <- yr + 1     # adjust row index ... first row is year 0, second row is year 1, etc.
    pN.rdyn[i,sim]=pN.rdyn[i-1,sim]*exp(rnorm(1,rmean,rsd))   # (N from previous year) * exp(draw r from dist)
    if(pN.rdyn[i,sim]<0) {pN.rdyn[i,sim] <- 0}                # if pop falls below 0, cut it off at 0 rather than go negative

    pN.rstat[i,sim]=pN.rstat[i-1,sim]*exp(rstat.sim)          # (N from previous year) * exp(constant r for each sim run into future)
    if(pN.rstat[i,sim]<0) {pN.rstat[i,sim] <- 0}              # if pop falls below 0, cut it off at 0 rather than go negative
  }
}

# view final year of simulation; sim left off on 10000 here if nsim=10000
# look at same simulation run for the two approaches... STATIC r gives WAY HIGHER UNCERTAINTY!!
mean(pN.rstat[ ,100])
mean(pN.rdyn[ ,100])


# calculate quantiles & means for simulations w/ DYNAMIC 'r' pulled from dist each year into future
# reminder: pN is matrix with YEARS as rows and SIMULATIONS as columns
quants.rdyn=t(apply(X=pN.rdyn,MARGIN=1,FUN=quantile,probs=c(0.025,0.5,0.975)))    # median + CI for each year; MARGIN=1 indicates rows for a matrix;
means.rdyn=t(apply(X=pN.rdyn,MARGIN=1,FUN=mean))             # mean for each year; MARGIN=1 is for rows
quants.rdyn.df <- as.data.frame(quants.rdyn)
quants.rdyn.df$Year <- 0:yrf
quants.rdyn.df$Mean <- as.vector(means.rdyn)
quants.rdyn.df <- quants.rdyn.df[ ,c("Year", "Mean", "2.5%", "50%", "97.5%")]  # reorder cols; rows=YEARS, cols=quantiles of sim runs (2.5%, 50%, 97.5%)
quants.rdyn.df

# calculate quantiles & means for simulations w/ STATIC 'r' pulled ONCE from dist and carried as constant in all future years
quants.rstat=t(apply(X=pN.rstat,MARGIN=1,FUN=quantile,probs=c(0.025,0.5,0.975)))    # median + CI for each year; MARGIN=1 indicates rows for a matrix;
means.rstat=t(apply(X=pN.rstat,MARGIN=1,FUN=mean))             # mean for each year; MARGIN=1 is for rows
quants.rstat.df <- as.data.frame(quants.rstat)
quants.rstat.df$Year <- 0:yrf
quants.rstat.df$Mean <- as.vector(means.rstat)
quants.rstat.df <- quants.rstat.df[ ,c("Year", "Mean", "2.5%", "50%",  "97.5%")]  # reorder cols; rows=YEARS, cols=quantiles of sim runs (2.5%, 50%, 97.5%)
quants.rstat.df

# if high uncertainty in 'r' leads to NEGATIVE turtle/nest values, that's biologically not feasible [log(.99)=-0.01 and log(-1)=NaN],
# so change projected value to 1 so log(1)=0 or else simulation stuff below won't work
# quants.rstat.df[quants.rstat.df < 1] <- 1
# quants.rdyn.df[quants.rdyn.df < 1] <- 1
# ----------------------------------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------
# plot SIMULATION QUANTILES = MEDIAN PROJECTION line plus shaded 2.5% and 97.5% lines to show distribution of outcomes
par(mfrow=c(1,1))
for (p in 1:2) {
  ocol=rgb(253,106,2,alpha=40,max=255)
  # 'r' static: plot from simulation: median predicted trajectory w/ 2.5% and 97.5% intervals (at each year, takes median, 2.5% & 97.5% quantiles of sim runs)
  if (p==1){
    matplot(x=0:yrf,y=log(quants.rstat.df[,3:5]), type="l",lty=c(2,1,2),lwd=c(1,2,1),col="black", xlab="Year", ylab=ylab.txt)
    polygon(c(0:yrf,yrf:0),c(log(quants.rstat.df[,3]),log(rev(quants.rstat.df[,5]))),col=ocol,border=FALSE)
    title(main="Simulation projection with static r")
  }
  # 'r' dynamic: plot from sim'n: median predicted trajectory w/ 2.5% and 97.5% intervals (at each year, takes median, 2.5% & 97.5% quantiles of sim runs)
  if (p==2) {
    matplot(x=0:yrf,y=log(quants.rdyn.df[,3:5]), type="l",lty=c(2,1,2),lwd=c(1,2,1),col="black", xlab="Year", ylab=ylab.txt)
    polygon(c(0:yrf,yrf:0),c(log(quants.rdyn.df[,3]),log(rev(quants.rdyn.df[,5]))),col=ocol,border=FALSE)
    title(main="Simulation projection with dynamic r")
  }
}
# ----------------------------------------------------------------------------------------------------------------------




# ----------------------------------------------------------------------------------------------------------------------
# DEFINE plotting FUNCTION to show all sim runs + abundance thresholds
# plot SIMULATION PROJECTION POINTS predicted by EACH SIMULATION RUN at each year in the future to show distribution of outcomes
# shows all simulation runs
require(ggplot2)
plot.simruns <- function(r.sel="rstat", thresh.lines=FALSE) {    # choose rdyn or rstat and whether to add lines for thresholds
  if(r.sel=="rstat") {pN <- pN.rstat; g.title <- "Simulation projection with static r"}     # 'r' static: runs from simulation with r static for each run
  if(r.sel=="rdyn")  {pN <- pN.rdyn; g.title <- "Simulation projection with dynamic r"}
  #head(as.vector(pN.))            # creates the vector going column by column (i.e., sim 1 yr by yr, then sim 2 yr by yr, etc.)
  rep4plot=cbind(rep(1:nrow(pN),times=ncol(pN)),as.vector(pN))  # Zach's fix to Rob's code (same as my guess below)
  colnames(rep4plot)=c("x","y")   # x = Year; y = Nests
  x=rep4plot[,1]
  y=log(rep4plot[,2])   # can change here to log(Nests) or regular exponential curve of Nests
  df <- data.frame(x = x, y = y, d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
  p <- ggplot(df) +
    geom_point(aes(x, y, col = d), size = 1) +
    scale_color_identity() +
    theme_bw() +
    labs(x = "Year",y=ylab.txt, title=g.title) +
    theme(axis.title.x = element_text(color = "blue", size = 16)) +
    theme(axis.title.y = element_text(color = "blue", size = 16))

  # add mean & 95% CI abundance thresholds as horizontal lines for 50%, 25%, and 12.5% abund.
  if(thresh.lines==TRUE) {
    thresh.vals50 <- 0.50*quantile(log(CurTotFem), probs=c(0.025, 0.50, 0.975))   # CurTotFem must be defined outside the function
    thresh.vals25 <- 0.25*quantile(log(CurTotFem), probs=c(0.025, 0.50, 0.975))
    thresh.vals12 <- 0.125*quantile(log(CurTotFem), probs=c(0.025, 0.50, 0.975))
    p <- p +  geom_hline(yintercept=thresh.vals50[1], colour="black", linetype="dashed") +  # 50% abund thresh lower 95% CI
      geom_hline(yintercept=thresh.vals50[2], colour="black", linetype="solid") +  # median
      geom_hline(yintercept=thresh.vals50[3], colour="black", linetype="dashed") +    # upper 95% CI
      geom_hline(yintercept=thresh.vals25[1], colour="black", linetype="dashed") +  # 25% abund thresh lower 95% CI
      geom_hline(yintercept=thresh.vals25[2], colour="black", linetype="solid") +  # median
      geom_hline(yintercept=thresh.vals25[3], colour="black", linetype="dashed") +    # upper 95% CI
      geom_hline(yintercept=thresh.vals12[1], colour="black", linetype="dashed") +  # 12.5% abund thresh lower 95% CI
      geom_hline(yintercept=thresh.vals12[2], colour="black", linetype="solid") +  # median
      geom_hline(yintercept=thresh.vals12[3], colour="black", linetype="dashed")    # upper 95% CI
  }
  print(p)
  if(r.sel=="rstat") p.rstat <- p    # save graph images outside of function
  if(r.sel=="rdyn") p.rdyn <- p
}
plot.simruns(r.sel="rstat", thresh.lines=FALSE)
plot.simruns(r.sel="rdyn", thresh.lines=FALSE)

# ----------------------------------------------------------------------------------------------------------------------







## ------------------------------------------------------------------------------------------------------------------
## EVALUATION PROBABILITIES OF REACHING ABUNDANCE THRESHOLDS:

# For a given year in the future (5, 10, 25, 50, 100)
# Calculate probability of reaching abundance thresholds
# As the proportion of simulation runs falling below the threshold for that year
# PIRO wants 95% CI associated with this probability.... currently do this as bootstrap of simulation runs for each evaluation year

# select whether to use the projections that used the 'r dynamic' or 'r static' approach
#rsel <- "rstat"
rsel <- "rdyn"

# based on selection above, pick the corresponding projection matrix
if(rsel=="rdyn")  pN.sel <- pN.rdyn
if(rsel=="rstat") pN.sel <- pN.rstat

# clutch frequency and remigration interval numbers factor in below if need to convert Annual Nests or Total Nests to Total Females
clutch.freq <- 3   # 3 nests per female in a nesting year
remig <- 3         # 3 year remigration interval = 3 year run sum period

# calculate "current" abundance
# can simply use N0 from simulation if data type is Total Females (already factors in run sum total earlier in model)
N0 <- pN.sel[1,]   #  first row of sim matrix (same for rdyn & rstat) is Year 0 draw of N in 2016 from dist (final data year for loggers)

# calculate "current" abundance using 3 year run sum of Annual Nests or Annual Nesters
if (data.type=="annual_Nests"|data.type=="annual_Nesters") {
  Nm1 <- rnorm(10000,Nmean.m1,Nsd.m1)    # 10,000 sim runs pulling random draw from "Final data year minus 1" = 2015 for loggers
  Nm1[Nm1<0] <- 0                        # if pop starting point is below 0, use 0 rather than go negative (not biologically feasible)
  Nm2 <- rnorm(10000,Nmean.m2,Nsd.m2)    # 10,000 sim runs pulling random draw from "Final data year minus 2" = 2014 for loggers
  Nm2[Nm2<0] <- 0                        # if pop starting point is below 0, use 0 rather than go negative (not biologically feasible)
}

# calculate Current TOTAL FEMALES from last data year (best estimate for current abundance for baseline)
# vector with length nsim -- different starting point of Current TOTAL FEMALES for each sim
if(data.type=="annual_Nests") CurTotFem <- (N0 + Nm1 + Nm2)/clutch.freq        # if data are Annual Nests, 3-yr run sum divided by 3 nests/female
if(data.type=="runsum_totalNests") CurTotFem  <- N0/clutch.freq                # if data are of TOTAL Nests (data as 3-yr Nest runsums), divide by 3 nests/female
if(data.type=="runsum_totalFemales") CurTotFem <- N0                 # if projections are of TOTAL Females (data as 3-yr Female runsums), just take final
length(CurTotFem)                              # should be number of sims

# Plot sim runs with abundance thresholds
plot.simruns(r.sel=rsel, thresh.lines=TRUE)
#plot.simruns(r.sel=rsel, thresh.lines=FALSE)  # plot without thresholds

# select here whether to use pN.rdyn or pN.rstat projected abundance matrix
dim(pN.sel)              # 101 rows because first row is N0, final data year (year 0) which is starting point for simulation years 1-100
popsize <- t(pN.sel[-1, ])   # transpose to get 10,000 sim rows and 100 year columns (same format as Boyd's Bayesian code); remove that first row (YEAR 0)
dim(popsize)              # now in proper format: rows = # sims (10,000) and cols = # years in future (1-100)

# set up array of 1 matrix per threshold to store 1 if pop is above a threshold, 0 if at or below it
nsim <- nsims
tmax <- yrf
pop.aboveT <- array(1, c(nsim, tmax, 3))  # order is: rows=sims, cols=yrs, matrices: 3 threshold matrices w/ nrows=10000 sims, ncols=100 yrs)
dim(pop.aboveT)

# Need projections in terms of Projected Total Females in the population to compare to prescribed abundance thresholds of interest
# create empty matrix ProjTotFem that has rows = # sims (10,000) and cols = # years in future (1-100)
ProjTotFem <- matrix(NA, nrow=nsim, ncol=tmax)

# define start year for projection time series
if(data.type=="annual_Nests") sy <- remig       # if projections are of Annual Nests, start on future year 3 to compute run sum backward
if(data.type=="runsum_totalNests") sy <- 1  # if projections are of TOTAL Nests (already runsummed), start on future year 1
if(data.type=="runsum_totalFemales") sy <- 1 # if projections are of TOTAL Females (already runsummed), start on future year 1

# Fill in matrix of PROJECTED TOTAL FEMALES for each year to prepare for compare each to ABUNDANCE THRESHOLDS
for (i in 1:nsim) {
  for(tt in sy:tmax) {                   # start at year 3 since projections are of ANNUAL NESTS (not a run sum of total abund each year)
    if(data.type=="annual_Nests") ProjTotFem[i,tt] <- sum(popsize[i,((tt-(remig-1)):tt)])/clutch.freq   # divide 3-yr Summed Annual Nests by clutch freq
    if(data.type=="runsum_totalNests") ProjTotFem[i,tt] <- popsize[i,tt]/clutch.freq        # divide Total Nests (projected as such) by clutch freq
    if(data.type=="runsum_totalFemales") ProjTotFem[i,tt] <- popsize[i,tt]    # no need to divide; already projected as TOTAL not Annual (no need to run sum)
  }
}
colMeans(ProjTotFem)  # view row 1 = sim 1 for all 100 years into future
plot(x=1:tmax, y=colMeans(ProjTotFem, na.rm=TRUE), ylab="Mean of Projections")
plot(x=1:tmax, y=apply(ProjTotFem, MARGIN=2, FUN=quantile, probs=0.5, na.rm=TRUE), ylab="Median of Projections")

# define abundance thresholds of interest: e.g., 50% of current abundance (total females in 2016 or sum annual nests or nesters from 2014-2016)
# Create an ABUNDANCE THRESHOLD matrix specific to each simulation run
# Each column is for a different abundance threshold (50% of current, 25%, or 12.5%)
ThreshTotFem <- matrix(NA, nrow=nsim, ncol=3)
ThreshTotFem[ ,1] <- 0.5*CurTotFem    # 50% thresh in column 1
ThreshTotFem[ ,2] <- 0.25*CurTotFem   # 25% thresh in column 2
ThreshTotFem[ ,3] <- 0.125*CurTotFem  # 12.5% thresh in column 3

# To Compare ProjTotFem to Specified Abundance thresholds (specific to each simulation run)
# go through each threshold "th", each sim "i", and projected year "y"
# leave as 1 if ProjTotFem is > ThreshTotFem and update to 0 if ProjTotFem <= ThreshTotFem
for(th in 1:3){
  for (i in 1:nsim) {
    for (y in 1:tmax) {
      if(!is.na(ProjTotFem[i, y]) & ProjTotFem[i, y] <= ThreshTotFem[i, th]) {pop.aboveT[i,y,th] <- 0}   # row i = sim num, col y = fut year, th = thresh mat
    }
  }
}



# Go through pop.aboveT output for each abundance threshold
# Check if any sim runs *END* below the threshold (disregard runs in which pop dips below thresh but recovers to end above thresh)
# For sim runs that *END* below threshold, calculate the years until it drops and stays below thresh
# PIRO (A) & (B): Calculate mean, median, and 95% CI limits for "Years to Threshold" using all the sims that *END* below threshold
# PIRO request (C): estimate the probability of the pop reaching those thresholds (50%, 25%, 12.5% of current abundance)
# in 5, 10, 25, 50, and 100 year time intervals with associated 95% confidence intervals

TIMEtoTHRESH <- data.frame(matrix(nrow=3, ncol=6), row.names=c("50% abund", "25% abund", "12.5% abund"))
names(TIMEtoTHRESH) <- c("probEndAbove", "probEndBelow", "MeanYrsToThresh", "2.5%", "50%", "97.5%")
TIMEtoTHRESH

eval.yrs <- c(5,10,25,50,100)   # year in future at which to evaluate prob of reaching the abundance thresholds
PROBatYEAR <- data.frame(matrix(nrow=3, ncol=15), row.names=c("50% abund", "25% abund", "12.5% abund"))
names(PROBatYEAR) <- c("Y5", "Y5l", "Y5u", "Y10", "Y10l", "Y10u", "Y25", "Y25l", "Y25u",
                       "Y50", "Y50l", "Y50u", "Y100", "Y100l", "Y100u")
PROBatYEAR

# THRESHOLDS: cycle through
for (th in 1:3){
  # for a threshold matrix (1= sim/yr projection above curr abund thresh %; 0= at or below thresh)
  # for all runs ending in 1, consider them a success -- pop did not fall/stay below the thresh
  simsSums <- rowSums(pop.aboveT[,,th])          # sum each row (sim); if a run stayed above thresh whole time, sum = 100
  simsDropBelow <- which(simsSums<tmax)          # identify rows that may DROP below thresh at some point (but may still recover and end up ABOVE thresh)
  simsEndBelow <- which(pop.aboveT[,tmax,th]!=1) # identify sim runs (rows) that DON'T end in 1; ended at tmax below thresh
  simsEndAbove <- which(pop.aboveT[,tmax,th]==1) # identify sim runs (rows) that DO end in 1; ended at tmax below thresh
  simsDropBelowEndAbove <- which(simsSums<tmax & pop.aboveT[,tmax,th]==1)  # sim runs that drop below thresh but then recover to end at tmax above thresh

  # if no sim runs end below thresh, then prob of reaching thresh within 100 yrs is 0 and "Years to reach thresh" is NA
  if(length(simsEndBelow)==0) {probEndBelow <- 0; probEndAbove <- 1; YrsToThresh <- NA; YrsToThresh.mean <- NA; YrsToThresh.quants <- c(NA, NA, NA)}

  # if some sim runs end below thresh, then look to see when they fall below and stay below (max year when pop.aboveT = 1 for that sim run)
  if(length(simsEndBelow)!=0) {
    probEndBelow <- length(simsEndBelow)/nsim
    probEndAbove <- length(simsEndAbove)/nsim
    YrsToThresh <- apply(X=pop.aboveT[simsEndBelow, , th], MARGIN=1, FUN=function(x){max(which(x==1))})
    YrsToThresh.mean <- mean(YrsToThresh)
    YrsToThresh.quants <- quantile(YrsToThresh, probs=c(0.025, 0.50, 0.975))
  }

  # Package the output for each threshold
  out <- c(probEndAbove, probEndBelow, round(YrsToThresh.mean, 1), round(YrsToThresh.quants,1))
  out.df <- as.data.frame(t(out))
  names(out.df) <- c("probEndAbove", "probEndBelow", "MeanYrsToThresh", "2.5%", "50%", "97.5%")
  TIMEtoTHRESH[th,] <- out.df

  # Prob of reach the threshold at year 5, 10, 25, 50, 100
  pop.aboveT.edit <- pop.aboveT
  pop.aboveT.edit[simsDropBelowEndAbove, , th] <- 1    # change 0s to 1s for sim runs that dropped below thresh but recovered to end above

  # EVALUATION YEARS: cycle through
  for (ey in 1:length(eval.yrs)) {
    fut.yr <- eval.yrs[ey]
    simsReach <- which(pop.aboveT.edit[ , fut.yr, th]==0)  # to get prob of reaching thresh in that fut.yr (don't later recover), calc nmbr sim rows=0
    PROBatYEAR[th, 1+(ey-1)*3] <- length(simsReach)/nsim

    # Bootstrap 95% CI from the nsim results for the YEAR and THRESH
    CI = 0.95
    nboot = nsim       # set bootstrap number of runs = nsim
    x = as.vector(pop.aboveT.edit[ , fut.yr, th])    # col of specified YEAR for the threshold;
    n = length(x)                                    # should be same as nsims
    xStat = length(which(x == 0))     # calculate the statistic on original sample: Prob of runs falling below an abund threshold
    # Generate 'nboot' number of bootstrap samples, i.e. an n x nboot array of  (rows=n ; cols=nboots)
    tmpdata = sample(x=x,size=n*nboot, replace=TRUE)         # random resamples from x
    bootstrapsample = matrix(tmpdata, nrow=n, ncol=nboot)
    # Compute the Probability stat for each resample (column of bootstrapsample)
    bsStats =  apply(X=bootstrapsample, MARGIN=2, FUN=function(simvec){length(which(simvec == 0))}) # MARGIN=2 for columns of matrix
    bsStats.mean = mean(bsStats)   # mean of the bootstrapped sample statistics = the sample statistic calculated on the original sample
    # Compute delta* for each bootstrap sample [delta* = resample stat - original sample stat (xbar)]
    deltastar = bsStats - xStat    # (integer) diff in number of runs falling at/below threshold between bootstrap sample and original sample
    # Find the Lower and Upper CI quantile for deltastar
    ciU = 1-(1-CI)/2      # lower/upper confidence interval limits e.g. 0.025 for 2.5% & 0.975 for 97.5% for a 95% CI
    ciL = 0+(1-CI)/2
    d = quantile(deltastar, c(ciL, ciU))  #  the quantile() function is sophisticated about choosing a quantile between two data points.
    # Calculate the specified confidence interval for original sample stat.
    ci = xStat - c(d[2], d[1])
    ci.prob <- ci/nsim
    # Add Bootstrapped CI to the storage matrix for this threshold row in appropriate year cols
    PROBatYEAR[th, (2+(ey-1)*3): (3+(ey-1)*3)] <- ci.prob

  } # END CYCLE THROUGH EVALUTION YEARS
}  # END CYCLE THROUGH THRESHOLDS

TIMEtoTHRESH  # Answers PIRO request (A) and (B)
PROBatYEAR    # Answers PIRO request (C)

rmean         # Answers PIRO request  (D) population's mean log growth rate
rvar          # Answers PIRO request  (D) variance of mean log growth rate
rci           # Answers PIRO request  (D) 95% confidence interval

lambda.mean   # Answers PIRO request  (E) finite rate of increase, lambda
lambda.ci     # Answers PIRO request  (D) 95% confidence interval for lambda


# Start writing to an output file
sink('analysis-output.txt')
print("(A & B) -------------- ")
print("TIMEtoTHRESH")
TIMEtoTHRESH  # Answers PIRO request (A) and (B)
cat("\n\n")

print("(C) -------------- ")
print("PROBatYEAR")
round(PROBatYEAR,3)    # Answers PIRO request (C)
cat("\n\n")

print("(D) -------------- ")
print("Population's mean log growth rate, variance, and 95% confidence interval")
as.numeric(rmean)         # Answers PIRO request  (D) population's mean log growth rate
rvar          # Answers PIRO request  (D) variance of mean log growth rate
rci           # Answers PIRO request  (D) 95% confidence interval
cat("\n\n")

print("(E) -------------- ")
print("Population's finite rate of increase (lambda) and 95% confidence interval")
as.numeric(lambda.mean)   # Answers PIRO request  (E) finite rate of increase, lambda
lambda.ci     # Answers PIRO request  (E) 95% confidence interval for lambda
# Stop writing to the file
sink()



## ------------------------------------------------------------------------------------------------------------------










































# ##==============================================================================
# ## Projections
# fy <- 24         # final year of data
# 1/3 * thedata[fy,2] # 1627.667; final year data, assuming mean clutch frequency of 3 per year
# exp(median(jags.model$sims.list$X[,fy])) # 1684.103 final year median posterior prediction [model doesn't believe 2008 estimate based on complete dataset]
# quantile(exp(jags.model$sims.list$X[ ,fy]), probs=c(0.025, 0.5, 0.975))  # SM added
#
# nsim <- mcmc.chains*n.samples
# tmax <- 100
#
# # Simulate single population process
# set.seed(205)
# log.delta.N <- matrix(NA, nrow=nsim, ncol=tmax)
# for (i in 1:nsim) log.delta.N[i,] <- rnorm(tmax, jags.model$sims.list$U[i], sqrt(jags.model$sims.list$Q[i]))
#
# # Apply population process to timeseries
# logpopsizes <- matrix(NA, nsim, tmax)
# for (i in 1:nsim) logpopsizes[i,] <- jags.model$sims.list$X[i,fy] + cumsum(log.delta.N[i,])
# popsize <- exp(logpopsizes) # convert back to normal space
#
# # ### Apply test decision rules to adult females returning over 3 years
#
# # find extirpated subpopulations (subpopulation is unlikely to contribute to DPS persistence if sum N < NET over 3 consecutive years)
# NET <- 125 # near-extinction threshold = 250 mature individuals or 125 mature females
# pop.extant <- matrix(1, nsim, tmax)
# for (i in 1:nsim) for(tt in 3:tmax) if(sum(popsize[i,((tt-2):tt)]) < NET) {pop.extant[i,tt:tmax] <- 0; break}
#
# # NET = 250; 125
# # P(N < 125 mature females within 50 years)
# print(1 - sum(pop.extant[,50])/nsim)  # 0.1905
# # P(N < 125 mature females within 100 years)
# print(1 - sum(pop.extant[,100])/nsim) # 0.3025
#
#
#
#
# # ------------------------------------------------------
# # ------------------------------------------------------
# ## NEW addition to by SM to fulfill specific PIRO request:
#
# ## (a) estimate *mean* time (number of years) until pop declines to 50%, 25%, and 12.5% of current abundance estimates
# ## (b) estimate *median* time (number of years) until pop declines to 50%, 25%, and 12.5% of current abundance estimates
# ## (c) the probability of the pop reaching those thresholds (50%, 25%, and 12.5% of current abundance estimates) in 5, 10, 15, 50 & 100 yrs w/ 95% CI
#
# ## (a) APPROACH: use 'popsize' matrix from simulation (nrows=10,000 sim runs; ncols=100 years projection; content = # of nesters N showing up that year)
# #  N for each year in a sim row = jags sim run est of N in final year + delta N from jags sim run est of pop growth rate & variance (new draw each future year)
# #  for each sim run, calc years with popsize N est <= 50%*starting N from that sim
#
# runSumF <- sum(exp(jags.model$sims.list$X[1,(fy-2):fy]))
#
# # looking at jags simulation outputs to see distribution on N (annual nesters) for last year of data (2016 for loggers)
# N.med <- exp(median(jags.model$sims.list$X[,fy]))
# N.lwr <- quantile(x=exp(jags.model$sims.list$X[,fy]), probs=c(0.025,0.50,0.975))[1]
# N.upr <- quantile(x=exp(jags.model$sims.list$X[,fy]), probs=c(0.025,0.50,0.975))[3]
#
#
# curr.N <- exp(jags.model$sims.list$X[ , fy])    # 10000 jags sim output estimates of N in 2016 (final data year for loggers)
# quantile(curr.N, prob=c(0.025,0.50,0.975))
# plot(density(curr.N))
#
# thresh.vec <- c(0.5, 0.25, 0.125)         # 50%, 25%, and 12.5% of "current" abund (model est for final data year from jags output * 3 for remig interval)
#
# pop.aboveT <- array(1, c(nsim, tmax, 3))  # order is: rows=sims, cols=yrs, matrices: 3 threshold matrices w/ nrows=10000 sims, ncols=100 yrs)
#
# for (th in 1:3) {
#   thr.pcnt <- thresh.vec[th]
#   for (i in 1:nsim) {
#     rs.CURF <- sum(exp(jags.model$sims.list$X[i,(fy-2):fy]))  # 3 year run sum of total females 2014-2016 (best estimate for current abundance for baseline)
#     ABUNDT <- thr.pcnt*rs.CURF   # abundance threshold 50% * 3 yr remig int * predicted Annual Nesters in final data year 2016
#     #ABUNDT <- thr.pcnt*3*exp(jags.model$sims.list$X[i,fy]) # abundance threshold 50% * 3 yr remig int * predicted Annual Nesters in final data year 2016
#
#     #for(tt in 3:tmax) if(sum(popsize[i,((tt-2):tt)]) <= ABUNDT) {pop.aboveT[i,tt:tmax,th] <- 0; break}
#     for(tt in 3:tmax) if(sum(popsize[i,((tt-2):tt)]) <= ABUNDT) {pop.aboveT[i,tt:tmax,th] <- 0}
#   }
# }
#
# # years to reach/fall below 50% current abundance threshold w/ 95% CI
# yrs.to.th50 <- apply(X=pop.aboveT[,,1], MARGIN=1, FUN=sum)  # for each row (sim run) add up the years it took to reach or fall below the abundance threshold
# mean(yrs.to.th50)
# quantile(yrs.to.th50, probs=c(0.025,0.50, 0.975))
# length(which(yrs.to.th50<100))  # number of years in which total females falls below threshold
#
# # years to reach/fall below 25% current abundance threshold  w/ 95% CI
# yrs.to.th25 <- apply(X=pop.aboveT[,,2], MARGIN=1, FUN=sum)  # for each row (sim run) add up the years it took to reach or fall below the abundance threshold
# mean(yrs.to.th25)
# quantile(yrs.to.th25, probs=c(0.025,0.50, 0.975))
#
# # years to reach/fall below 12.5% current abundance threshold  w/ 95% CI
# yrs.to.th12.5 <- apply(X=pop.aboveT[,,3], MARGIN=1, FUN=sum)  # for each row (sim run) add up the years it took to reach or fall below the abundance threshold
# mean(yrs.to.th12.5)
# quantile(yrs.to.th12.5, probs=c(0.025,0.50, 0.975))
#
#
#



##==============================================================================