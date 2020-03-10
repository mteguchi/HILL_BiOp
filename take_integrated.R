#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#	PIFSC Hawaiian SSLL Sea Turtle Take Model
#		
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Date Initialized: August 3, 2019
# Compiler: Zach Siders
# Authors: Tomo Eguchi, Summer Martin, Zach Siders, T. Todd Jones

###-----------------------------------------------------
#		Initialization
###-----------------------------------------------------
	#devtools::install_github( "James-Thorson/Conway-Maxwell-Poisson" )
	library(CMP)
	library(compoisson)
	library(mvtnorm)
	library(truncnorm)
	library(doParallel)  #Foreach Parallel Adaptor 
	library(foreach)     #Provides foreach looping construct
	library(abind)
	library(jagsUI)
	library(coda)
	library(tidyverse) #must be loaded after compoisson due to masking of MASS::select
	library(loo)
	library(rstan)
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())
	

	main.folder <- "/Volumes/HDD/Users/Zach/Documents/Projects/NOAA_BioOp_Take_LH_LB/Models/Integrated/"

	#####
	#####
	redo <- FALSE #FLAG TO RERUN ALL PREVIOUS OPS
	#####
	#####

	data.path <- paste0(main.folder,"Data/")
	imput.path <- paste0(main.folder,"Imputation/")
	trend.path <- paste0(main.folder,"Trend/")
	table.path <- paste0(main.folder,"Output/Tables/")
	fig.path <- paste0(main.folder,"Output/Figures/")
	sim.path <- paste0(main.folder,"Output/Simulations/")

	if(!dir.exists(data.path)) dir.create(data.path)
	if(!dir.exists(imput.path)) dir.create(imput.path)
	if(!dir.exists(trend.path)) dir.create(trend.path)
	if(!dir.exists(paste0(main.folder,"Output"))) dir.create(paste0(main.folder,"Output"))
	if(!dir.exists(table.path)) dir.create(table.path)
	if(!dir.exists(fig.path)) dir.create(fig.path)
	if(!dir.exists(sim.path)) dir.create(sim.path)

	source(paste0(main.folder,"take_helper_Fn.R"))

###-----------------------------------------------------
#		LOGGERHEAD GROWTH (led by Zach Siders)
###-----------------------------------------------------
	#Collaborators: Rob Ahrens, Nicholas Ducharme-Barth, T. Todd Jones, Summer Martin, Calandra N. Turner Tomaszewicz
	
	check.prerun <- file.exists(paste0(data.path,"VBGM_HMC_sims.Rdata"))
	if(!check.prerun | redo){
		CC.growth <- read.csv(paste0(data.path, 'NPCc CCL at est age.csv'))

		#Hatchling size (Nishimura 1967)
		CC.hatchling <- read.csv(paste0(data.path, "Nishimura_1967_Hatchling.csv"))

		hatchling_size <- mean(CC.hatchling$Length_cm)
		hatchling_sd <- sd(CC.hatchling$Length_cm)

		CC.nester <- read.csv(paste0(data.path,"Hatase_2002_Nesters.csv"))

		gen_nester <- unlist(apply(CC.nester, 1, function(x) {rnorm(x[3], x[1], x[2])}))
		nester_size <- mean(gen_nester)/10
		nester_sd <- sd(gen_nester)/10

		CC.growth.stan.dat <- list(n_obs = nrow(CC.growth),
		                 age = CC.growth$Est_Age+0.001,
		                 l = CC.growth$CCL,
		                 seq_ages = seq(0, max(CC.growth$Est_Age), by=0.1),
		                 nseq = length(seq(0, max(CC.growth$Est_Age), by=0.1)),
		                 nester_size = nester_size,
		                 nester_sd = nester_sd,
		                 n_hatchlings = nrow(CC.hatchling),
		                 hatchling_size = CC.hatchling$Length_cm)

		CC.growth.stan.fit.l0 <- stan(model_code=vbgm_stan_l0_ln,
		                 data = CC.growth.stan.dat,
		                 init = init.vbgm.l0,
		                 warmup = 5000,
		                 iter = 7500)
		
		CC.growth.stan.sims.l0 <- rstan::extract(CC.growth.stan.fit.l0)

		CC.VBGF <- with(CC.growth.stan.sims.l0, list(model="Lknot", Linf = median(Linf), K = median(K), Lknot = median(Lknot), Amat=median(Amat)))

		save(CC.growth.stan.sims.l0, file=paste0(data.path,"VBGM_HMC_sims.Rdata"))
	}else{
		load(file=paste0(data.path,"VBGM_HMC_sims.Rdata"))
		CC.VBGF <- with(CC.growth.stan.sims.l0, list(model="Lknot", Linf = median(Linf), K = median(K), Lknot = median(Lknot), Amat=median(Amat)))
	}
###-----------------------------------------------------
#		IMPUTATION (led by Tomo Eguchi)
###-----------------------------------------------------
	check.prerun <- file.exists(paste0(imput.path, "N_imput.Rdata"))
	if(!check.prerun | redo){
		###--------------------------------------------------
		#		Process Raw Observations of Nests (DC)
		###--------------------------------------------------

			period.JM <- 12
			period.W <- 6
			maxN <- 10000

			all.years <- 2001:2017
			idx <- 1:length(all.years)

			year.begin.JM <- 2001
			year.end <- 2017
			data.jags.JM <- data.extract(location = "JM", 
			                             year.begin = year.begin.JM, 
			                             year.end = year.end,
			                             file.path = data.path)

			JM.keep <- 2001:2017
			idx.JM <- idx[all.years %in% JM.keep]
			n.keep.JM <- length(idx.JM)
			dt.JM <- idx.JM[2:length(idx.JM)] - idx.JM[1:(length(idx.JM)-1)]

			year.begin.W <- 2006
			data.jags.W <- data.extract(location = "W", 
			                            year.begin = year.begin.W, 
			                            year.end = year.end,
			                             file.path = data.path)
			W.keep <- c(2006, 2007, 2008, 2009, 2010, 2011, 2012, 2016, 2017)
			idx.W <- idx[all.years %in% W.keep]
			n.keep.W <- length(idx.W)
			dt.W <- idx.W[2:length(idx.W)] - idx.W[1:(length(idx.W)-1)]

			#Combine datasets for analysis
			# JM has more data than W, so we need to pad W data
			y.W <- rbind(array(data = NA, 
			                   dim = c(nrow(data.jags.JM$jags.data2$y) - nrow(data.jags.W$jags.data2$y),
			                           ncol(data.jags.JM$jags.data2$y))),
			             data.jags.W$jags.data2$y)

			y <- cbind(as.vector(t(data.jags.JM$jags.data2$y)), as.vector(t(y.W)))
			y.raw <- as.data.frame(apply(y,2,exp))
			colnames(y.raw) <- c("JM","W")
			y.raw$year.frac <- rep(2001:2017, each=12) + seq(1,12)/12
			save(y.raw, file=paste0(data.path,"DC_raw_nests_month.Rdata"))

			years <- rep(2001:2017, each = 12)
			n.years <- 17

			# for estimating U  ######
			n.timeseries <- ncol(y)
		###--------------------------------------------------
		#		JAGS Implementation
		###--------------------------------------------------

			jags.data <- list(y = y,
	                  m = rep(1:12, times = n.years),
	                  n.steps = nrow(y),
	                  n.months = 12,
	                  pi = pi,
	                  period = c(period.JM, period.W),
	                  n.timeseries = n.timeseries,
	                  n.years = n.years)

			jags.params <- c("c", "beta.cos", "beta.sin",
	                 'sigma.X', "sigma.y", "N", 
	                  "y", "X", "deviance")

			MCMC.params <- list(n.chains = 5,
                    n.samples = 100000,
                    n.burnin = 50000,
                    n.thin = 5)

			jm <- jags(jags.data,
		           inits = NULL,
		           parameters.to.save= jags.params,
		           model.file = paste0(main.folder,'model_norm_norm_Four_imputation.txt'),
		           n.chains = MCMC.params$n.chains,
		           n.burnin = MCMC.params$n.burnin,
		           n.thin = MCMC.params$n.thin,
		           n.iter = MCMC.params$n.samples,
		           DIC = T, parallel=T)

			save(jm, file=paste0(imput.path,"JAGS_model_run.Rdata"))
		###--------------------------------------------------
		#		Summary of JAGS model
		###--------------------------------------------------

			Ns.stats.JM <- data.frame(time = year.begin.JM:year.end,
		                          low = as.vector(t(jm$q2.5$N[,1])),
		                          median = as.vector(t(jm$q50$N[,1])),
		                          high = as.vector(t(jm$q97.5$N[,1])))
			Ns.stats.JM$location <- "Jamursba-Medi"

			Ns.stats.W <- data.frame(time = year.begin.JM:year.end,
			                         low = as.vector(t(jm$q2.5$N[,2])),
			                         median = as.vector(t(jm$q50$N[,2])),
			                         high = as.vector(t(jm$q97.5$N[,2])))
			Ns.stats.W$location <- "Wermon"

			Ns.stats <- rbind(Ns.stats.JM, Ns.stats.W)

			save(Ns.stats.JM, Ns.stats.W, file=paste0(imput.path, "N_imput.Rdata"))
	}else{
		load(paste0(imput.path, "N_imput.Rdata"))
	}
###-----------------------------------------------------
#		Historical ANE (led by Zach Siders)
###-----------------------------------------------------
	check.prerun <- file.exists(paste0(data.path,"historical_ANE.csv"))
	if(!check.prerun | redo){
		###--------------------------------------------------
		#		Take Demographics Data
		###--------------------------------------------------
			#-------
			# Historical Take	
			
			td.dat <- read.csv(paste0(data.path,"turtle_l_m.csv"), stringsAsFactors=FALSE)

			colnames(td.dat) <- c("Year", "Spp", "M.low", "M.high", "Len")

			td.dat <- td.dat[td.dat$Spp!="" & td.dat$Year != 2019,]
			td.dat <- td.dat[td.dat$Spp %in% c("Leatherback", "Loggerhead"),]

			DC.VBGF <- list(Linf = 142.7, K = 0.2262, tknot=-0.17, Amat = 16.1)
			CC.VBGF <- list(Linf = 80.4473850, K = 0.1396317, Lknot=4.7363329, Amat = 26.4950786)

			DC.Pj <- 0.81
			DC.Pa <- 0.893
			CC.Pj <- 0.8
			CC.Pa <- 0.895

			DC.PF <- 0.73
			CC.PF <- 0.65
		###--------------------------------------------------
		#		Individual Characteristics
		###--------------------------------------------------
		
			pred.Age <- seq(0,100, by=0.01) #sequence of ages to back calculate over
			DC.Lpred <- with(DC.VBGF, Linf * (1-exp(-K*(pred.Age-tknot))))
			CC.Lpred <- with(CC.VBGF, Linf - (Linf-Lknot)*exp(-K*pred.Age))

			td.dat$Len[is.na(td.dat$Len) & td.dat$Spp == "Leatherback"] <- median(td.dat$Len[td.dat$Spp == "Leatherback"], na.rm=T)
			td.dat$Len[is.na(td.dat$Len) & td.dat$Spp == "Loggerhead"] <- median(td.dat$Len[td.dat$Spp == "Loggerhead"], na.rm=T)
			
			td.dat$M.mu <- rowMeans(td.dat[,c("M.low","M.high")], na.rm=T)
			td.dat$M.mu[is.nan(td.dat$M.mu)] <- mean(td.dat$M.mu[!is.nan(td.dat$M.mu) & td.dat$Spp == td.dat$Spp[is.nan(td.dat$M.mu)]])

			td.dat$Age <- NA
			td.dat$Age[td.dat$Spp == "Leatherback"] <- sapply(td.dat$Len[td.dat$Spp == "Leatherback"], function(x) {pred.Age[which.min(abs(x - DC.Lpred))]})
			td.dat$Age[td.dat$Spp == "Loggerhead"] <- sapply(td.dat$Len[td.dat$Spp == "Loggerhead"], function(x) {pred.Age[which.min(abs(x - CC.Lpred))]})

			#determine the animals stage based on maturity
			td.dat$Stage <- NA
			td.dat$Stage[td.dat$Spp == "Leatherback"] <- ifelse(td.dat$Age[td.dat$Spp == "Leatherback"] > DC.VBGF$Amat, "A", "J")
			td.dat$Stage[td.dat$Spp == "Loggerhead"] <- ifelse(td.dat$Age[td.dat$Spp == "Loggerhead"] > CC.VBGF$Amat, "A", "J")

			#calculate the number of years until maturity
			td.dat$YatLarge <- NA
			td.dat$YatLarge[td.dat$Spp == "Leatherback"] <- DC.VBGF$Amat - td.dat$Age[td.dat$Spp == "Leatherback"]
			td.dat$YatLarge[td.dat$Spp == "Loggerhead"] <- CC.VBGF$Amat - td.dat$Age[td.dat$Spp == "Loggerhead"]
		###--------------------------------------------------
		#		back project for the Remigration Interval
		###--------------------------------------------------
		
			#back project the RI
			td.dat$YatLargeR <- round(td.dat$YatLarge)
			td.dat$YatLargeR[td.dat$YatLargeR < 0] <- 0

			td.dat$Nest1 <- td.dat$Year + td.dat$YatLargeR
			td.dat$Nest2 <- td.dat$Nest1 + 3
			td.dat$Nest3 <- td.dat$Nest2 + 3
			td.dat$Nest4 <- td.dat$Nest3 + 3
			td.dat$Nest5 <- td.dat$Nest4 + 3
			td.dat$Nest6 <- td.dat$Nest5 + 3
			td.dat$Nest7 <- td.dat$Nest6 + 3
		###--------------------------------------------------
		#		Convert to ANEs
		###--------------------------------------------------
		
			ANEyr <- td.dat[,paste0("Nest",1:7)] - td.dat$Year
			ANEproj <- ANEyr
			ANEproj$Nest1 <- (ifelse(td.dat$Spp == "Leatherback", DC.Pj, CC.Pj) ^ ANEproj$Nest1) * ifelse(td.dat$Spp == "Leatherback", DC.PF, CC.PF) * td.dat$M.mu
			ANEproj[,-1] <- ifelse(td.dat$Spp == "Leatherback", DC.Pa, CC.Pa) ^ ANEproj[,-1] * ANEproj$Nest1

			ANEproj.nFD <- ANEyr
			ANEproj.nFD$Nest1 <- (ifelse(td.dat$Spp == "Leatherback", DC.Pj, CC.Pj) ^ ANEproj.nFD$Nest1)
			ANEproj.nFD[,-1] <- ifelse(td.dat$Spp == "Leatherback", DC.Pa, CC.Pa) ^ ANEproj.nFD[,-1] * ANEproj.nFD$Nest1

			ANEp <- data.frame(nestyr = unlist(td.dat[,paste0("Nest",1:7)]),
			                   ANE = unlist(ANEproj), 
			                   Spp = rep(td.dat$Spp, 7))

			ANEtab <- aggregate(ANE ~ Spp + nestyr, data=ANEp, FUN=sum)

			write.csv(ANEtab, file=paste0(data.path,"historical_ANE.csv"), row.names=FALSE)
	}else{
		ANEtab <- read.csv(paste0(data.path,"historical_ANE.csv"))
	}


	check.prerun <- file.exists(paste0(data.path, "take_demo_mvn.Rdata"))

	if(!check.prerun | redo){
		###--------------------------------------------------
		#		Fit multivariate take demographics
		###--------------------------------------------------
			tab.spp <- as.data.frame(with(td.dat, table(Year, Spp)))

			DC.td.df <- na.omit(td.dat[td.dat$Spp=="Leatherback",])

			DC.td.tab <- tab.spp[tab.spp$Spp=="Leatherback",c("Year","Spp","Freq")]

			DC.td.dat <- list(N = nrow(DC.td.df),
			               x = cbind(log(DC.td.df$L), boot::logit(DC.td.df$M.mu)),
			               nyear = nrow(DC.td.tab),
			               year = as.integer(factor(DC.td.df$Year)),
			               rtl = DC.td.tab$Freq)
			DC.td.fit <- stan(model_code = mvnorm,
			               data = DC.td.dat,
			               init = mv.DC.init,
			               warmup = 5000,
			               iter = 7500)
			DC.td.sims <- rstan::extract(DC.td.fit)

			###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			# LOGGERHEADS

			CC.td.df <- na.omit(td.dat[td.dat$Spp=="Loggerhead",])
			CC.td.df$M.mu[CC.td.df$M.mu==1] <- 0.999

			CC.td.tab <- tab.spp[tab.spp$Spp=="Loggerhead",c("Year","Spp","Freq")]

			CC.td.dat <- list(N = nrow(CC.td.df),
			               x = cbind(log(CC.td.df$L), boot::logit(CC.td.df$M.mu)),
			               nyear = nrow(CC.td.tab),
			               year = as.integer(factor(CC.td.df$Year)),
			               rtl = CC.td.tab$Freq)		

			CC.td.fit <- stan(model_code = mvnorm,
			               data = CC.td.dat,
			               init = mv.CC.init,
			               warmup = 5000,
			               iter = 7500)
			CC.td.sims <- rstan::extract(CC.td.fit)
		###--------------------------------------------------
		#		Summarize take Demographics
		###--------------------------------------------------
			DC.mu0 <- median(DC.td.sims$mu0)
			DC.beta0 <- median(DC.td.sims$beta0)
			DC.beta1 <- median(DC.td.sims$beta1)
			DC.sigma <- apply(DC.td.sims$sigma, 2, median)
			DC.rho <- median(DC.td.sims$rho)
			DC.cov <- matrix(c(DC.sigma[1]^2, DC.sigma[1]*DC.sigma[2]*DC.rho, DC.sigma[1]*DC.sigma[2]*DC.rho, DC.sigma[2]^2),2,2)

			CC.mu0 <- median(CC.td.sims$mu0)
			CC.beta0 <- median(CC.td.sims$beta0)
			CC.beta1 <- median(CC.td.sims$beta1)
			CC.sigma <- apply(CC.td.sims$sigma, 2, median)
			CC.rho <- median(CC.td.sims$rho)
			CC.cov <- matrix(c(CC.sigma[1]^2, CC.sigma[1]*CC.sigma[2]*CC.rho, CC.sigma[1]*CC.sigma[2]*CC.rho, CC.sigma[2]^2),2,2)

			
			save(DC.mu0, DC.beta0, DC.beta1, DC.cov, CC.mu0, CC.beta0, CC.beta1, CC.cov, file=paste0(data.path,"take_demo_mvn.Rdata"))
	}else{
		load(paste0(data.path,"take_demo_mvn.Rdata"))
	}
###-----------------------------------------------------
#		Trend Analysis (led by Summer Martin)
###-----------------------------------------------------
	
	check.prerun <- all(file.exists(c(paste0(trend.path,"Historical.ANE.Scenarios/Dc_JM&W_MEDIAN_sUQ"),
	                              paste0(trend.path,"Historical.ANE.Scenarios/Dc_JM&W_LOW_sUQ"),
	                              paste0(trend.path,"Historical.ANE.Scenarios/Dc_JM&W_HIGH_sUQ"),
	                              paste0(trend.path,"Historical.ANE.Scenarios/Cc_Yakushima_sUQ"))))
	if(!check.prerun | redo){ 
		###--------------------------------------------------
		#		 Data Read In
		###--------------------------------------------------
			#-------
			# Yakushima 3 beaches (Loggerheads)

			CC.dat <- read.csv(paste0(data.path,"Yakushima_data_for_BiOp.csv"), header=T)
			
			clutch.freq <- 3             # clutch frequency = nests per female in a season (for converting nests to females)
			CC.dat[ ,c("Females_Mae", "Females_Inak", "Females_Yotsu")] <- CC.dat[ , 2:4]/clutch.freq   # use CF to convert from nests to females that nest each year

			#-------
			# Jamursba-Medi & Wermon (Leatherbacks)

			#-------
			# Jamursba-Medi
			
			JM.dat <- Ns.stats.JM   # only keep specific data years per group discussions
			JM.dat[,2:4] <- exp(JM.dat[,2:4])
			clutch.freq <- 5.5    # clutch frequency =nests per female in a season; Tapilatu etal 2013 = mean CF 5.5 +/- 1.6
			JM.dat$Females_median <- JM.dat$median/clutch.freq   # use CF to convert from nests to females that nest each year
			JM.dat$Females_low <- JM.dat$low/clutch.freq   # use CF to convert from nests to females that nest each year
			JM.dat$Females_high <- JM.dat$high/clutch.freq   # use CF to convert from nests to females that nest each year

			#-------
			# Wermon
			analysis.yrs <- c(2006:2017)
			W.dat <- Ns.stats.W   # only keep specific data years per group discussions
			W.dat[,2:4] <- exp(W.dat[,2:4])
			W.dat <- W.dat[W.dat$time %in% analysis.yrs,]
			W.dat$Females_median <- W.dat$median/clutch.freq   # use CF to convert from nests to females that nest each year
			W.dat$Females_low <- W.dat$low/clutch.freq   # use CF to convert from nests to females that nest each year
			W.dat$Females_high <- W.dat$high/clutch.freq   # use CF to convert from nests to females that nest each year

			colnames(JM.dat)[1] <- colnames(W.dat)[1] <- "Season"
		
			#-------
			# Combine

			DC.dat.med <- merge(x=JM.dat[,c("Season", "Females_median")], y=W.dat[,c("Season", "Females_median")], by="Season", all.x=TRUE, all.y=TRUE)

			DC.dat.low <- merge(x=JM.dat[,c("Season", "Females_low")], y=W.dat[,c("Season", "Females_low")], by="Season", all.x=TRUE, all.y=TRUE)

			DC.dat.high <- merge(x=JM.dat[,c("Season", "Females_high")], y=W.dat[,c("Season", "Females_high")], by="Season", all.x=TRUE, all.y=TRUE)

			DC.raw.dat <- merge(x=JM.dat[,c("Season","low","median","high")], y = W.dat[,c("Season","low","median","high")], by="Season", all.x=TRUE, all.y=TRUE, suffixes=c(".JM",".W"))


			CC.dat.iT <- CC.dat
			DC.dat.med.iT <- DC.dat.med
			DC.dat.low.iT <- DC.dat.low
			DC.dat.high.iT <- DC.dat.high
		###--------------------------------------------------
		#		Adding the historical ANE
		###--------------------------------------------------

			#-------
			# ANE Loggerheads	

				CC.ANE <- ANEtab[ANEtab$Spp == "Loggerhead",]
				#proportion by beach
				prop.CC1 <- CC.dat[,5]/rowSums(CC.dat[5:7], na.rm=T)
				prop.CC2 <- CC.dat[,6]/rowSums(CC.dat[5:7], na.rm=T)
				prop.CC3 <- CC.dat[,7]/rowSums(CC.dat[5:7], na.rm=T)

				#tack it back in 
				CC.dat[CC.dat$Year %in% CC.ANE$nestyr, 5] <- CC.dat[CC.dat$Year %in% CC.ANE$nestyr, 5] + CC.ANE[CC.ANE$nestyr %in% CC.dat$Year,"ANE"] * prop.CC1[CC.dat$Year %in% CC.ANE$nestyr]
				CC.dat[CC.dat$Year %in% CC.ANE$nestyr, 6] <- CC.dat[CC.dat$Year %in% CC.ANE$nestyr, 6] + CC.ANE[CC.ANE$nestyr %in% CC.dat$Year,"ANE"] * prop.CC2[CC.dat$Year %in% CC.ANE$nestyr]
				CC.dat[CC.dat$Year %in% CC.ANE$nestyr, 7] <- CC.dat[CC.dat$Year %in% CC.ANE$nestyr, 7] + CC.ANE[CC.ANE$nestyr %in% CC.dat$Year,"ANE"] * prop.CC3[CC.dat$Year %in% CC.ANE$nestyr]

			#-------
			# ANE Leatherbacks
				DC.ANE <- ANEtab[ANEtab$Spp == "Leatherback",]

				prop.DC <- DC.dat.med[,2]/rowSums(DC.dat.med[,2:3], na.rm=T)

				DC.dat.med[DC.dat.med$Season %in% DC.ANE$nestyr, 2] <- DC.dat.med[DC.dat.med$Season %in% DC.ANE$nestyr, 2] + DC.ANE[DC.ANE$nestyr %in% DC.dat.med$Season,"ANE"] * prop.DC[DC.dat.med$Season %in% DC.ANE$nestyr]
				DC.dat.med[DC.dat.med$Season %in% DC.ANE$nestyr, 3] <- DC.dat.med[DC.dat.med$Season %in% DC.ANE$nestyr, 3] + DC.ANE[DC.ANE$nestyr %in% DC.dat.med$Season,"ANE"] * (1-prop.DC[DC.dat.med$Season %in% DC.ANE$nestyr])

				DC.dat.low[DC.dat.low$Season %in% DC.ANE$nestyr, 2] <- DC.dat.low[DC.dat.low$Season %in% DC.ANE$nestyr, 2] + DC.ANE[DC.ANE$nestyr %in% DC.dat.low$Season,"ANE"] * prop.DC[DC.dat.low$Season %in% DC.ANE$nestyr]
				DC.dat.low[DC.dat.low$Season %in% DC.ANE$nestyr, 3] <- DC.dat.low[DC.dat.low$Season %in% DC.ANE$nestyr, 3] + DC.ANE[DC.ANE$nestyr %in% DC.dat.low$Season,"ANE"] * (1-prop.DC[DC.dat.low$Season %in% DC.ANE$nestyr])

				DC.dat.high[DC.dat.high$Season %in% DC.ANE$nestyr, 2] <- DC.dat.high[DC.dat.high$Season %in% DC.ANE$nestyr, 2] + DC.ANE[DC.ANE$nestyr %in% DC.dat.high$Season,"ANE"] * prop.DC[DC.dat.high$Season %in% DC.ANE$nestyr]
				DC.dat.high[DC.dat.high$Season %in% DC.ANE$nestyr, 3] <- DC.dat.high[DC.dat.high$Season %in% DC.ANE$nestyr, 3] + DC.ANE[DC.ANE$nestyr %in% DC.dat.high$Season,"ANE"] * (1-prop.DC[DC.dat.high$Season %in% DC.ANE$nestyr])
		###--------------------------------------------------
		#		Running the Trend analysis
		###--------------------------------------------------
			
			###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			### Historical ANE added in

			hist.scen.path <- paste0(trend.path,"Historical.ANE.Scenarios/")
			if(!dir.exists(hist.scen.path)) dir.create(hist.scen.path)
			# ----------------------------------------------
			# 1 - Loggerheads - using 3 beaches from Yakushima, Japan 
			# ----------------------------------------------
			scenario <- "Cc_Yakushima_sUQ"
			thedata.loggers <- CC.dat
			remigLH <- 3.3
			clutch.freqLH <- 4.6 
			CC.file.tag <- paste(scenario,"_", Sys.Date(), sep="")
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(hist.scen.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))
			# ----------------------------------------------

			# ----------------------------------------------
			# 2A - Leatherbacks. Using MEDIAN annual nest count estimates Tomo's imputation model output.
			# ----------------------------------------------
			scenario <- "Dc_JM&W_MEDIAN_sUQ"
			thedata.leathers <- DC.dat.med
			remigLB <- 3.06
			clutch.freqLB <- 5.5
			DC.file.tag <- paste(scenario,"_", Sys.Date(), sep="") 
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(hist.scen.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))
			# ----------------------------------------------

			# ----------------------------------------------
			# 2B - Leatherbacks. Using LOW (lower 95% CI value) annual nest count estimates from Tomo's imputation model output.
			# ----------------------------------------------
			scenario <- "Dc_JM&W_LOW_sUQ"
			thedata.leathers <- DC.dat.low
			remigLB <- 3.06
			clutch.freqLB <- 5.5 
			DC.file.tag.low <- paste(scenario,"_", Sys.Date(), sep="")
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(hist.scen.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))
			# ----------------------------------------------

			# ----------------------------------------------
			# 2C - Leatherbacks. Using HIGH (upper 95% CI value) annual nest count estimates from Tomo's imputation model output.
			# ----------------------------------------------
			scenario <- "Dc_JM&W_HIGH_sUQ"
			thedata.leathers <- DC.dat.high
			remigLB <- 3
			clutch.freqLB <- 5.5 
			DC.file.tag.high <- paste(scenario,"_", Sys.Date(), sep="")
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(hist.scen.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))

			###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			### No Historical ANE added in

			include.hist.take.path <- paste0(trend.path,"Include.Hist.Take.Scenarios/")
			if(!dir.exists(include.hist.take.path)) dir.create(include.hist.take.path)
			# ----------------------------------------------
			# 1 - Loggerheads - using 3 beaches from Yakushima, Japan 
			# ----------------------------------------------
			scenario <- "Cc_Yakushima_sUQ"
			thedata.loggers <- CC.dat.iT
			remigLH <- 3.3
			clutch.freqLH <- 4.6 
			CC.file.tag <- paste(scenario,"_", Sys.Date(), sep="")
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(include.hist.take.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))
			# ----------------------------------------------

			# ----------------------------------------------
			# 2A - Leatherbacks. Using MEDIAN annual nest count estimates Tomo's imputation model output.
			# ----------------------------------------------
			scenario <- "Dc_JM&W_MEDIAN_sUQ"
			thedata.leathers <- DC.dat.med.iT
			remigLB <- 3.06
			clutch.freqLB <- 5.5
			DC.file.tag <- paste(scenario,"_", Sys.Date(), sep="") 
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(include.hist.take.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))
			# ----------------------------------------------

			# ----------------------------------------------
			# 2B - Leatherbacks. Using LOW (lower 95% CI value) annual nest count estimates from Tomo's imputation model output.
			# ----------------------------------------------
			scenario <- "Dc_JM&W_LOW_sUQ"
			thedata.leathers <- DC.dat.low.iT
			remigLB <- 3.06
			clutch.freqLB <- 5.5 
			DC.file.tag.low <- paste(scenario,"_", Sys.Date(), sep="")
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(include.hist.take.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))
			# ----------------------------------------------

			# ----------------------------------------------
			# 2C - Leatherbacks. Using HIGH (upper 95% CI value) annual nest count estimates from Tomo's imputation model output.
			# ----------------------------------------------
			scenario <- "Dc_JM&W_HIGH_sUQ"
			thedata.leathers <- DC.dat.high.iT
			remigLB <- 3
			clutch.freqLB <- 5.5 
			DC.file.tag.high <- paste(scenario,"_", Sys.Date(), sep="")
			rsel <- "rdyn"   # future predictions: choose "rdyn" (dynamic) or "rstat" (static) for future predictions
			save.dir <- paste0(include.hist.take.path, scenario)
			source(paste0(trend.path,"singleUQ_indeptUQs_PROJECTIONS.R"))
	}else{
			hist.scen.path <- paste0(trend.path,"Historical.ANE.Scenarios/")
			scenario <- "Cc_Yakushima_sUQ"
			file.date <- gsub("_0_data_used.csv","",gsub(paste0(scenario,"_"),"",grep("_data_used.csv",list.files(paste0(hist.scen.path,scenario)), value=T)))
			file.date <- file.date[length(file.date)]
			CC.file.tag <- paste(scenario,"_", file.date, sep="")
			scenario <- "Dc_JM&W_MEDIAN_sUQ"
			DC.file.tag <- paste(scenario,"_", file.date, sep="") 
			scenario <- "Dc_JM&W_LOW_sUQ"
			DC.file.tag.low <- paste(scenario,"_", file.date, sep="")
			scenario <- "Dc_JM&W_HIGH_sUQ"
			DC.file.tag.high <- paste(scenario,"_", file.date, sep="")
			include.hist.take.path <- paste0(trend.path,"Include.Hist.Take.Scenarios/")
	}
	
###-----------------------------------------------------
#		Take Model (led by Zach Siders)
###-----------------------------------------------------
	###-----------------------------------------------------
	#		Initialization
	###-----------------------------------------------------
		#store DC MVN Take Demographic Parameters
		TD_DC_MVN = list(mu0 = DC.mu0, beta0 = DC.beta0, beta1 = DC.beta1, cov = DC.cov)
		#store CC MVN Take Demographic Parameters
		TD_CC_MVN = list(mu0 = CC.mu0, beta0 = CC.beta0, beta1 = CC.beta1, cov = CC.cov)

		#Conway-Maxwell Poisson parameters for 3 year segements
		#Anticipated Take Level from Marti McCracken's model for Dc
		DC_ATL <- list(mu1 = 2.124568, nu1 = 0.4805365, mu2 = 2.344938, nu2 = 0.141262, mu3 = 0.03930914, nu3 = 0.1149811)
		#Anticipated Take Level from Marti McCracken's model for Cc
		CC_ATL = list(mu1 = 3.444104, nu1 = 0.06451325, mu2 = 7.506372e-05, nu2 =  0.01686626, mu3 = NULL, nu3 = NULL)

		#Lontoh (2014) Dissertation on Dc; # of animals with a given remigration interval
		DC.RI <- c(1, 13, 7, 6, 4, 1)
		#Fit a Conway-Maxwell Poisson to the data
		DC.CMP <- com.fit(matrix(c(1:6, DC.RI), ncol=2))
	###-----------------------------------------------------
	#		Loading the simulations
	###-----------------------------------------------------
		#read in the U, Q, and final years N
		DC.trend.sims <- read.csv(paste0(hist.scen.path,"Dc_JM&W_MEDIAN_sUQ/Posteriors_U_and_AbundFinalYrs_",DC.file.tag,".csv"), header=T)
		DC.trend.sims.low <- read.csv(paste0(hist.scen.path,"Dc_JM&W_LOW_sUQ/Posteriors_U_and_AbundFinalYrs_",DC.file.tag.low,".csv"), header=T)
		DC.trend.sims.high <- read.csv(paste0(hist.scen.path,"Dc_JM&W_HIGH_sUQ/Posteriors_U_and_AbundFinalYrs_",DC.file.tag.high,".csv"), header=T)
		#includes historical take in U estimate
		DC.trend.sims.iT <- read.csv(paste0(include.hist.take.path,"Dc_JM&W_MEDIAN_sUQ/Posteriors_U_and_AbundFinalYrs_",DC.file.tag,".csv"), header=T)
		DC.trend.sims.iT.low <- read.csv(paste0(include.hist.take.path,"Dc_JM&W_LOW_sUQ/Posteriors_U_and_AbundFinalYrs_",DC.file.tag.low,".csv"), header=T)
		DC.trend.sims.iT.high <- read.csv(paste0(include.hist.take.path,"Dc_JM&W_HIGH_sUQ/Posteriors_U_and_AbundFinalYrs_",DC.file.tag.high,".csv"), header=T)

		CC.trend.sims <- read.csv(paste0(hist.scen.path,"Cc_Yakushima_sUQ/Posteriors_U_and_AbundFinalYrs_",CC.file.tag,".csv"), header=T)
		CC.trend.sims.iT <- read.csv(paste0(include.hist.take.path,"Cc_Yakushima_sUQ/Posteriors_U_and_AbundFinalYrs_",CC.file.tag,".csv"), header=T)

	###-----------------------------------------------------
	#		Simulation Summary
	###-----------------------------------------------------
		nsims <- nrow(DC.trend.sims)
		dc.trend <- DC.trend.sims[,c("U","Q","N_fym0")]
		colnames(dc.trend)[3] <- "N0"
		dc.trend.low <- DC.trend.sims.low[,c("U","Q","N_fym0")]
		colnames(dc.trend.low)[3] <- "N0"
		dc.trend.high <- DC.trend.sims.high[,c("U","Q","N_fym0")]
		colnames(dc.trend.high)[3] <- "N0"

		cc.trend <- CC.trend.sims[,c("U","Q","N_fym0")]
		colnames(cc.trend)[3] <- "N0"

		#-------
		# Write out Current Abundance	
		
			write.csv(curr.abund.fn(DC.trend.sims, 3.06), file=paste0(table.path,"DC.Curr.Abund.csv"))
			write.csv(curr.abund.fn(DC.trend.sims.low, 3.06), file=paste0(table.path,"DC.Curr.Abund.LOW.csv"))
			write.csv(curr.abund.fn(DC.trend.sims.high, 3.06), file=paste0(table.path,"DC.Curr.Abund.HIGH.csv"))

			write.csv(curr.abund.fn(CC.trend.sims, 3.3), file=paste0(table.path,"CC.Curr.Abund.csv"))

		#-------
		# Write out U summary
			dc.u.sum <- data.frame(notake = sum.fn(dc.trend$U),
			                    take = sum.fn(DC.trend.sims.iT$U))
			dc.u.sum.low <- data.frame(notake = sum.fn(dc.trend.low$U),
			                    take = sum.fn(DC.trend.sims.iT.low$U))
			dc.u.sum.high <- data.frame(notake = sum.fn(dc.trend.high$U),
			                    take = sum.fn(DC.trend.sims.iT.high$U))

			cc.u.sum <- data.frame(notake = sum.fn(cc.trend$U),
			                    take = sum.fn(CC.trend.sims.iT$U))

			write.csv(round(dc.u.sum,3), file=paste0(table.path,"DC.U_summary.csv"))
			write.csv(round(dc.u.sum.low,3), file=paste0(table.path,"DC.U_summary.LOW.csv"))
			write.csv(round(dc.u.sum.high,3), file=paste0(table.path,"DC.U_summary.HIGH.csv"))
			write.csv(round(cc.u.sum,3), file=paste0(table.path,"CC.U_summary.csv"))

		#-------
		# Plot overlapping U's
			
			pdf(file=paste0(fig.path, "U_Leatherbacks.pdf"), width=4, height=4)
			par(mar=c(4,4,1,1))
				u.den(dc.trend$U, DC.trend.sims.iT$U)
				fig.lab("W. Pac.\nLeatherbacks", xscale=0.8, yscale=0.9, cex=1)
				abline(v=0, col='grey80')
			dev.off()

			pdf(file=paste0(fig.path, "U_Leatherbacks.LOW.pdf"), width=4, height=4)
			par(mar=c(4,4,1,1))
				u.den(dc.trend.low$U, DC.trend.sims.iT.low$U)
				fig.lab("W. Pac.\nLeatherbacks", xscale=0.8, yscale=0.9, cex=1)
				abline(v=0, col='grey80')
			dev.off()

			pdf(file=paste0(fig.path, "U_Leatherbacks.HIGH.pdf"), width=4, height=4)
			par(mar=c(4,4,1,1))
				u.den(dc.trend.high$U, DC.trend.sims.iT.high$U)
				fig.lab("W. Pac.\nLeatherbacks", xscale=0.8, yscale=0.9, cex=1)
				abline(v=0, col='grey80')
			dev.off()

			pdf(file=paste0(fig.path, "U_Loggerheads.pdf"), width=4, height=4)
			par(mar=c(4,4,1,1))
				u.den(cc.trend$U, CC.trend.sims.iT$U)
				fig.lab("N. Pac.\nLoggerheads", xscale=0.8, yscale=0.9, cex=1)
				abline(v=0, col='grey80')
			dev.off()

		#-------
		# Plot Joint
			plot.joint(dc.trend, "DC", 6)
			plot.joint(dc.trend.low, "DC.LOW", 6)
			plot.joint(dc.trend.high, "DC.HIGH", 6)
			plot.joint(cc.trend, "CC", 6)
	###-----------------------------------------------------
	#		Define Scenarios (Parameter Read In)
	###-----------------------------------------------------
		DC.scenario <- list(
	                 n_y = 100, #number of years in scenario
	                 VBGF = list(model="tknot",Linf = 142.7, K = 0.2262, tknot=-0.17, Amat = 16.1), #VBGF parameters
	                 TD_MVN = TD_DC_MVN, #take demographics MVN parameters
	                 PjV = list(mu = 0.81, sd = 0.03), #juvenile survival
	                 PaV = list(mu = 0.893, sd = 0.013), #adult survival
	                 RI.dist = c("CMP"), #distribution of RI
	                 RV = list(RI = com.mean(DC.CMP$lambda, DC.CMP$nu), RI_sd = DC.CMP$nu, CF = 5.5, CF_sd = 1.6, CS = 77.9, CS_sd = 2.35, PF = 0.73, PF_sd = 0), #reproductive values
	                 ATLp = list(mu1 = 2.124568, nu1 = 0.4805365, mu2 = 2.344938, nu2 = 0.141262, mu3 = 0.03930914, nu3 = 0.1149811),
	                 ANE_dyn = FALSE, #whether ANE is dynamic
	                 Surv_dyn = FALSE,
	                 dynUQ = TRUE,
	                 ATL_scale = FALSE,
	                 grim.reaper=FALSE)

		CC.scenario <- list(
	                 n_y = 100, #number of years in scenario
	                 VBGF = list(model="Lknot",Linf = 80.4473850, K = 0.1396317, Lknot=4.7363329, Amat = 26.4950786), #VBGF parameters
	                 TD_MVN = TD_CC_MVN, #take demographics MVN parameters
	                 PjV = list(mu = 0.8, sd = 0.031), #juvenile survival
	                 PaV = list(mu = 0.895, sd = 0.028), #adult survival
	                 RI.dist = "normal",
	                 RV = list(RI = 3.3, RI_sd = 2.3, CF = 4.6, CF_sd = 1.1, CS = 122, CS_sd = 18.4, PF = 0.65, PF_sd = 0), #reproductive values
	                 ATLp = list(mu1 = 3.444104, nu1 = 0.06451325, mu2 = 7.506372e-05, nu2 =  0.01686626, mu3 = NULL, nu3 = NULL),
	                 ANE_dyn = FALSE, #whether ANE is dynamic
	                 Surv_dyn = FALSE,
	                 dynUQ = FALSE,
	                 ATL_scale = FALSE,
	                 grim.reaper=FALSE)
		DC.scenario.GR <- list(
	                 n_y = 100, #number of years in scenario
	                 VBGF = list(model="tknot",Linf = 142.7, K = 0.2262, tknot=-0.17, Amat = 16.1), #VBGF parameters
	                 TD_MVN = TD_DC_MVN, #take demographics MVN parameters
	                 PjV = list(mu = 0.81, sd = 0.03), #juvenile survival
	                 PaV = list(mu = 0.893, sd = 0.013), #adult survival
	                 RI.dist = c("CMP"), #distribution of RI
	                 RV = list(RI = com.mean(DC.CMP$lambda, DC.CMP$nu), RI_sd = DC.CMP$nu, CF = 5.5, CF_sd = 1.6, CS = 77.9, CS_sd = 2.35, PF = 0.73, PF_sd = 0), #reproductive values
	                 ATLp = list(mu1 = 2.124568, nu1 = 0.4805365, mu2 = 2.344938, nu2 = 0.141262, mu3 = 0.03930914, nu3 = 0.1149811),
	                 ANE_dyn = FALSE, #whether ANE is dynamic
	                 Surv_dyn = FALSE,
	                 dynUQ = TRUE,
	                 ATL_scale = FALSE,
	                 grim.reaper=TRUE)

		CC.scenario.GR <- list(
	                 n_y = 100, #number of years in scenario
	                 VBGF = list(model="Lknot",Linf = 80.4473850, K = 0.1396317, Lknot=4.7363329, Amat = 26.4950786), #VBGF parameters
	                 TD_MVN = TD_CC_MVN, #take demographics MVN parameters
	                 PjV = list(mu = 0.8, sd = 0.031), #juvenile survival
	                 PaV = list(mu = 0.895, sd = 0.028), #adult survival
	                 RI.dist = "normal",
	                 RV = list(RI = 3.3, RI_sd = 2.3, CF = 4.6, CF_sd = 1.1, CS = 122, CS_sd = 18.4, PF = 0.65, PF_sd = 0), #reproductive values
	                 ATLp = list(mu1 = 3.444104, nu1 = 0.06451325, mu2 = 7.506372e-05, nu2 =  0.01686626, mu3 = NULL, nu3 = NULL),
	                 ANE_dyn = FALSE, #whether ANE is dynamic
	                 Surv_dyn = FALSE,
	                 dynUQ = FALSE,
	                 ATL_scale = FALSE,
	                 grim.reaper=TRUE)
		
	###-----------------------------------------------------
	#		projection function
	###-----------------------------------------------------
			
		# Toggles:
			# dynUQ: toggles whether a static U and Q to a dynamic (annual) U and Q is used
			# ATL_scale: toggles whether the ATL scales along with the population

		# Included modes:
			# det. ... deterministic runs
			# sto. ... stochastic runs (RI, Pj)
			# by 
			# nost ... sex and discard mortality are drawn with binomial

		
		### DC
			DC.l <- sim.fn(dc.trend, DC.scenario, 10000, fn=pva.proj)
			save(DC.l, file=paste0(sim.path,"DC_median_sim.Rdata"))
			DC.l.low <- sim.fn(dc.trend.low, DC.scenario, 10000, fn=pva.proj)
			save(DC.l.low, file=paste0(sim.path,"DC_low_sim.Rdata"))
			DC.l.high <- sim.fn(dc.trend.high, DC.scenario, 10000, fn=pva.proj)
			save(DC.l.high, file=paste0(sim.path,"DC_high_sim.Rdata"))
		####

		### CC
			CC.l <- sim.fn(cc.trend, CC.scenario, 10000, fn=pva.proj)
			save(CC.l, file=paste0(sim.path,"CC_sim.Rdata"))
		#####

		#-------
		# Grim reaper
		
		### DC 
			DC.l <- sim.fn(dc.trend, DC.scenario.GR, 10000, fn=pva.proj)
			save(DC.l, file=paste0(sim.path,"DC_median_sim.GR.Rdata"))
			DC.l.low <- sim.fn(dc.trend.low, DC.scenario.GR, 10000, fn=pva.proj)
			save(DC.l.low, file=paste0(sim.path,"DC_low_sim.GR.Rdata"))
			DC.l.high <- sim.fn(dc.trend.high, DC.scenario.GR, 10000, fn=pva.proj)
			save(DC.l.high, file=paste0(sim.path,"DC_high_sim.GR.Rdata"))
		####

		### CC
			CC.l <- sim.fn(cc.trend, CC.scenario.GR, 10000, fn=pva.proj)
			save(CC.l, file=paste0(sim.path,"CC_sim.GR.Rdata"))
		#####
	###-----------------------------------------------------
	#		projection summary
	###-----------------------------------------------------
		thres.pop <- c(.5, .25, .125)
		thres.yr <- c(5, 10, 25, 50, 100)

		
		DC.sum <- proj.summ.fn(DC.l, dc.trend, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="DC", mode="MEDIAN")
		DC.sum.low <- proj.summ.fn(DC.l.low, dc.trend.low, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="DC", mode="LOW")
		DC.sum.high <- proj.summ.fn(DC.l.high, dc.trend.high, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="DC", mode="HIGH")
		CC.sum <- proj.summ.fn(CC.l, cc.trend, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="CC")

		#grim reaper
		DC.sum <- proj.summ.fn(DC.l, dc.trend, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="DC", mode="MEDIAN.GR")
		DC.sum.low <- proj.summ.fn(DC.l.low, dc.trend.low, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="DC", mode="LOW.GR")
		DC.sum.high <- proj.summ.fn(DC.l.high, dc.trend.high, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="DC", mode="HIGH.GR")
		CC.sum <- proj.summ.fn(CC.l, cc.trend, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp="CC", mode="GR")

	###-----------------------------------------------------
	#		plot raw
	###-----------------------------------------------------
		load(paste0(data.path,"DC_raw_nests_month.Rdata"))

		pdf(file=paste0(fig.path,"Raw_nest_data_monthly.pdf"))
			layout(matrix(c(rep(1:3,each=3),4,4,4,5,0,6,6,6,7),ncol=2))

			par(mar=c(2,5,0.1,1),cex.axis=1.2, oma=c(1,1,2,3))
			with(CC.dat, plot(Year, Maehama.Beach, ylim=range(pretty(range(CC.dat[,2:4], na.rm=T))), xlab="", ylab="", type="b", pch=16, las=1, xaxt='n'))
			axis(1, pretty(CC.dat$Year), labels=FALSE)
			mtext("N. Pac. Loggerheads", side=3, line=0.5)
			fig.lab("Maehama", xscale=0.025, yscale=0.975, adj=c(0,1))
			with(CC.dat, plot(Year, Inakahama.Beach, ylim=range(pretty(range(CC.dat[,2:4], na.rm=T))), xlab="", ylab="", type="b", pch=16, las=1, xaxt='n'))
			fig.lab("Inakahama", xscale=0.025, yscale=0.975, adj=c(0,1))
			axis(1, pretty(CC.dat$Year), labels=FALSE)
			with(CC.dat, plot(Year, Yakushima.Yotsuse, ylim=range(pretty(range(CC.dat[,2:4], na.rm=T))), xlab="", ylab="", type="b", pch=16, las=1, xaxt='n'))
			fig.lab("Yotsuse", xscale=0.025, yscale=0.975, adj=c(0,1))
			axis(1, pretty(CC.dat$Year), labels=TRUE, xpd=NA)

			### DC
			with(DC.raw.dat, {
				plot(Season, median.JM, ylim=range(pretty(range(DC.raw.dat[,2:7], na.rm=T))), xlab="", ylab="", type="n", pch=16, las=1, xaxt='n', xlim=range(Season))
				polygon(c(Season, rev(Season)),
				        c(low.JM, rev(high.JM)),
				        col=col2rgbA('gray50',0.25), border='gray50')
				points(Season, median.JM, type='b', pch=16)})
			fig.lab("Jamursba-Medi", xscale=0.975, yscale=0.975, adj=c(1,1))
			mtext("W. Pac. Leatherbacks", side=3, line=0.5)
			axis(1, pretty(DC.raw.dat$Season), labels=FALSE)


			plot(y.raw$year.frac, y.raw$JM, ylim=range(pretty(range(y.raw[,1:2],na.rm=T))), type='l', col='black', pch=16, las=1, xlab="", ylab="", xlim=range(DC.raw.dat$Season), xaxt='n', yaxt='n', lwd=1.2)
			axis(1, pretty(DC.raw.dat$Season), labels=FALSE)
			axis(2, pretty(range(y.raw[,1:2],na.rm=T))[seq(1,6,by=2)], las=1)
			fig.lab("Monthly", xscale=0.99, yscale=0.95, adj=c(1,1), cex=1)

			
			with(DC.raw.dat, {
				plot(Season, median.W, ylim=range(pretty(range(DC.raw.dat[,2:7], na.rm=T))), xlab="", ylab="", type="n", pch=16, las=1, xaxt='n')
				polygon(c(Season, rev(Season)),
				        c(low.W, rev(high.W)),
				        col=col2rgbA('gray50',0.25), border='gray50')
				points(Season, median.W, type='b', pch=16)})
			fig.lab("Wermon", xscale=0.975, yscale=0.975, adj=c(1,1))
			axis(1, pretty(DC.raw.dat$Season), labels=FALSE)

			plot(y.raw$year.frac, y.raw$W, ylim=range(pretty(range(y.raw[,1:2],na.rm=T))), type='l', col='black', pch=16, las=1, xlab="", ylab="", xlim=range(DC.raw.dat$Season), xaxt='n', yaxt='n', lwd=1.2)
			fig.lab("Monthly", xscale=0.99, yscale=0.95, adj=c(1,1), cex=1)
			axis(1, pretty(DC.raw.dat$Season), labels=TRUE)
			axis(2, pretty(range(y.raw[,1:2],na.rm=T))[seq(1,6,by=2)], las=1)

			mtext("# of Nests", side=2, outer=T, line=-0.75, cex=1.4)
		dev.off()















