data.extract <- function(location, year.begin,year.end, season.begin = year.begin, season.end = year.end, file.path = NULL){
	# In March 2019, we received new data for 2018. So, the raw data file
	# has been updated.  
	# On 16 April 2019, the last few data points for 2019 were received
	# so the data files have been updated. 
	if (is.null(season.begin)) season.begin <- year.begin
	if (is.null(season.end)) season.end <- year.end
	
	if (location == "JM"){
		data.0 <- read.csv(paste0(file.path,"JM_nests_April2019.csv"))
		
		data.0 %>% 
			select(Year_begin, Month_begin, JM_Nests) %>%
			mutate(Nests = JM_Nests) -> data.0
		
	} else if (location == "W"){
		data.0 <- read.csv(paste0(file.path,"W_nests_April2019.csv"))
		data.0 %>% 
			select(Year_begin, Month_begin, W_Nests) %>%
			mutate(Nests = W_Nests) -> data.0
	}
	
	# create regularly spaced time series:
	data.2 <- data.frame(Year = rep(min(data.0$Year_begin, na.rm = T):max(data.0$Year_begin, na.rm = T),each = 12),
											 Month_begin = rep(1:12, max(data.0$Year_begin, na.rm = T) - min(data.0$Year_begin, na.rm = T) + 1)) %>%
		mutate(begin_date = as.Date(paste(Year, Month_begin, '01', sep = "-"), format = "%Y-%m-%d"), Frac.Year = Year + (Month_begin-0.5)/12) %>%
		select(Year, Month_begin, begin_date, Frac.Year)
	
	# also make "nesting season" that starts April and ends March
	
	data.0 %>% mutate(begin_date = as.Date(paste(Year_begin, Month_begin, '01', sep = "-"),format = "%Y-%m-%d")) %>%
		mutate(Year = Year_begin,
					 Month = Month_begin,
					 f_month = as.factor(Month),
					 f_year = as.factor(Year),
					 Frac.Year = Year + (Month_begin-0.5)/12) %>%
		select(Year, Month, Frac.Year, begin_date, Nests) %>%
		na.omit() %>%
		right_join(.,data.2, by = "begin_date") %>%
		transmute(Year = Year.y,
							Month = Month_begin,
							Frac.Year = Frac.Year.y,
							Nests = Nests,
							Season = ifelse(Month < 4, Year-1, Year),
							Seq.Month = ifelse(Month < 4, Month + 9, Month - 3)) %>%
		reshape::sort_df(.,vars = "Frac.Year") %>%
		filter(Season >= season.begin & Season <= season.end) -> data.1
	
	data.1 %>% filter(Month > 3 & Month < 10) -> data.summer
	data.1 %>% filter(Month > 9 | Month < 4) %>%
		mutate(Seq.Month = Seq.Month - 6) -> data.winter
	
	jags.data <- list(y = log(data.1$Nests),
										m = data.1$Seq.Month,
										T = nrow(data.1))
	
	y <- matrix(log(data.1$Nests), ncol = 12, byrow = TRUE)
	
	jags.data2 <- list(y = y,
										 m = matrix(data.1$Seq.Month, ncol = 12, byrow = TRUE),
										 n.years = nrow(y))
	
	y.summer <- matrix(log(data.summer$Nests),
										 ncol = 6, byrow = TRUE)
	
	y.winter <- matrix(log(data.winter$Nests),
										 ncol = 6, byrow = TRUE)
	
	jags.data2.summer <- list(y = y.summer,
														m = matrix(data.summer$Seq.Month,ncol = 6, byrow = TRUE),
														n.years = nrow(y.summer))
	
	jags.data2.winter <- list(y = y.winter,
														m = matrix(data.winter$Seq.Month, ncol = 6, byrow = TRUE),
														n.years = nrow(y.winter))
	
	out <- list(jags.data = jags.data,
							jags.data2 = jags.data2,
							jags.data.summer = jags.data2.summer,
							jags.data.winter = jags.data2.winter,
							data.1 = data.1,
							data.summer = data.summer,
							data.winter = data.winter)
	return(out)
}
sum.fn <- function(U){
	meanU <- mean(U)
	varU <- var(U)
	qU <- quantile(U, probs=c(0.5,0.025, 0.975))

	lamb <- exp(U)
	meanL <- mean(lamb)
	varL <- var(lamb)
	qL <- quantile(lamb, probs=c(0.5, 0.025, 0.975))

	summ <- matrix(c(meanU, qU[1], varU, qU[2:3], meanL, qL[1], varL, qL[2:3]), ncol=1)
	rownames(summ) <- c("meanU", "medianU", "varU", "L95U", "U95U", "meanl", "medianl", "varl", "L95l", "U95l")

	return(summ)
}
pva.proj <- function(sim, trend, scenario){
	with(scenario, {
	

	#convert U to lambda
		lambda <- rep(1, n_y) #exp of U
		Q.sd <- rep(NA, n_y) #std. dev. of Q (is variance)
		#toggle for dynamic or static U & Q
		if(dynUQ){
			samp.sim <- sample(1:length(trend$U), n_y-1, replace = FALSE)
			lambda[1] <- exp(trend$U[sim])
			lambda[2:n_y] <- exp(trend$U[samp.sim])
			Q.sd[1] <- sqrt(trend$Q[sim])
			Q.sd[2:n_y] <- sqrt(trend$Q[samp.sim])
		}else{
			lambda[1:n_y] <- exp(trend$U[sim])
			Q.sd[1:n_y] <- sqrt(trend$Q[sim])
		}
		
	###-----------------------------------------------------
	#draw a ATL
		#ATL_scale is a logical flag for scaling the ATL at the same rate as the population growth rate
		if(ATL_scale){
			ATL_segments <- matrix(0, n_y, 3)
			for(i in 1:n_y){
				if(!is.null(ATLp$mu1)){
				ATL_segments[i,1] <- rCMP(1, mu = ATLp$mu1*cumprod(lambda)[i], nu = ATLp$nu1)
				}
				if(!is.null(ATLp$mu2)){
					ATL_segments[i,2] <- rCMP(1, mu = ATLp$mu2*cumprod(lambda)[i], nu = ATLp$nu2)
				}
				if(!is.null(ATLp$mu3)){
					ATL_segments[i,3] <- rCMP(1, mu = ATLp$mu3*cumprod(lambda)[i], nu = ATLp$nu3)
				}
			}
			
			ATL <- rowSums(ATL_segments)
		}else{
			ATL_segments <- matrix(0, n_y, 3)
			if(!is.null(ATLp$mu1)){
				ATL_segments[,1] <- rCMP(n_y, mu = ATLp$mu1, nu = ATLp$nu1)
			}
			if(!is.null(ATLp$mu2)){
				ATL_segments[,2] <- rCMP(n_y, mu = ATLp$mu2, nu = ATLp$nu2)
			}
			if(!is.null(ATLp$mu3)){
				ATL_segments[,3] <- rCMP(n_y, mu = ATLp$mu3, nu = ATLp$nu3)
			}
			ATL <- rowSums(ATL_segments)
		}

	#draw a Size and Mortality based on ATL
		mu.l <- TD_MVN$beta0 + TD_MVN$beta1*ATL #linear model for mean lengths in take (log)
		mu.m <- rep(TD_MVN$mu0, n_y) # mean discard mortality in take (logit)

		pred.Age <- seq(0,100, by=0.01) #sequence of ages to back calculate over
		Lpred <- ifelse(VBGF$model =="tknot",with(VBGF, Linf * (1-exp(-K*(pred.Age-tknot)))), with(VBGF, Linf - (Linf - Lknot)*(1-exp(-K*pred.Age)))) #lengths over those ages

		det.TD.draw.raw <- rep(list(NA),n_y) #storage list of MVN draws
		det.TD.draw <- rep(list(NA), n_y) #storage list
		det.Take <- rep(0, n_y) #storage vector of take
		det.Take.nost <- rep(0, n_y) #storage vector of take

		sto.TD.draw.raw <- rep(list(NA),n_y) #storage list of MVN draws
		sto.TD.draw <- rep(list(NA), n_y) #storage list
		sto.Take <- rep(0, n_y) #storage vector of take
		sto.Take.nost <- rep(0, n_y) #storage vector of take

		#draw a proportion female in the Take for each year
		PF.sto <- rnorm(n_y, RV$PF, RV$PF_sd)
		PF.det <- RV$PF

		for(i in 1:n_y){

			if(ATL[i] > 0){
				#draw from MVN log(length) and logit(discard mortality)
				det.TD.draw.raw[[i]] <- sto.TD.draw.raw[[i]] <- rmvnorm(ATL[i], c(mu.l[i], mu.m[i]), TD_MVN$cov)
				#convert to length and discard mortality (prob. space)
				det.TD.draw[[i]] <- sto.TD.draw[[i]] <- data.frame(l = exp(sto.TD.draw.raw[[i]][,1]), m = boot::inv.logit(sto.TD.draw.raw[[i]][,2]))
				#get an anticpated Age from the VBGF
				det.TD.draw[[i]]$Age <- sto.TD.draw[[i]]$Age <- sapply(sto.TD.draw[[i]]$l, function(x) {pred.Age[which.min(abs(x - Lpred))]})
				#calculate the number of years until maturity
				det.TD.draw[[i]]$YatLarge <- sto.TD.draw[[i]]$YatLarge <- VBGF$Amat - sto.TD.draw[[i]]$Age
				#determine the animals stage based on maturity
				det.TD.draw[[i]]$Stage <- sto.TD.draw[[i]]$Stage <- ifelse(sto.TD.draw[[i]]$Age > VBGF$Amat, "A", "J")

				#### Stochastic
				#draw a remigration interval for every animal based on a specified distribution (normal or CMP)
				if(RI.dist =="normal"){
					sto.TD.draw[[i]]$RI <- rtruncnorm(nrow(sto.TD.draw[[i]]), a=0, mean=RV$RI, sd=RV$RI_sd) #remigration interval
				}else if(RI.dist =="CMP"){
					sto.TD.draw[[i]]$RI <- rCMP(nrow(sto.TD.draw[[i]]), mu=RV$RI, nu=RV$RI_sd) #remigration interval
					while(any(sto.TD.draw[[i]]$RI==0)){
						sto.TD.draw[[i]]$RI[sto.TD.draw[[i]]$RI==0] <- rCMP(sum(sto.TD.draw[[i]]$RI==0), mu = RV$RI, nu = RV$RI_sd)
					}
				}else{
					sto.TD.draw[[i]]$RI <- RV$RI
					#error handling
					print("Distribution was not specified as normal or CMP")
				}
				#draw a juvenile survival for every animal
				sto.TD.draw[[i]]$Pj <- rnorm(nrow(sto.TD.draw[[i]]), PjV$mu, PjV$sd) #juvenile survival
				#Calculate the Adult Nester Equivalency
				sto.TD.draw[[i]]$ANEj <- (sto.TD.draw[[i]]$Pj^sto.TD.draw[[i]]$YatLarge) * (1/sto.TD.draw[[i]]$RI) #ANE
				#Draw whether the animal is M/F based on yearly proportion
				sto.TD.draw[[i]]$Sex <- rbinom(nrow(sto.TD.draw[[i]]), 1, PF.sto[i])
				#Draw whether the animal lived or died based on discard mortality
				if(grim.reaper){
					sto.TD.draw[[i]]$Fdead <- 1
				}else{
					sto.TD.draw[[i]]$Fdead <- rbinom(nrow(sto.TD.draw[[i]]), 1, sto.TD.draw[[i]]$m)
				}

				#Fill in Adult ANE
				sto.TD.draw[[i]]$ANEa <- 1
				#Calculate the realized ANE
				sto.TD.draw[[i]]$rANE <- ((sto.TD.draw[[i]]$ANEj * as.integer(sto.TD.draw[[i]]$Stage=="J")) + (sto.TD.draw[[i]]$ANEa * as.integer(sto.TD.draw[[i]]$Stage=="A"))) * (sto.TD.draw[[i]]$Sex * sto.TD.draw[[i]]$Fdead)
				#Calculate the realized ANE w/o stochastic Sex or Discard Mortality draws
				sto.TD.draw[[i]]$rANE.nost <- ((sto.TD.draw[[i]]$ANEj * as.integer(sto.TD.draw[[i]]$Stage=="J")) + (sto.TD.draw[[i]]$ANEa * as.integer(sto.TD.draw[[i]]$Stage=="A"))) * (PF.sto[i] * sto.TD.draw[[i]]$m)
				
					
				
				#### Deterministic
				det.TD.draw[[i]]$RI <- RV$RI
				det.TD.draw[[i]]$Pj <- PjV$mu
				det.TD.draw[[i]]$Pa <- PaV$mu
				det.TD.draw[[i]]$ANEj <- (det.TD.draw[[i]]$Pj^det.TD.draw[[i]]$YatLarge) * (1/det.TD.draw[[i]]$RI)
				det.TD.draw[[i]]$Sex <- rbinom(length(det.TD.draw[[i]]$Stage), 1, PF.det)
				if(grim.reaper){
					det.TD.draw[[i]]$Fdead <- 1
				}else{
					det.TD.draw[[i]]$Fdead <- rbinom(nrow(det.TD.draw[[i]]), 1, det.TD.draw[[i]]$m)
				}
				
				det.TD.draw[[i]]$ANEa <- 1
				
				det.TD.draw[[i]]$rANE <- ((det.TD.draw[[i]]$ANEj * as.integer(det.TD.draw[[i]]$Stage=="J")) + (det.TD.draw[[i]]$ANEa * as.integer(det.TD.draw[[i]]$Stage=="A"))) * (det.TD.draw[[i]]$Sex * det.TD.draw[[i]]$Fdead)
				det.TD.draw[[i]]$rANE.nost <- ((det.TD.draw[[i]]$ANEj * as.integer(det.TD.draw[[i]]$Stage=="J")) + (det.TD.draw[[i]]$ANEa * as.integer(det.TD.draw[[i]]$Stage=="A"))) * (PF.det * det.TD.draw[[i]]$m)
				
				

				det.Take[i] <- sum(det.TD.draw[[i]]$rANE)
				det.Take.nost[i] <- sum(det.TD.draw[[i]]$rANE.nost)
				sto.Take[i] <- sum(sto.TD.draw[[i]]$rANE)
				sto.Take.nost[i] <- sum(sto.TD.draw[[i]]$rANE.nost)
			}
		}
		
		### deterministic
		det.nost.Nt_take <- det.Nt_take <- det.Nt <- det.nost.Nt_take.iQ <- det.Nt_take.iQ <- det.Nt.iQ <- rep(0, n_y)
		
		det.Nt_take.iQ[1] <- rnorm(1,(trend$N0[sim] - det.Take[1])*lambda[1], Q.sd[1])
		det.nost.Nt_take.iQ[1] <- rnorm(1,(trend$N0[sim] - det.Take.nost[1])*lambda[1], Q.sd[1])
		det.Nt.iQ[1] <- rnorm(1, trend$N0[sim] * lambda[1], Q.sd[1])
		det.Nt_take[1] <- (trend$N0[sim] - det.Take[1])*lambda[1]
		det.nost.Nt_take[1] <- (trend$N0[sim] - det.Take.nost[1])*lambda[1]
		det.Nt[1] <-  trend$N0[sim] * lambda[1]
		for(i in 2:n_y){
			#include process variance Q
			det.Nt_take.iQ[i] <- rnorm(1,(det.Nt_take.iQ[i-1] - det.Take[i])*lambda[i], Q.sd[i])
			det.nost.Nt_take.iQ[i] <- rnorm(1,(det.Nt_take.iQ[i-1] - det.Take.nost[i])*lambda[i], Q.sd[i])
			det.Nt.iQ[i] <- rnorm(1, det.Nt.iQ[i-1] * lambda[i], Q.sd[i])
			#drop process variance Q
			det.Nt_take[i] <- (det.Nt_take[i-1] - det.Take[i])*lambda[i]
			det.nost.Nt_take[i] <- (det.Nt_take[i-1] - det.Take.nost[i])*lambda[i]
			det.Nt[i] <-  det.Nt[i-1] * lambda[i]
		}

		### stochastic
		sto.nost.Nt_take <- sto.nost.Nt <- sto.Nt_take <- sto.Nt <- sto.nost.Nt_take.iQ <- sto.nost.Nt.iQ <-sto.Nt_take.iQ <- sto.Nt.iQ <- rep(0, n_y)

		sto.Nt_take.iQ[1] <- rnorm(1,(trend$N0[sim] - sto.Take[1])*lambda[1], Q.sd[1])
		sto.nost.Nt_take.iQ[1] <- rnorm(1,(trend$N0[sim] - sto.Take.nost[1])*lambda[1], Q.sd[1])
		sto.Nt.iQ[1] <- rnorm(1, trend$N0[sim] * lambda[1], Q.sd[1])
		sto.Nt_take[1] <- (trend$N0[sim] - sto.Take[1])*lambda[1]
		sto.nost.Nt_take[1] <- (trend$N0[sim] - sto.Take.nost[1])*lambda[1]
		sto.Nt[1] <-  trend$N0[sim] * lambda[1]
		for(i in 2:n_y){
			#include process variance Q
			sto.Nt_take.iQ[i] <- rnorm(1,(sto.Nt_take.iQ[i-1] - sto.Take[i])*lambda[i], Q.sd[i])
			sto.nost.Nt_take.iQ[i] <- rnorm(1,(sto.Nt_take.iQ[i-1] - sto.Take.nost[i])*lambda[i], Q.sd[i])
			sto.Nt.iQ[i] <- rnorm(1, sto.Nt.iQ[i-1] * lambda[i], Q.sd[i])
			#drop process variance Q
			sto.Nt_take[i] <- (sto.Nt_take[i-1] - sto.Take[i])*lambda[i]
			sto.nost.Nt_take[i] <- (sto.Nt_take[i-1] - sto.Take.nost[i])*lambda[i]
			sto.Nt[i] <-  sto.Nt[i-1] * lambda[i]
		}

	return(as.matrix(
	 data.frame(det.nost.Nt_take = det.nost.Nt_take,
	            det.Nt_take = det.Nt_take,
	            det.Nt = det.Nt, 
	            det.nost.Nt_take.iQ = det.nost.Nt_take.iQ,
	            det.Nt_take.iQ = det.Nt_take.iQ,
	            det.Nt.iQ = det.Nt.iQ,
	            sto.nost.Nt_take = sto.nost.Nt_take,
	            sto.Nt_take = sto.Nt_take,
	            sto.Nt = sto.Nt, 
	            sto.nost.Nt_take.iQ = sto.nost.Nt_take.iQ,
	            sto.Nt_take.iQ = sto.Nt_take.iQ,
	            sto.Nt.iQ = sto.Nt.iQ)))
	})
}
curr.abund.fn <- function(dat, RI, round=TRUE){
	curr.abund <- matrix(c(quantile(dat$N_fym0, probs=c(0.5,0.025,0.975)), quantile(dat$N_fym1, probs=c(0.5,0.025,0.975)), quantile(dat$N_fym2, probs=c(0.5,0.025,0.975)), quantile(dat$N_fym3, probs=c(0.5,0.025,0.975)), quantile(rowSums(dat[,c("N_fym0","N_fym1","N_fym2","N_fym3")]),probs=c(0.5,0.025,0.975))*(RI/4)), nrow=5, byrow=T)
	rownames(curr.abund) <- c("N0","N-1","N-2","N-3","Sum")
	colnames(curr.abund) <- c("Median","L95%","U95%")
	if(round) curr.abund <- round(curr.abund)
	return(curr.abund)
}
sim.fn <- function(trend, scenario, nsim, fn){
	pva.proj <- fn
	UseCores <- detectCores() - 2
	cl <- makeCluster(UseCores)
	registerDoParallel(cl)
	nst <- ceiling(seq(from=1, to = nsim, length.out=(UseCores+1)))[-(UseCores+1)]
	nen <- ceiling(seq(from=1, to = nsim, length.out=(UseCores+1)))[-1] - c(rep(1,(UseCores-1)),0)
	
	df <- foreach(i = 1:length(nst), .packages=c("CMP", "mvtnorm","truncnorm")) %dopar% {

		iseq <- seq(nst[i],nen[i])
		arr <- array(NA, dim=c(100, 12, length(iseq)))

		for(j in 1:length(iseq)){
			arr[,,j] <- pva.proj(iseq[j], trend, scenario)
		}

		return(arr)
	}
	stopCluster(cl)
	return(df)
}
bootCI <- function(dt, R){
	samps <- replicate(R, sample(1:length(dt), floor(0.9*length(dt))))
	q <- apply(samps, 2, function(v) {sum(dt[v]/length(v))})
	return(q)
}
fig.lab <- function(i, xscale=0.05, yscale, cex=1.4, adj=c(0.5,0.5)){
	text(x = par('usr')[1] + abs(diff(par('usr')[1:2]))*xscale,
	     y = par('usr')[3] + abs(diff(par('usr')[3:4]))*yscale,
	     ifelse(is.numeric(i),LETTERS[i],i), xpd=NA, cex=cex, adj=adj)
}
u.den <- function(notake, take){
	d1 <- density(notake, adj=2)
	d2 <- density(take, adj=2)
	d1$y <- d1$y/max(d1$y)
	d2$y <- d2$y/max(d2$y)

	plot(d1$x, d1$y, xlab="r", ylab="Relative Density", ylim=c(0,1.01), las=1, col='chartreuse4', type='l', lwd=3)
	lines(d2$x, d2$y, col='dodgerblue4', lwd=3, lty=2)
	legend("topleft", legend=c(expression(N[j]), expression(N[j]-F)), lwd=3, col=c('chartreuse4', 'dodgerblue4'), bty='n')
}
col2rgbA<-function(color,transparency)
{
  rgb(t(col2rgb(color))/255,alpha=transparency)
}
proj.summ.fn <- function(sim.l, trend, thres.pop, thres.yr, alpha=0.05, keepers=c(5,6,11,12), spp, mode=NULL){

	lci <- alpha/2
	uci <- 1-(alpha/2)

	arr <- abind(sim.l, along=3)

	### Summary arrays
	#dims = length(thres.pop), 12 modes from pva.proj, 4 summary metrics (mean, median, LCI, UCI)
	tab.yr <- array(NA, dim=c(length(thres.pop),12, 4))
	#dims = length(thres.yr), length(thres.pop), 12 modes from pva.proj, 3 summary metrics (median, LCI, UCI)
	tab.prob <- array(NA, dim=c(length(thres.yr), length(thres.pop), 12, 3))

	### Storage arrays
	#dims = length(thres.pop), 12 modes from pva.proj, nsims
	yr <- array(NA, dim=c(length(thres.pop), 12, dim(arr)[3]))
	#dims = length(thres.yr), length(thres.pop), 12 modes from pva.proj, nsims
	prob <- array(NA, dim=c(length(thres.yr), length(thres.pop), 12, dim(arr)[3]))

	#duplicate array to fill in extinct years with 0's
	arr2 <- arr

	for(i in 1:dim(arr)[3]){
		N0 <- trend$N0[i]
		N0.thres <- N0*thres.pop

		#fills in zeros for extinct runs
		arr2[,,i] <- apply(arr2[,,i], 2, function(x) {if(any(x < 0)){ d <- min(which(x < 0)); x[d:100] <- 0}; return(x)})

		#first year to fall below threshold
		yr[,,i] <- apply(arr2[,,i], 2, FUN = function(x) {sapply(N0.thres, function(v) ifelse(min(which((x < v)==1))==Inf, NA, min(which((x < v)==1))))})

		# prob[,,,i] <- array(apply(arr2[,,i], 2, function(x) {sapply(N0.thres, function(v){ sapply(thres.yr, function(y) sum(x[1:y] < v)/y)})}), dim=c(5,3,12))
		prob[,,,i] <- array(apply(arr2[,,i], 2, function(x) {sapply(N0.thres, function(v){ sapply(thres.yr, function(y) as.integer(any(x[1:y] < v)))})}), dim=c(5,3,12))
	}

	#Summarize Years to threshold
	tab.yr[,,1] <- apply(yr, c(1,2), mean, na.rm=T)
	tab.yr[,,2] <- apply(yr, c(1,2), median, na.rm=T)
	tab.yr[,,3] <- apply(yr, c(1,2), quantile, prob = lci, na.rm=T)
	tab.yr[,,4] <- apply(yr, c(1,2), quantile, prob = uci, na.rm=T)

	#Summarize probability of falling below threshold
	tab.prob[,,,1] <- apply(prob, c(1,2,3), function(x) {sum(x)/length(x)})
	bootci.prob <- apply(prob,c(1,2,3), bootCI, R=1000)
	tab.prob[,,,2] <- apply(bootci.prob, c(2,3,4), quantile, prob=lci)
	tab.prob[,,,3] <- apply(bootci.prob, c(2,3,4), quantile, prob=uci)

	# Table 1
		tab1.det.NT <- matrix(c(1-aperm(tab.prob[,,keepers[2],], c(2,1,3))[,5,1], aperm(tab.prob[,,6,], c(2,1,3))[,5,1], tab.yr[,6,]), nrow=3)
		tab1.det.T <- matrix(c(1-aperm(tab.prob[,,keepers[1],], c(2,1,3))[,5,1], aperm(tab.prob[,,5,], c(2,1,3))[,5,1], tab.yr[,5,]), nrow=3)
		tab1.sto.NT <- matrix(c(1-aperm(tab.prob[,,12,], c(2,1,3))[,5,1], aperm(tab.prob[,,keepers[4],], c(2,1,3))[,5,1], tab.yr[,12,]), nrow=3)
		tab1.sto.T <- matrix(c(1-aperm(tab.prob[,,11,], c(2,1,3))[,5,1], aperm(tab.prob[,,keepers[3],], c(2,1,3))[,5,1], tab.yr[,11,]), nrow=3)

		rownames(tab1.det.NT) <- rownames(tab1.det.T) <- rownames(tab1.sto.NT) <- rownames(tab1.sto.T) <- paste0(thres.pop*100, "%")
		colnames(tab1.det.NT) <- colnames(tab1.det.T) <- colnames(tab1.sto.NT) <- colnames(tab1.sto.T) <- c("Prob.Above","Prob.Below","MeanYr","MedYr","L95Yr","U95Yr")

		write.csv(round(tab1.det.NT,2), file=paste0(table.path,"Table1.",spp, mode,".det.NT.csv"))
		write.csv(round(tab1.det.T,2), file=paste0(table.path,"Table1.",spp, mode,".det.T.csv"))
		write.csv(round(tab1.sto.NT,2), file=paste0(table.path,"Table1.",spp, mode,".sto.NT.csv"))
		write.csv(round(tab1.sto.T,2), file=paste0(table.path,"Table1.",spp, mode,".sto.T.csv"))

	# Table 2
		tp.det.NT <- apply(aperm(tab.prob[,,keepers[2],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]
		tp.det.T <- apply(aperm(tab.prob[,,keepers[1],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]
		tp.sto.NT <- apply(aperm(tab.prob[,,keepers[4],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]
		tp.sto.T <- apply(aperm(tab.prob[,,keepers[3],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]

		rownames(tp.det.NT) <- rownames(tp.det.T) <- rownames(tp.sto.NT) <- rownames(tp.sto.T) <- paste0(rep(paste0(thres.pop*100,"%"),each=3),  c("", "L95","U95"))
		colnames(tp.det.NT) <- colnames(tp.det.T) <- colnames(tp.sto.NT) <- colnames(tp.sto.T) <- paste0(thres.yr,'yr')

		write.csv(round(tp.det.NT,3), file=paste0(table.path,"Table2.",spp, mode,".det.NT.csv"))
		write.csv(round(tp.det.T,3), file=paste0(table.path,"Table2.",spp, mode,".det.T.csv"))
		write.csv(round(tp.sto.NT,3), file=paste0(table.path,"Table2.",spp, mode,".sto.NT.csv"))
		write.csv(round(tp.sto.T,3), file=paste0(table.path,"Table2.",spp, mode,".sto.T.csv"))

	#Median Projections
		proj.med <- apply(arr2, c(1,2), quantile, probs=c(lci, 0.5, uci))
		proj.log.med <- apply(arr2, c(1,2), function(x) {quantile(log(x), probs=c(lci, 0.5, uci))})
		proj.log.med[is.infinite(proj.log.med)] <- min(proj.log.med[!is.infinite(proj.log.med)]) * 5
		diff.det <- arr2[,keepers[1],] - arr2[,keepers[2],]
		diff.sto <- arr2[,keepers[3],] - arr2[,keepers[4],]
		proj.diff.det <- apply(diff.det, 1, quantile, probs=c(lci, 0.5, uci))
		proj.diff.sto <- apply(diff.sto, 1, quantile, probs=c(lci, 0.5, uci))
		if(spp=="DC"){
			pdf(file=paste0(fig.path,spp,mode,"proj100.pdf"))
				###PROJECT 100
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2018:2117, proj.med[2,,keepers[1]], 
				     xlim=c(2017,2118), 
				     ylim=c(0,max(proj.med[,,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.med[1,,keepers[1]], rev(proj.med[3,,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.med[1,,keepers[2]], rev(proj.med[3,,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2018:2117, proj.med[2,,keepers[1]], col="dodgerblue4", lwd=3)
				lines(2018:2117, proj.med[2,,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(2018:2118)[pretty(2018:2118)>2018 & pretty(2018:2118) < 2118]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2018:2117, proj.med[2,,keepers[3]], 
				     xlim=c(2017,2118), 
				     ylim=c(0,max(proj.med[,,keepers])), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.med[1,,keepers[3]], rev(proj.med[3,,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.med[1,,keepers[4]], rev(proj.med[3,,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2018:2117, proj.med[2,,keepers[3]], col="dodgerblue4", lwd=3)
				lines(2018:2117, proj.med[2,,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				legend("bottomright", legend=c(expression(paste(N[t] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"logproj100.pdf"))
				###PROJECT 100
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2018:2117, proj.log.med[2,,keepers[1]], 
				     xlim=c(2017,2118), 
				     ylim=c(0.0001,max(proj.log.med[,,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.log.med[1,,keepers[1]], rev(proj.log.med[3,,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.log.med[1,,keepers[2]], rev(proj.log.med[3,,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2018:2117, proj.log.med[2,,keepers[1]], col="dodgerblue4", lwd=3)
				lines(2018:2117, proj.log.med[2,,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(2018:2118)[pretty(2018:2118)>2018 & pretty(2018:2118) < 2118]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2018:2117, proj.log.med[2,,keepers[3]], 
				     xlim=c(2017,2118), 
				     ylim=c(0.0001,max(proj.log.med[,,keepers])), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.log.med[1,,keepers[3]], rev(proj.log.med[3,,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.log.med[1,,keepers[4]], rev(proj.log.med[3,,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2018:2117, proj.log.med[2,,keepers[3]], col="dodgerblue4", lwd=3)
				lines(2018:2117, proj.log.med[2,,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				legend("topright", legend=c(expression(paste(N[j] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"proj10.pdf"))
				#### PROJECT 10
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2018:2027, proj.med[2,1:10,keepers[1]], 
				     xlim=c(2017.75,2027.25), 
				     ylim=c(0,max(proj.med[,1:10,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2027, rev(2018:2027)), y = c(proj.med[1,1:10,keepers[1]], rev(proj.med[3,1:10,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2018:2027, rev(2018:2027)), y = c(proj.med[1,1:10,keepers[2]], rev(proj.med[3,1:10,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2018:2027, proj.med[2,1:10,keepers[1]], col="dodgerblue4", lwd=3)
				lines(2018:2027, proj.med[2,1:10,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(2018:2028)[pretty(2018:2028)>2018 & pretty(2018:2028) < 2028]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2018:2027, proj.med[2,1:10,keepers[3]], 
				     xlim=c(2017.75,2027.25), 
				     ylim=c(0,max(proj.med[,1:10,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2027, rev(2018:2027)), y = c(proj.med[1,1:10,keepers[3]], rev(proj.med[3,1:10,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2018:2027, rev(2018:2027)), y = c(proj.med[1,1:10,keepers[4]], rev(proj.med[3,1:10,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2018:2027, proj.med[2,1:10,keepers[3]], col="dodgerblue4", lwd=3)
				lines(2018:2027, proj.med[2,1:10,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				legend("bottomleft", legend=c(expression(paste(N[j] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2)
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"diff100.pdf"))
				###PROJECT 100
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2018:2117, proj.diff.det[2,], 
				     xlim=c(2017,2118), 
				     ylim=range(pretty(range(proj.diff.det))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.diff.det[1,], rev(proj.diff.det[3,])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2018:2117, proj.diff.det[2,], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2018:2118)[pretty(2018:2118)>2018 & pretty(2018:2118) < 2118]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2018:2117, proj.diff.sto[2,], 
				     xlim=c(2017,2118), 
				     ylim=range(pretty(range(proj.diff.sto))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2117, rev(2018:2117)), y = c(proj.diff.sto[1,], rev(proj.diff.sto[3,])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2018:2117, proj.diff.sto[2,], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2018:2118)[pretty(2018:2118)>2018 & pretty(2018:2118) < 2118]
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				legend("bottomleft", legend=c(expression(paste(N[j] - (N[j]-F)," Median")), expression(paste(N[j] - (N[j]-F)," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("darkorchid4",col2rgbA("darkorchid",0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"diff10.pdf"))
				#### PROJECT 10
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2018:2027, proj.diff.det[2,1:10], 
				     xlim=c(2017.75,2027.25), 
				     ylim=range(pretty(range(proj.diff.det[,1:10]))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2027, rev(2018:2027)), y = c(proj.diff.det[1,1:10], rev(proj.diff.det[3,1:10])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2018:2027, proj.diff.det[2,1:10], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2018:2028)[pretty(2018:2028)>2018 & pretty(2018:2028) < 2028]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2018:2027, proj.diff.sto[2,1:10], 
				     xlim=c(2017.75,2027.25), 
				     ylim=range(pretty(range(proj.diff.sto[,1:10]))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2018:2027, rev(2018:2027)), y = c(proj.diff.sto[1,1:10], rev(proj.diff.sto[3,1:10])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2018:2027, proj.diff.sto[2,1:10], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2018:2028)[pretty(2018:2028)>2018 & pretty(2018:2028) < 2028]
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				legend("bottomleft", legend=c(expression(paste(N[j] - (N[j]-F)," Median")), expression(paste(N[j] - (N[j]-F)," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("darkorchid4",col2rgbA("darkorchid",0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()

		}else{
			pdf(file=paste0(fig.path,spp,mode,"proj100.pdf"))
				###PROJECT 100
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2016:2115, proj.med[2,,keepers[1]], 
				     xlim=c(2015,2116), 
				     ylim=c(0,max(proj.med[2,,keepers])), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.med[1,,keepers[1]], rev(proj.med[3,,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.med[1,,keepers[2]], rev(proj.med[3,,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2016:2115, proj.med[2,,keepers[1]], col="dodgerblue4", lwd=3)
				lines(2016:2115, proj.med[2,,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(2016:2116)[pretty(2016:2116)>2016 & pretty(2016:2116) < 2116]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2016:2115, proj.med[2,,keepers[3]], 
				     xlim=c(2015,2116), 
				     ylim=c(0,max(proj.med[2,,keepers])), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.med[1,,keepers[3]], rev(proj.med[3,,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.med[1,,keepers[4]], rev(proj.med[3,,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2016:2115, proj.med[2,,keepers[3]], col="dodgerblue4", lwd=3)
				lines(2016:2115, proj.med[2,,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)
				legend("topleft", legend=c(expression(paste(N[j] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"logproj100.pdf"))
				###PROJECT 100
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2016:2115, proj.log.med[2,,keepers[1]], 
				     xlim=c(2015,2116), 
				     ylim=c(0.0001,max(proj.log.med[,,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.log.med[1,,keepers[1]], rev(proj.log.med[3,,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.log.med[1,,keepers[2]], rev(proj.log.med[3,,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2016:2115, proj.log.med[2,,keepers[1]], col="dodgerblue4", lwd=3)
				lines(2016:2115, proj.log.med[2,,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(2016:2116)[pretty(2016:2116)>2016 & pretty(2016:2116) < 2116]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2016:2115, proj.log.med[2,,keepers[3]], 
				     xlim=c(2015,2116), 
				     ylim=c(0.0001,max(proj.log.med[,,keepers])), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.log.med[1,,keepers[3]], rev(proj.log.med[3,,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.log.med[1,,keepers[4]], rev(proj.log.med[3,,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2016:2115, proj.log.med[2,,keepers[3]], col="dodgerblue4", lwd=3)
				lines(2016:2115, proj.log.med[2,,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)
				legend("topleft", legend=c(expression(paste(N[j] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2, inset=c(0,0.1))
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"proj10.pdf"))
				#### PROJECT 10
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2016:2025, proj.med[2,1:10,keepers[1]], 
				     xlim=c(2015.75,2025.25), 
				     ylim=c(0,max(proj.med[,1:10,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2025, rev(2016:2025)), y = c(proj.med[1,1:10,keepers[1]], rev(proj.med[3,1:10,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2016:2025, rev(2016:2025)), y = c(proj.med[1,1:10,keepers[2]], rev(proj.med[3,1:10,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2016:2025, proj.med[2,1:10,keepers[1]], col="dodgerblue4", lwd=3)
				lines(2016:2025, proj.med[2,1:10,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(2016:2026)[pretty(2016:2026)>2016 & pretty(2016:2026) < 2026]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2016:2025, proj.med[2,1:10,keepers[3]], 
				     xlim=c(2015.75,2025.25), 
				     ylim=c(0,max(proj.med[,1:10,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2025, rev(2016:2025)), y = c(proj.med[1,1:10,keepers[3]], rev(proj.med[3,1:10,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(2016:2025, rev(2016:2025)), y = c(proj.med[1,1:10,keepers[4]], rev(proj.med[3,1:10,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(2016:2025, proj.med[2,1:10,keepers[3]], col="dodgerblue4", lwd=3)
				lines(2016:2025, proj.med[2,1:10,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)
				legend("topleft", legend=c(expression(paste(N[j] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"diff100.pdf"))
				###PROJECT 100
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2016:2115, proj.diff.det[2,], 
				     xlim=c(2015,2116), 
				     ylim=range(pretty(range(proj.diff.det))), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.diff.det[1,], rev(proj.diff.det[3,])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2016:2115, proj.diff.det[2,], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2016:2116)[pretty(2016:2116)>2016 & pretty(2016:2116) < 2116]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2016:2115, proj.diff.sto[2,], 
				     xlim=c(2015,2116), 
				     ylim=range(pretty(range(proj.diff.sto))), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2115, rev(2016:2115)), y = c(proj.diff.sto[1,], rev(proj.diff.sto[3,])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2016:2115, proj.diff.sto[2,], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2016:2116)[pretty(2016:2116)>2016 & pretty(2016:2116) < 2116]
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)
				legend("topleft", legend=c(expression(paste(N[j] - (N[j]-F)," Median")), expression(paste(N[j] - (N[j]-F)," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("darkorchid4",col2rgbA("darkorchid",0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()

			pdf(file=paste0(fig.path,spp,mode,"diff10.pdf"))
				#### PROJECT 10
				par(mfrow=c(2,1), mar=c(3,4,2,1))
				# Deterministic
				plot(2016:2025, proj.diff.det[2,1:10], 
				     xlim=c(2015.75,2025.25), 
				     ylim=range(pretty(range(proj.diff.det[,1:10]))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2025, rev(2016:2025)), y = c(proj.diff.det[1,1:10], rev(proj.diff.det[3,1:10])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2016:2025, proj.diff.det[2,1:10], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2016:2026)[pretty(2016:2026)>2016 & pretty(2016:2026) < 2026]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)

				# Stochastic
				par(mar=c(4,4,1,1))
				plot(2016:2025, proj.diff.sto[2,1:10], 
				     xlim=c(2015.75,2025.25), 
				     ylim=range(pretty(range(proj.diff.sto[,1:10]))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(2016:2025, rev(2016:2025)), y = c(proj.diff.sto[1,1:10], rev(proj.diff.sto[3,1:10])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(2016:2025, proj.diff.sto[2,1:10], col="darkorchid4", lwd=3)
				yrpretty <- pretty(2016:2026)[pretty(2016:2026)>2016 & pretty(2016:2026) < 2026]
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("N. Pac. Loggerheads", xscale=0.15, yscale=0.925, cex = 1)
				legend("topleft", legend=c(expression(paste(N[j] - (N[j]-F)," Median")), expression(paste(N[j] - (N[j]-F)," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("darkorchid4",col2rgbA("darkorchid",0.5)), bty='n', pt.cex = 2, inset = c(0,0.1))
			dev.off()
		}
	
	return(list(tab.yr = tab.yr,
	            tab.prob = tab.prob))
}
cor.lab <- function(x, y, xscale, yscale, cex=1.2, adj=c(0.5,0.5)){
	xpt <- par('usr')[1]+abs(diff(par('usr')[1:2]))*xscale
	ypt <- par('usr')[3]+abs(diff(par('usr')[3:4]))*yscale

	cor <- round(cor(x,y),2)

	text(xpt, ypt, cor, xpd=NA, cex=cex, adj=adj)
}
kde.plot <- function(x, y, n=100){
	kde <- MASS::kde2d(x, y, n=n)
	kde$z <- kde$z/max(kde$z)
	contour(kde, levels=c(0.05, 0.5), col='black', drawlabels=FALSE, add=TRUE, lwd=c(1,3))
}
plot.joint <- function(trend, spp, graph){
	d.u <- density(trend$U, adj=2)
	d.N0 <- density(trend$N0, adj=2)
	d.Q <- density(trend$Q, adj=2)

	d.u$y <- d.u$y/max(d.u$y)
	d.N0$y <- d.N0$y/max(d.N0$y)
	d.Q$y <- d.Q$y/max(d.Q$y)

	samp <- sample(1:nrow(trend), floor(0.5*nrow(trend)))
	pdf(file=paste0(fig.path,spp,"joint_post.pdf"))
		par(mfrow=c(3,3), mar=c(4,4,1,1), cex.axis=1.1, cex.lab=1.1)
	for(i in 1:9){
		if(i ==1) {
			#1
			plot(d.u$x, d.u$y, type='l', xlab='U', ylab='Rel. Density', ylim=c(0,1), yaxt='n', las=1)
			axis(2, at=pretty(c(0,1)), las=1)
		}
		if(i==2){
			#2
			plot(trend$U, trend$N0, pch=".", col='grey80', xlab='U', ylab=expression(N[0]), las=1)
			kde.plot(trend$U[samp], trend$N0[samp], n=100)
		}
		if(i==3){
			#3
			plot(trend$U[samp], trend$Q[samp], pch=".", col='grey80',xlab='U', ylab='Q', las=1)
			kde.plot(trend$U, trend$Q, n=100)
		}
		if(i==4){
			#4
			plot.new()
			box()
			cor.lab(trend$U, trend$N0, xscale=0.5, yscale=0.5, cex=1.4)
			fig.lab(expression(paste("U, ",N[0])), xscale=0.15, yscale=0.1, cex=1.1)
		}
		if(i==5){
			#5
			plot(d.N0$x, d.N0$y, type='l', xlab=expression(N[0]), ylab='Rel. Density', ylim=c(0,1), yaxt='n', las=1)
			axis(2, at=pretty(c(0,1)), las=1)
		}
		if(i==6){
			#6
			plot(trend$N0[samp], trend$Q[samp], pch=".", col='grey80',xlab=expression(N[0]), ylab='Q', las=1)
			kde.plot(trend$N0, trend$Q, n=100)
		}
		if(i==7){
			#7
			plot.new()
			box()
			cor.lab(trend$U, trend$Q, xscale=0.5, yscale=0.5, cex=1.4)
			fig.lab("U, Q", xscale=0.15, yscale=0.1, cex=1.1)
		}
		if(i==8){
			#8
			plot.new()
			box()
			cor.lab(trend$N0, trend$Q, xscale=0.5, yscale=0.5, cex=1.4)
			fig.lab(expression(paste(N[0],", Q")), xscale=0.15, yscale=0.1, cex=1.1)
		}
		if(i==9){
			#9
			plot(d.Q$x, d.Q$y, type='l', xlab='Q', ylab='Rel. Density', ylim=c(0,1), yaxt='n', las=1)
			axis(2, at=pretty(c(0,1)), las=1)
		}
		if(i==graph) legend("topright", legend=c("0.05", "0.5"), lwd=c(1,3), bty='n')
	}
	dev.off()
}
mvnorm <- "
	functions {
		matrix cov_matrix_2d(vector sigma, real rho) {
		    matrix[2,2] Sigma;
		    Sigma[1,1] = square(sigma[1]);
		    Sigma[2,2] = square(sigma[2]);
		    Sigma[1,2] = sigma[1] * sigma[2] * rho;
		    Sigma[2,1] = Sigma[1,2];
		    return Sigma;
		  }
	}
	data {
		int<lower=1> N; //number of obs
		vector[2] x[N]; //lengths and mortality
		int<lower=1> nyear; //number of years
		int<lower=1> year[N]; //year pointer
		vector[nyear] rtl; // realized take
	}
	parameters {
		real<lower=-1, upper=1> rho; //correlation
		vector<lower=0>[2] sigma; //sigma of mus
		real beta0; //int for rtl
		real beta1; //slope for rtl
		real mu0; //mu for m
	}
	transformed parameters{
		vector[2] mu[nyear];
		for(y in 1:nyear){
			mu[y,1] = beta0 + beta1*rtl[y];
			mu[y,2] = mu0;
		}
	}
	model {
		mu0 ~ normal(0,2);
		beta0 ~ normal(0,2);
		beta1 ~ normal(0,2);
		sigma ~ normal(0,2);
		(rho + 1) / 2 ~ beta(2, 2);

		for(n in 1:N){
			x[n] ~ multi_normal(mu[year[n]], cov_matrix_2d(sigma, rho));
		}
	}"
mv.DC.init <- function(chain_id){
	lm <- lm(log(Len)~Year, data=DC.td.df)

	beta0 <- rnorm(1, coef(lm)[1], abs(coef(lm)[1]*0.1))
	beta1 <- rnorm(1, coef(lm)[2], abs(coef(lm)[2]*0.1))

	mu0 <- rnorm(1, 0, 1)

	sigma <- rtruncnorm(2, a=0, mean = apply(DC.td.dat$x, 2, sd), sd = apply(DC.td.dat$x, 2, sd)*0.1)
	rho <- rtruncnorm(1, a=-1, b=1, cor(DC.td.dat$x)[1,2], abs(cor(DC.td.dat$x)[1,2]*0.1))

	return(list(beta0 = beta0,
	            beta1 = beta1,
	            mu0 = mu0,
	            sigma = sigma,
	            rho = rho))
}
mv.CC.init <- function(chain_id){
	lm <- lm(log(Len)~Year, data=CC.td.df)

	beta0 <- rnorm(1, coef(lm)[1], abs(coef(lm)[1]*0.1))
	beta1 <- rnorm(1, coef(lm)[2], abs(coef(lm)[2]*0.1))

	mu0 <- rnorm(1, 0, 1)

	sigma <- rtruncnorm(2, a=0, mean = apply(CC.td.dat$x, 2, sd), sd = apply(CC.td.dat$x, 2, sd)*0.1)
	rho <- rtruncnorm(1, a=-1, b=1, cor(CC.td.dat$x)[1,2], abs(cor(CC.td.dat$x)[1,2]*0.1))

	return(list(beta0 = beta0,
	            beta1 = beta1,
	            mu0 = mu0,
	            sigma = sigma,
	            rho = rho))
}
vbgm_stan_l0_ln <- "
	data{
		int<lower=1> n_obs; //number of observations
		vector<lower=0>[n_obs] age; //ages of fish
		vector<lower=0>[n_obs] l; //length of fish
		int<lower=1> nseq;
		vector<lower=0>[nseq] seq_ages;
		int<lower=1> n_hatchlings;
		vector<lower=0>[n_hatchlings] hatchling_size;
		real<lower=0> nester_size;
		real<lower=0> nester_sd;
	}
	parameters{
		real<lower=0> Linf; //L infinity
		real<lower=0> K; //vb K
		real<lower=0> Lknot; //t knot

		real<lower=0> sigma_obs; //sigma_obs
		real<lower=0> sigma_Lknot; //sigma_obs
	}
	model{
		vector[n_obs] Lpred;
		Linf ~ normal(nester_size, nester_sd);
		K ~ normal(0.1, 0.5);
		Lknot ~ normal(4, 0.2);

		sigma_obs ~ normal(0, 1);
		sigma_Lknot ~ normal(0, 1);

		for(i in 1:n_obs){
			Lpred[i] = Linf - (Linf-Lknot)*exp(-K*age[i]);
			 target += lognormal_lpdf(l[i]|log(Lpred[i]), sigma_obs);
		}
		for(i in 1:n_hatchlings){
			target += lognormal_lpdf(hatchling_size[i]|log(Lknot),sigma_Lknot);
		}
	}
	generated quantities{
		real Amat;
		real tknot;
		vector[nseq] Lage;
		Amat = (1/K)*log((Linf-Lknot)/(Linf*(1-0.975)));
		tknot = (1/K)*log((Linf-Lknot)/Linf);

		for(i in 1:nseq){
			Lage[i] = Linf - (Linf-Lknot)*exp(-K*seq_ages[i]);
		}
}"
init.vbgm.l0 <- function(chain_id){
	return(list(Linf = rtruncnorm(1, a=0, mean=nester_size, sd=nester_sd),
	       		K = rtruncnorm(1, a=0, mean=0.1, sd=0.01),
	       		Lknot = rtruncnorm(1,a=0, mean=hatchling_size,sd=hatchling_sd),
	       		sigma_obs = rtruncnorm(1, a=0, mean=1, sd=0.01),
				sigma_Lknot = rtruncnorm(1, a=0, mean=1, sd=0.01)))
}
