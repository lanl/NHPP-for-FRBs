# © (or copyright) 2016. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001 
# for Los Alamos National Laboratory (LANL), which is operated by Triad 
# National Security, LLC for the U.S. Department of Energy/National Nuclear 
# Security Administration.
#
# All rights in the program are reserved by Triad National Security, LLC, 
# and the U.S. Department of Energy/National Nuclear Security Administration. 
# The Government is granted for itself and others acting on its behalf a 
# nonexclusive, paid-up, irrevocable worldwide license in this material to 
# reproduce, prepare derivative works, distribute copies to the public, 
# perform publicly and display publicly, and to permit others to do so.
# 
# Additionally, redistribution and use in source and binary forms, with or 
# without modification, are permitted provided that the following conditions 
# are met:
# 1. Redistributions of source code must retain the above copyright 
# notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
# 3. Neither the name of Triad National Security, LLC, Los Alamos 
# National Laboratory, LANL, the U.S. Government, nor the names of its 
# contributors may be used to endorse or promote products derived from this 
# software without specific prior written permission.
 
# THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS 
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Estimate the rate of FRBs.
rm(list=ls())

library(RColorBrewer)
library(MASS)

# Scott's code to fit the model
source('likelihood.R')

# Data. Chrono ordered.
surveys <- c('Lorimer 2007', 'Deneva 2009', 'Keane 2010', 'Siemion 2011', 'Burgay 2012', 
             'Petroff 2014', 'Spitler 2014', 'Burke-Spolaor 2014',
             'Ravi 2015', 'Petroff 2015', 'Law 2015', 'Champion 2016')
n <- c(1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 10)
# Sensitivity at FWHM divided by 2
S <- c(0.590, 0.365, 0.447, 118, 0.529, 0.868, 0.210, 0.615, 0.555, 0.555, 0.240, 0.560) / 2
# FWHM diameter in arcminutes divided by 2 to get radius divide by 60 to get degrees
R <- c(14, 4, 14, 150, 14, 14, 4, 14, 14, 14, 60, 14)/(2*60)
# Number of beams.
beams <- c(13, 7, 13, 30, 13, 13, 7, 13, 13, 13, 1, 13)
# Time per beam
tpb <- c(490.5, 459.1, 1557.5, 135.8, 532.6, 926, 11503.3, 917.3, 68.4, 85.5, 166, 2786.5)
time <- tpb*beams
# observed flux
flux <- c(30, # Lorimer 2007
          0.4, 	# Spitler 2014
          0.3,	# Burke-Spolaor 2014
          2.1,	# Ravi et al. (estimated from Figure 3, 1432 MHz)
          0.47, # Petroff 2015
          1.3, 0.4, 0.5, 0.5, # Chamption 2016 (these values are actually from Thornton 2013)
          0.87, 0.42, 0.69, 1.67, 0.18, # Champion 2016 (estimated from Figure 1)
          2.2)  # The last Champion entry that Scott found somewhere (FRB Cat?)
          
llClosure <- makeNHPPlikelihoodFunction(n, S, time, R, flux)
fit <- optim(par=c(15,2.5), fn=llClosure$loglikelihood, gr=llClosure$gradient, 
			lower=c(1e-8, 0), method='L-BFGS-B', control=list(fnscale=-1))
			
# Convert to the regular parameters
fit$b <- fit$par[2] - 1
fit$a <- fit$par[1] / fit$b

# Convert to per sky per day
fit$a <- fit$a * 24 * 41253

## get confidence intervals on cumulative rates
s0 <- exp(seq(log(.1),log(1000), length=51))
conf <- llClosure$cumulativeRate(fit$par[1], fit$par[2], s0)
conf <- conf * 24 * 41253

# Replace with parametric bootstrap
## Parametric bootstrap CI
B <- 5000
V <- -solve(llClosure$hessian(fit$par))
parSample <- mvrnorm(n=B, mu=fit$par, Sigma=V)
# Some constraints here
i = which((parSample[,1] < 0) | (parSample[,2] < 1))
while(length(i) > 0) {
	parSample[i,] <- mvrnorm(n=length(i), mu=fit$par, Sigma=V)
	i = which((parSample[,1] < 0) | (parSample[,2] < 1))
}
LambdaSample <- outer(s0, -(parSample[,2]-1), FUN='^')
LambdaSample <- apply(LambdaSample, 1, '*', parSample[,1]/(parSample[,2]-1))
LambdaSample <- LambdaSample * 24 * 41253
conf[,2:3] <- t(apply(LambdaSample, 2, quantile, c(.025, .975)))

# Get the fit assuming a Euclidean scaling.
llEuclid <- function(alpha) { llClosure$loglikelihood( c(alpha, 2.5) ) }
fitEuclid <- optimize(f=llEuclid, interval=c(0,10), maximum=TRUE)
fitEuclid$a <- fitEuclid$maximum / 1.5
fitEuclid$a <- fitEuclid$a * 24 * 41253

# Plot the line, with uncertainty, and a line with -3/2 slope.
pdf(file='fit.pdf', height=6, width=6)
plot(s0, conf[,1], log='xy', type='l',
     xlim = range(s0), ylim=c(1e-1, 1e8),
     xlab='S=flux sensitivity limit (Jy)', 
     ylab=expression(plain('rate with flux > S ')*(sky^{-1}*day^{-1})))
matlines(s0, conf[,2:3], lty=2, col=1)
grid()
abline(log10(fitEuclid$a), -1.5, col='grey')

# Need length(surveys) colors
# 14 right now
# Start with Set1
col.survey <- c(brewer.pal(9, 'Set1'))
# Get rid of the yellow
col.survey <- col.survey[-6]
# Add some from Set2
col.survey <- c(brewer.pal(6, 'Dark2'), col.survey)
pch.survey <- 1:length(surveys)
denom <- time*pi*R^2
denom <- denom / ( 24 * 41253)
points(2*S, n/denom, col=col.survey, pch=pch.survey)
segments(2*S, 0.5*qchisq(0.025, 2*n)/denom, 2*S, 0.5*qchisq(0.975, 2*n+2)/denom, col=col.survey)
i <- which(n==0)
if (length(i)){ 
	yi <- 0.5*qchisq(0.975, 2*n[i]+2)/denom[i]
	arrows(2*S[i], yi/100, 2*S[i], yi, code=1, cex=.1, col=col.survey[i])
}
legend(x='topleft', legend=surveys, col=col.survey, pch=pch.survey, lty=1, ncol=3, xjust=1, cex=.8)
dev.off()

## Form Survivor function of flux given N and use it for a 
## log-log probability plot of the flux values.
## assume S is sorted
Pfun <- function(ss,b, S, time, area){
	# survivor function for observed flux
	ctas <- rev(c(0, cumsum(rev(time*area*S^-b))))
	cta <- c(0, cumsum(time*area))
	#browser()
	i <- 1+findInterval(ss, S)
	out <- (ctas[i] + cta[i]*ss^(-b)) / max(ctas)
	return(out)
}

Qfun <- function(p, b, S, time, area){
	# quantile function corresponding to Pfun
	ctas <- rev(c(0, cumsum(rev(time*area*S^-b))))
	cta <- c(0, cumsum(time*area))
	#browser()
	i <- 1+findInterval(-p, -Pfun(S, b, S, time, area))
	#browser()
	out <- ((p*max(ctas) - ctas[i]) / cta[i])^(-1/b)
	return(out)
}

## bootstrapped PIs
pe <- function(p,u){
	# empirical quantiles of u vector with elements in (0,1)
	N <- length(u)
	if (N==0) return(p)
	u <- c(0, sort(u), 1)
	ppoints <- c(0,(seq(N)-0.5)/N, 1)
	i <- findInterval(p, ppoints)
	out <- u[i] + (u[i+1]-u[i])*(p-ppoints[i])/(ppoints[i+1]-ppoints[i])
}

area <- pi*R^2
bhat <- fit$par[2]-1
p <- exp(seq(log(.01),log(1), length=101))
B <- 5000 # number of parametric bootstrap samples
Squants <- matrix(0, length(p), B) # to hold bootstrap output
SquantsEuclid <- matrix(0, length(p), B) # to hold bootstrap output
sdb <- sqrt(-solve(llClosure$hessian(fit$par))[2,2])
oo <- order(S)
for (i in seq(B)){
	bi <- -1
	while(bi < 0 ) {
		bi <- rnorm(1, bhat, sdb)
	}
	N <- rpois(1, sum(n))
	u <- runif(N)
	#s <- Qfun(u, ppoints(N), bi, S, time, area) # bootstrapped fluxes
	Squants[,i] <- Qfun(1-pe(1-p,u), bi, S[oo], time[oo], area[oo])
	SquantsEuclid[,i] <- Qfun(1-pe(1-p,u), 3/2, S[oo], time[oo], area[oo])
	#browser()
}
bounds <- t(apply(Squants, 1, quantile, p=c(0.025,0.975), na.rm=TRUE))
boundsEuclid <- t(apply(SquantsEuclid, 1, quantile, p=c(0.025,0.975), na.rm=TRUE))

# Make col.flux from col.survey
col.flux <- character(length(flux))
pch.flux <- numeric(length(flux))
com <- 1
for(i in 1:length(n)) {
	fin <- com+n[i]-1
	if(fin >= com) { 
		col.flux[seq(com,fin,1)] <- col.survey[i]
		pch.flux[seq(com,fin,1)] <- pch.survey[i]
		com <- fin+1
	}
}

flux.sort <- sort(flux, index.return=TRUE)

pdf(file='surv.pdf', height=6, width=6)
plot(flux.sort$x, rev(ppoints(flux)), log='xy', col=col.flux[flux.sort$ix], pch=pch.flux[flux.sort$ix],
	xlim=c(min(S), max(flux,S,1e3)), ylim=c(.01,1.05),
	xlab='S (Jy)', ylab='P{flux > S}')
lines(Qfun(p, bhat, S[oo], time[oo], area[oo]), p, col='black')
matlines(bounds, p, lty=2, col='black')
matlines(boundsEuclid, p, lty=2, col='grey')
rug(x=S)
grid()
legend(x='topright', legend=surveys[n>0], col=col.survey[n>0], pch=pch.survey[n>0], cex=.85)
dev.off()

# Generate a few results to compare to Champion et al and Rane et al.
sComp <- c(2/3, 4/3)
confComp95 <- llClosure$cumulativeRate(fit$par[1], fit$par[2], sComp)
confComp95 <- confComp95 * 24 * 41253

confComp99 <- llClosure$cumulativeRate(fit$par[1], fit$par[2], sComp, confLevel=.99)
confComp99 <- confComp99 * 24 * 41253
