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

### Loglikelihood of power-law NHPP associated with single-dish
### telescope observations of burst events.  Input data for each 
### of several surveys consists of 
### -- time on sky
### -- minimum detectable flux
### -- telescope effective area function
### The effective area function for a telescope gives the b-scaled effective
### area defined as 
### 	A(b) = 2*pi* integral_0^infty r*(p(r))^b dr
### where p(r) is the relative sensitivity at radius r (degrees) with p(0)=1.
### The effective area function must return a matrix with rows corresponding 
### to surveys and three columns, A=A(b), Bd=A'/A and Bdd=A"/A.
### 

makeNHPPlikelihoodFunction <- function(n, sensitivity, time, radius, flux){
  ## Vector 'n' has the number of observations in each study.
  ## Vector 'sensitivity' has minimum detectible fluxes for each study.
  ## Vector 'time' has time on sky for each study.
  ## Vector 'radius' has half-max radii for each study.
  ## Vector 'flux' of length sum(n) has flux values for all detections.
  ##
  ## This version assumes Gaussian primary beams.
  ##
  	N <- sum(n)
  	M <- length(n)
  	if (length(flux)!=N) 
  		stop(paste('length of flux (', length(flux), ') must equal sum of n (', N, sep=''))
	if (length(sensitivity)!=M)
		stop(paste('length of sensitivity (', length(sensitivity), ') must equal lenth of  n (', M, sep=''))
	if (length(time)!=M)
		stop(paste('length of time (', length(time), ') must equal lenth of  n (', M, sep=''))
 	if (length(radius)!=M) 
		stop(paste('length of radius (', length(radius), ') must equal lenth of  n (', M, sep=''))
 
	## Gaussian effective area function
	area <- function(radius, b){
		A <- pi*radius^2/(b*log(2))
		Bdot <- -1/b
		Bdotdot <- 2/b^2
		out <- cbind(A=A, Bd=Bdot, Bdd=Bdotdot)
		return(out)
	}
		
	## power integral function 
	powerIntegral <- function(sensitivity, beta){
		II <- sensitivity^-(beta-1) /(beta-1)
		IIdot <- -II*(1/(beta-1) + log(sensitivity))
		IIdotdot <- -2*IIdot/(beta-1) + II*log(sensitivity)^2
		out <- cbind(I=II, Id=IIdot, Idd=IIdotdot)
		return(out)
	}

	## loglikelihood function
	loglikelihoodList <- function(alpha, beta){
		## Return a list with NHPP loglikelihood, gradient and hessian.
		
		## values related to area
		aa <- area(radius, beta-1)
		A <- aa[,'A']
		Bd <- aa[,'Bd']
		Bdd <- aa[,'Bdd']
		taa <- time * A * alpha
		## power integrals
		pow <- powerIntegral(sensitivity, beta)
		I <- pow[,'I']
		Id <- pow[,'Id']
		Idd <- pow[,'Idd']

		## loglikelihood
		ll <- sum(-taa*I + n*log(taa)) - beta*sum(log(flux))
		## gradient
		g1 <- sum(-taa*I/alpha + n/alpha)
		g2 <- sum(-taa*(Bd*I+Id) + n*Bd) - sum(log(flux))
		grad <- c(g1, g2)
		## hessian
		h11 <- sum(-n/alpha^2)
		h12 <- sum(-taa*(Bd*I+Id)/alpha)
		h22 <- sum(-taa*(Bdd*I+2*Bd*Id+Idd) + n*(Bdd-Bd^2))
		hess <- matrix(c(h11, h12, h12, h22), 2, 2)
		
		## return 
		out <- list(loglikelihood=ll, gradient=grad, hessian=hess)
		return(out)
	}
	
	loglikelihood <- function(par) loglikelihoodList(par[1], par[2])$loglikelihood
	gradient <- function(par) loglikelihoodList(par[1], par[2])$gradient
	hessian <- function(par) loglikelihoodList(par[1], par[2])$hessian
	
	## confidence bands for Lambda(s)
	cumulativeRate <- function(alpha, beta, s, confLevel=0.95){
		## For scalar power paramteres alpha, beta, and vector
		## s of flux values return a three column matrix with 
		## rows for s and columns Lhat=Lambda(s), Llo=lower bound, 
		## Lhi=upper bound, using Wald-type variance estimates
		## when input alpha, beta are the MLEs
		
		lll <- loglikelihoodList(alpha, beta)
		L <- lll$loglikelihood
		Ld <- lll$gradient
		Ldd <- lll$hessian
		xs <- (1/(beta-1) + log(s))
		J <- Ka <- Kb <- H <- array(0, dim=c(2,2,length(s)))
		J[1,1,] <- alpha
		J[1,2,] <- alpha*xs
		J[2,2,] <- 1
		Ka[1,2,] <- alpha*xs
		Ka[2,1,] <- alpha*xs
		Ka[2,2,] <- alpha*xs*xs
		Kb[2,2,] <- 1
		
		sdl <- rep(0, length(s)) # std dev of log(Lambda(s))
		for (i in seq(length(s))){
			H[,,i] <- t(J[,,i])%*%Ldd%*%J[,,i] + Ld[1]*Ka[,,i] + Ld[2]*Kb[,,i]
			sdl[i] <- sqrt(solve(-H[,,i])[1,1])
		}
		ell <- log(alpha/(beta-1)) -(beta-1)*log(s)
		z <- -qnorm((1-confLevel)/2)
		Lhat <- exp(ell)
		Llo <- exp(ell-z*sdl)
		Lhi <- exp(ell+z*sdl)
		out <- cbind(Lhat=Lhat, Llo=Llo, Lhi=Lhi)
		return(out)
	}
	
	## return likelihood-related functions and cumulative rate function
	return(list(loglikelihoodList=loglikelihoodList, 
				loglikelihood=loglikelihood,
				gradient=gradient,
				hessian=hessian,
				cumulativeRate=cumulativeRate))
  }


if (FALSE){
### Simulate detections and estimate cumulative rate function.

a <- 10 # rate per sky*day greater than 1 Jy 
b <- 3/2 # Euclidian flux scaling
alpha0 <- a*b
beta0 <- b+1
M <- 4 # number of surveys
S <- round(exp(seq(log(1), log(25), length=M))) # log-spaced sensitivities (Jy)
R <- rep(1, length=M) # equal beam radii (deg^2)
En <- rep(3, length=M) # expected numbers of detections
area <- pi*R^2/(b*log(2))
time <- En*S^b/(a*area) # time to reach EN detections (days)

## 1) example
n <- c(3,3,2,1)
fluxList <- list(
	c(2.238773, 2.271227, 5.082403),
	c(24.751881,  6.375682, 15.886933),
	c(13.20281, 18.58447),
	c(48.76736))

## 2) random data
#n <- rpois(M, En) # numbers of detections
#fluxList <- vector('list', length=M)
#for (i in seq(M)) fluxList[[i]] <- S[i]*runif(n[i])^(-1/b)

flux <- unlist(fluxList)

llClosure <- makeNHPPlikelihoodFunction(n, S, time, R, flux)
ll <- llClosure$loglikelihoodList
ll(alpha0,beta0)

## check gradient and hessian against numerical derivatives
# lval <- ll(alpha0, beta0)
# lval$loglikelihood
# lval$gradient / gradient(llClosure$loglikelihood, c(alpha0,beta0))
# lval$hessian / hessian(llClosure$loglikelihood, c(alpha0,beta0))
## Bingo!

## get MLEs
fit <- optim(par=c(alpha0,beta0), fn=llClosure$loglikelihood, gr=llClosure$gradient, 
			lower=0, method='L-BFGS-B', control=list(fnscale=-1))

## get confidence intervals on cumulative rates
s0 <- exp(seq(log(.1),log(1000), length=51))
conf <- llClosure$cumulativeRate(fit$par[1], fit$par[2], s0)

## plot estimated rates
#pdf('rates.pdf', width=6, height=5)
#par(mar=c(5,5,0,1)+.1)
plot(1,1, type='n', 	xlim=range(s0), ylim=c(.0005, 100), log='xy',
	xlab='S=flux sensitivity limit (Jy)', 
	ylab=expression(plain('rate of FRBs with flux greater than S ')*(day^{-1}*deg^{-2})))
abline(a=log10(a), b=-b, col='grey', lwd=4)
matlines(s0, conf, lty=c(1,2,2), col=1, 
	lwd=c(2,1,1), type='l')

## add intervals for each study, assuming known b
denom <- time*area
points(S, n/denom)
segments(S, 0.5*qchisq(0.025, 2*n)/denom, S, 0.5*qchisq(0.975, 2*n+2)/denom)
i <- which(n==0)
if (length(i)){ 
	yi <- 0.5*qchisq(0.975, 2*n[i]+2)/denom[i]
	arrows(S[i], yi/100, S[i], yi, code=1, cex=.1)
	points(S[i], 0.5*qchisq(0.5, 2*n[i]+2)/denom[i])
	}
#dev.off()

## Form Survivor function of flux given N and use it for a 
## log-log probability plot of the flux values.
## assume S is sorted
Pfun <- function(s,b, S, time, area){
	# survivor function for observed flux
	ctas <- rev(c(0, cumsum(rev(time*area*S^-b))))
	cta <- c(0, cumsum(time*area))
	i <- 1+findInterval(s, S)
	out <- (ctas[i] + cta[i]*s^(-b)) / max(ctas)
	return(out)
}

Qfun <- function(p, b, S, time, area){
	# quantile function corresponding to Pfun
	ctas <- rev(c(0, cumsum(rev(time*area*S^-b))))
	cta <- c(0, cumsum(time*area))
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

bhat <- fit$par[2]-1
p <- exp(seq(log(.01),log(1), length=101))
B <- 5000 # number of parametric bootstrap samples
Squants <- matrix(0, length(p), B) # to hold bootstrap output
SquantsEuclid <- matrix(0, length(p), B) # to hold bootstrap output
sdb <- sqrt(-solve(llClosure$hessian(fit$par))[2,2])
for (i in seq(B)){
	bi <- rnorm(1, bhat, sdb)
	N <- rpois(1, sum(n))
	u <- runif(N)
	#s <- Qfun(u, ppoints(N), bi, S, time, area) # bootstrapped fluxes
	Squants[,i] <- Qfun(1-pe(1-p,u), bi, S, time, area)
	SquantsEuclid[,i] <- Qfun(1-pe(1-p,u), 3/2, S, time, area)
	#browser()
}
bounds <- t(apply(Squants, 1, quantile, p=c(0.025,0.975), na.rm=TRUE))
boundsEuclid <- t(apply(SquantsEuclid, 1, quantile, p=c(0.025,0.975), na.rm=TRUE))

#pdf('pplot.pdf', width=6, height=6)
plot(sort(flux), rev(ppoints(flux)), log='xy', 
	xlim=c(min(S), max(flux,S,1e3)), ylim=c(.01,1),
	xlab='observed flux', ylab='probability of a greater flux')
lines(Qfun(p, bhat, S, time, area), p, col='black')
#lines(Qfun(p, b, S, time, area), p, col='red')
matlines(bounds, p, lty=2, col='black')
#matlines(boundsEuclid, p, lty=2, col='grey')
rug(x=S)
#dev.off()


}
