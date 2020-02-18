library(deSolve)

#======================================================================
#Model without treatment: Run to verify the equilibrium
#   Parameter:
#
#		Y is the proportion of triatomines infected 
#		X is the proportion of dogs infected 
#		a is expected number of bites on dogs per triatomine (1/14)
#		m is equilibrium triatomine density per dog
#		n is length of incubation period in insects (45), 
#		b is the transmission efficiency from infectious triatomines to susceptible dog, through biting (0.00068)
#		c is the probability of infection of an uninfected triatomine by biting an infectious dog (0.49, or 0.28- halfed to account for cyclic parasitemia)
#   r is daily force of infected dog mortality 1/(365*3)=0.0009
#   g is the daily probability of bug mortality
#======================================================================

#=============================================

RM <- function(time, state, parameters) 
{
  with(
    as.list(c(state, parameters)), 
    {
      dX <- m*a*b*Y*(1-X)-r*X   
      dY <- a*c*X*(exp(-g*n)-Y)-(g*Y)
      dm <- (R*(1-m/K)*m )
      return(list(c(dX, dY, dm)))
    }
  )
}

#=============================================
#  Run the model and Plot
#=============================================

init <- c(X = 0.01, Y = 0, m= 40)
parameters <- c(a=1/14, b=0.00068, m=40, n=45, g=0.005, c=0.28, r=1/(3*365), K = 40, R= 0.09)
times <- seq(0, 20000, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = RM, parms = parameters))
RESULTS<-data.frame(out$X,out$Y)
RESULTSm <-data.frame(out$m)
#========Proportion dogs and bugs infected===
matplot(times, RESULTS, type = "l", xlab = "Time", ylab = "Proportion (X,Y)", main = "Pre-Tx Model", lwd = 1, lty = 1, bty = "l", col = c("black","red"))
legend(8000, .4, c("Dogs", "Triatomines"), pch = 1, col = c("black","red"))

#=============================================
#  Run the model
#=============================================
tail(RESULTS)
#       out.X     out.Y
#20001 0.5368288 0.5448028


#=============================================
#  II. TRY EQUILIBRIUM FUNCTION
#=============================================

getEq<-function(a, b, m, n, g, c, r)
{
  R0 <- (m*(a^2)*b*c*exp(-g*n)) / (r*g)
  S<-a/g
  Xeq<-(R0-1)/(R0+(c*S))
  Yeq<-((a*c*Xeq)/(g+a*c*Xeq))*exp((-g*n))
  togo<-c(R0,Xeq,Yeq)
  return(c(R0, Xeq, Yeq))
}

RMequil <- getEq(a=1/14, b=.00068, m=40, n=45, g=.005, c=.28, r=1/(3*365))

Xeq <- RMequil[2]
# 0.5368288
Yeq <- RMequil[3]
# 0.5448028


#==========================================================================
#  II. Model Two: Bugs and Dogs transmission dynamics after fluralaner tx
#       -Pull in the Equilibrium values (Xeq and Yeq) to use as starting points in treatment function (RMTx2)
#==========================================================================

#=============================================
#   estimate parameter z, percent of bugs that die after feeding on treated dogs, from data in Laino et al. 2019 https://www.ncbi.nlm.nih.gov/pubmed/30981313
#	 	x = days post fluralaner treatment
#		y = percentage of bugs that died after feeding on treated dogs (5th stage pyrethroid resistant instars) NOT USED
#		z = percentage of bugs that died after feeding on treated dogs (5th stage suceptible instars) USED 
#=============================================


#Initial vectors for days post treatment, % killed
x <- c(4, 30, 60, 90, 120, 210, 360) #days post treatment
y <- c(0.99, 1.0, 1.0, 0.47, 0.49, 0, 0)  #% killed
z <- c(1.0, 0.99, 0.99, 0.79, 0.7, 0.02, 0) #% killed

plot(z ~ x)

#=============================================
#  fit data with logistic curve
#       -using % dead 5th stage suceptible bugs against days post fluralaner tx
# 		-fit2 stores fitted values
#		-plot to check fit
#		-extract fit values using equation y = Asym / (1 + exp((xmid - input) / scal))
#=============================================
fit2 <- nls(z ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, z))
summary(fit2)

lines(seq(0, 400, length.out = 400),
      predict(fit2, newdata = data.frame(x = seq(0.5, 400, length.out = 400))))

Asym<-summary(fit2)$parameters[1,1]
xmid<-summary(fit2)$parameters[2,1]
scal<-summary(fit2)$parameters[3,1]


#=============================================
#     Parameters for model with treatment
#
#		Y - the proportion of triatomines infected 
#		X - the proportion of dogs infected 
#		a is expected number of bites on dogs per triatomine (1/14)
#		m is equilibrium triatomine density per dog
#		n is length of incubation period in insects (45), 
#		g is the daily force of triatomine mortality, mzl at K, the carrying capacity  
#		b is the transmission efficiency from infectious triatomines to susceptible dog, through biting (0.00068)
#		c is the probability of infection of an uninfected triatomine by biting an infectious dog
#   k is the probabilty of tranmission through oral ingestion of vectors (0.1)
#   r is daily force of infected dog mortality:   1.38/1460= 0.0009
#   p is probability dog eats the bug (to be varied)
#   R is maximum birthrate (estimate from T. brasiliensis)
#   K is carrying capacity of vectors per dog
#=============================================


#=============================================
#   force the time dependent covariate 
#       z is the probability of death from Fluralaner per bite
#       jr added code to apply more than one treatment
#       jr added function for implementing treatments at different times
#=============================================

times <- seq(0, 20000, by = 1)
signal <- data.frame(times = times, import = rep(0, length(times)))


trt2.function <- function(fit2, trt.days, signal) {
  Asym<-summary(fit2)$parameters[1,1]
  xmid<-summary(fit2)$parameters[2,1]
  scal<-summary(fit2)$parameters[3,1]
  trt.start <- trt.days
  trt.end <- c(trt.start + 400)
  trt.segments <- unlist(Map(':', trt.start, trt.end))
  time.segments <- rep(0:401, times=length(trt.start)) 
  signal$import[trt.segments] = (Asym / (1 + exp((xmid - times[time.segments]) / scal)))  
  return(signal)
}

#Input number of days and times to implement treament
trt.days <- c(10000) 
out <-trt2.function(fit2, trt.days, signal)


input <- approxfun(out, rule = 2)


#=============================================
#   model importing z from the above
#.    -added density dependent mortality g*(1-m/K) for bugs
#     -Code without density dependent death rate and birthrate for dogs 
#         
#=============================================

RMTx2 <- function(times, stateTx2, parametersTx2)   
{
  with(
    as.list(c(stateTx2, parametersTx2)), 
    {
      z <- input(times)
      dX <- ((m*a*b*Y)+(p*k*(a*m*z*Y)))*(1-X)-r*X #jr added density dependent mortality   #dx: change in proportion of infected dogs.  Vectorial transmission is the first term (m*a*b*Y), oral transmission the second (k*q*Y)
      dY <- a*c*X*(exp(-g*n)-Y)-((g*(1-m/K)*Y)+(m*a*z*Y)) #dx: change in proportion of infected bugs.  mzl added density dependent mortality
      dm <- ((R*(1-m/K)*m )+(-m*a*z)) #where the ratio of bugs: dog decreases by (-maz) and increases by logistic growth for bug population (R*(1-m/K)*m )
      return(list(c(dX, dY, dm)))
    }
  )
}

#=============================================
#  Run the model
#   -proportion of infected dogs (X) and bugs (Y) are taken from the proportions at the equilibrium of the model pre-treatment
#   -in 
#=============================================

initTx2 <- c(X = Xeq, Y= Yeq, m=40) 
parametersTx2 <- c(a=1/14, b=0.00068, n=45, g=0.005, c=0.28, k=.10, r= 1/(3*365), p=0.8, K=40, R= 0.09)
outTx2 <- as.data.frame(ode(y = initTx2, times = times, func = RMTx2, parms = parametersTx2))
RESULTS2 <- data.frame(outTx2$X, outTx2$Y)
#times <- seq(0, 10000, by = 1)
#timesTx2 <- seq(0, 10000, by = 1)


#=============================================
#   Plot Results
#   m     Bugs per Dog
#   X     Dog Prevalence
#   Y*m   Positive Bugs Per Dog
#   Y     Bug Prevalence
#=============================================
par(mfrow=c(1,2))

matplot(times, RESULTS2, type = "l", xlab = "Time", ylab = "Proportion (X,Y)", main = "Tx Model Using Equilibrium", lwd = 1, lty = 1, bty = "l", col = c("black","red"))
legend(8000, .4, c("Dogs", "Triatomines"), pch = 1, col = c("black","red"))


plot(outTx2$m, type='l',main='Bugs per Dog', ylab = "# bugs", xlab = "Days", cex.main=0.8)
abline(v=c(9995, 10120, 10240), col="green", lty=2)
plot(outTx2$Y*outTx2$m, type='l', main = 'Positive Bugs Per Dog', ylab = "# bugs", xlab = "Days", cex.main=0.8)
abline(v=9995, col="green", lty=2)
plot(outTx2$X, type='l', main= 'Dog Prevalence, 1 tx per year per 3 years', ylab = "Proportion infected", xlab = "Days", cex.main=0.8, xlim=c(0,20000))
abline(v=c(10000, 10400, 10800), col="green", lty=2)
plot(outTx2$Y, type='l', main = 'Bug prevalence, 1 tx per year for 3 years ', ylab = "Proportion infected", xlab = "Days", cex.main=0.8, xlim = c(0, 19000))
abline(v=c(10000, 10400, 10800), col="green", lty=2)

