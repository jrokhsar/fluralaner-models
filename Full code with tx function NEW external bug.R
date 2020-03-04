
library(deSolve)

#======================================================================
#  I. Model One: Bugs and Dogs transmission dynamics pre-fluralaner tx:
#       - jr added in dm to see how model 2 would look prior to treatment and determine when X and Y reach equilibrium
#       - X doesn't reach equilibrium until after 10000 days
#       - Model one and has been combined with model 2, just here for reference.
#======================================================================


#=============================================
# Transmission dynamics for dogs and bugs pre-fluralaner tx
#     Parameters
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
#  Run the model
#=============================================

init <- c(X = 0.01, Y = 0, m= 40)
parameters <- c(a=1/14, b=0.00068, m=10, n=45, g=0.005, c=0.28, r=1/(3*365), K = 40, R= 0.09)
times <- seq(0, 40000, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = RM, parms = parameters))
RESULTS<-data.frame(out$X,out$Y)
RESULTSm <-data.frame(out$m)



#=============================================
#  Plot results
#=============================================

#========Proportion dogs and bugs infected===
matplot(times, RESULTS, type = "l", xlab = "Time", ylab = "Proportion (X,Y)", main = "Pre-Tx Model", lwd = 1, lty = 1, bty = "l", col = c("black","red"))
legend(8000, .4, c("Dogs", "Triatomines"), pch = 1, col = c("black","red"))

#========Ratio of bugs:dogs===
matplot(times, RESULTSm, type = "l", xlab = "Time", ylab = "Ratio bugs:dogs", main = "Pre-Tx Model, Bugs:Dogs", lwd = 1, lty = 1, bty = "l", col = "green")
legend(8000, .4, c("Bugs:Dogs"), pch = 1, col = c("green"))


#==========================================================================
#  II. Model Two: Bugs and Dogs transmission dynamics after fluralaner tx
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
#     Parameters
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

#trt.days <- c(10000, 10365, 10730, 11095)
#trt.days <- c(10000, 10090, 10180, 10270, 10365)
out <-trt2.function(fit2, trt.days, signal)


input <- approxfun(out, rule = 2)


#=============================================
#   model importing z from the above
#      -jr added in an external source of bugs not affected by tx MM
#         
#=============================================

RMTx2 <- function(times, stateTx2, parametersTx2)   
{
  with(
    as.list(c(stateTx2, parametersTx2)), 
    {
      z <- input(times)
      dX <- (((m*a*b*Y)+MM)+(p*k*(a*m*z*Y)))*(1-X)-r*X 
      #dX <- ((m*a*b*Y)+(p*k*(a*m*z*Y)))*(1-X)-r*X 
      dY <- a*c*X*(exp(-g*n)-Y)-((g*Y)+(m*a*z*Y))
      dm <- (R*(1-m/K)*m )+(-m*a*z)
      return(list(c(dX, dY, dm)))
    }
  )
}

#=============================================
#  Run the model
#   -proportion of infected dogs (X) and bugs (Y) are taken from the proportions at the equilibrium of the model pre-treatment
#   -in 
#=============================================

initTx2 <- c(X = 0.01, Y= 0, m=40) 
parametersTx2 <- c(a=1/14, b=0.00068, n=45, g= 0.005, c=0.28, k= 0.10, r= 1/(3*365), p=0.8, K=40, R= 0.09, MM=0.004)
outTx2 <- as.data.frame(ode(y = initTx2, times = times, func = RMTx2, parms = parametersTx2))
RESULTS2<-data.frame(outTx2$X,outTx2$Y)
RESULTS2m <-data.frame(outTx2$m, outTx2$Y*outTx2$Y)
#timesTx2 <- seq(0, 10000, by = 1)


#=============================================
#   Plot Results
#   m     Bugs per Dog
#   X     Dog Prevalence
#   Y*m   Positive Bugs Per Dog
#   Y     Bug Prevalence
#=============================================
par(mfrow=c(1,3))
plot(outTx2$m, type='l',main='Bugs per Dog - dm', ylab = "# bugs", xlab = "Days", cex.main=0.8)

plot(outTx2$Y*outTx2$m, type='l', main = 'Positive Bugs Per Dog', ylab = "# bugs", xlab = "Days", cex.main=0.8)
abline(v=c(10000, 10084, 10168, 10252), col="green", lty=2)
plot(times/365, outTx2$X, type='l', main= 'Dog Prevalence, Treating q 12 weeks for One Year', ylab = "Proportion infected", xlab = "Years", 
     cex.main=1, col= "blue", xlim= c(9990/365, 12500/365))
abline(v=c(10000/365, 10084/365, 10168/365, 10252/365), col="green", lty=2, lwd=2)
plot(times/365, outTx2$Y, type='l', main = 'Bug Prevalence, Treating q 12 weeks for One Year', ylab = "Proportion infected", xlab = "Days", 
     cex.main=1, col= "red")
abline(v=c(10000/365, 10084/365, 10168/365, 10252/365), col="green", lty=2, lwd=2)


#=============================================
#  Plot results with bugs and dogs together on one graph
#  Bugs per dog and positive bugs per dogs on the second graph
#=============================================

matplot(times/365, RESULTS2, type = "l", xlab = "Time (Years)", ylab = "Proportion Infected (%)", main= "External Bug Constant set to 0.004",
        lwd = 3, lty = 1, bty = "l", col = c("blue","red"))
       # xlim=c(9900/365, 13300/365))
legend(12300/365, 0.1, c("Dogs", "Bugs"), pch = 16, col = c("blue","red"), bty = "l")
#abline(v=c(10000/365), col="green", lty=2, lwd=3)


matplot(times/365, RESULTS2, type = "l", xlab = "Time (Years)", ylab = "Proportion Infected (%)", main= "One Treatment/Year, 4 Years", lwd = 3, lty = 1, bty = "l", col = c("blue","red"),
        xlim=c(9900/365, 14000/365))
legend(12500/365, 0.1, c("Dogs", "Bugs", "Treatment"), pch = 16, col = c("blue","red", "green"), bty = "n")
abline(v=c(10000/365, 10365/365, 10730/365, 11095/365), col="green", lty=2, lwd=3)



#========Proportion dogs and bugs infected===
matplot(times/365, RESULTS2, type = "l", xlab = "Time (Years)", ylab = "Proportion Infected (%)", main= "Treatment Every 90 Days per 1 Year", lwd = 3, lty = 1, bty = "l", col = c("blue","red"),
    xlim=c(9900/365, 14000/365))
legend(12500/365, 0.1, c("Dogs", "Bugs", "Treatment"), pch = 16, col = c("blue","red", "green"), bty = "n")
abline(v=c(10000/365, 10090/365, 10180/365, 10270/365), col="green", lty=2, lwd=3)
#========Ratio of bugs:dogs===
matplot(times, RESULTS2m, type = "l", xlab = "Time", ylab = "Ratio bugs:dogs", main = "Pre-Tx Model, Bugs:Dogs", lwd = 1, lty = 1, bty = "l", col = "green")
legend(8000, .4, c("Bugs:Dogs"), pch = 1, col = c("green"))
