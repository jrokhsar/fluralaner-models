
#==========================================================================
#  II. Model Two: Bugs and Dogs transmission dynamics after fluralaner tx
#==========================================================================



#=============================================
#   estimate parameter z, percent of bugs that die after feeding on treated dogs, from data in Laino et al. 2019 https://www.ncbi.nlm.nih.gov/pubmed/30981313
#	 	x = days post fluralaner treatment
#		y = percentage of bugs that died after feeding on treated dogs (5th stage pyrethroid resistant instars) NOT USED
#		z = percentage of bugs that died after feeding on treated dogs (5th stage suceptible instars) USED 
#   For x, y, and z, tried to implement the treatment after 365 days
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
#   force the time dependent covariate 
#       z is the probability of death from Fluralaner per bite
#=============================================

#=============================================
#   TRY TO MAKE THE FUNCTION FOR INPUTTING DIFFERENT TREATMENT TIMES
#=============================================

times <- seq(0, 20000, by = 1)
signal <- data.frame(times = times, import = rep(0, length(times)))


######TEST DIFFERENT INPUTS FOR FUNCTION
# inputs
Asym<-summary(fit2)$parameters[1,1]
xmid<-summary(fit2)$parameters[2,1]
scal<-summary(fit2)$parameters[3,1]
trt.start <- c(10000,10090,10180,10270)  #input for start times of treatments 
trt.end <- c(trt.start + 400) #input for the end times of treatments
trt.segments <- unlist(Map(':', trt.start, trt.end)) #Create the treatment segments
time.segments <- rep(0:401, times=length(trt.start)) #create the corresponding time segments 

signal$import[trt.segments] = (Asym / (1 + exp((xmid - times[time.segments]) / scal)))  
 

#Check that the lengths of the time segments and treatment segments are the same
length(signal$import[trt.segments])  
length(signal$times[time.segments])

#Plot to verify correct - this is correct
(plot(times,signal$import, xlim = c(9900, 14000),
      main="(Asym / (1 + exp((xmid - times)", cex.main=0.8, cex=0.5))


#=============================================
#   Put all the above information in the function below
#=============================================

times <- seq(0, 20000, by = 1)
signal <- data.frame(times = times, import = rep(0, length(times)))

trt.function <- function(model_fit, trt_days, signal) {
  #inputs
  Asym<-summary(fit2)$parameters[1,1]
  xmid<-summary(fit2)$parameters[2,1]
  scal<-summary(fit2)$parameters[3,1]
  trt.start = trt_days
  trt.end = c(trt.start + 400)
  trt.segments <- unlist(Map(':', trt.start, trt.end))
  time.segments <- rep(0:401, times=length(trt.start))
  signal$import[trt.segments] = (Asym / (1 + exp((xmid - times[time.segments] / scal))))
  return(signal)
}


#Test function

trt.days =c(10000,10090,10180,10270) #Vector of treatment days

out <-trt.function(fit2, trt.days, signal)



#Plot to verify correct - this is correct
(plot(times, out$import, xlim = c(9900, 14000),
      main="(Asym / (1 + exp((xmid - times)", cex.main=0.8, cex=0.5))


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

trt.days <- c(10000,10090,10180,10270) 

out <-trt2.function(fit2, trt.days, signal)

input <- approxfun(out, rule = 2)


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

initTx2 <- c(X = 0.01, Y= 0, m=40) 
parametersTx2 <- c(a=1/14, b=0.00068, n=45, g=0.01, c=0.28, k=.10, r= 0.0009, p=0.8, K=40, R= 0.09)
outTx2 <- as.data.frame(ode(y = initTx2, times = times, func = RMTx2, parms = parametersTx2))
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
plot(outTx2$m, type='l',main='Bugs per Dog', ylab = "# bugs", xlab = "Days", cex.main=0.8)
abline(v=c(9995, 10120, 10240), col="green", lty=2)
plot(outTx2$Y*outTx2$m, type='l', main = 'Positive Bugs Per Dog', ylab = "# bugs", xlab = "Days", cex.main=0.8)
abline(v=9995, col="green", lty=2)
plot(outTx2$X, type='l', main= 'Dog Prevalence, 1 tx per year per 3 years', ylab = "Proportion infected", xlab = "Days", cex.main=0.8, xlim=c(9500,20000))
abline(v=c(10000, 10400, 10800), col="green", lty=2)
plot(outTx2$Y, type='l', main = 'Bug prevalence, 1 tx per year for 3 years ', ylab = "Proportion infected", xlab = "Days", cex.main=0.8, xlim = c(9900, 10400))
abline(v=c(10000, 10400, 10800), col="green", lty=2)


