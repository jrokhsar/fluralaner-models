#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$distPlot1 <- renderPlot({
        
        # Set parameter values
        #    ndays <- input$years
        #    nyrs  <- input$days
        #    doglifespan <-input$(1/(years*365)) 
        nyrs <- input$years
        r <- input$r        # dog lifespan in years
        m <- input$m        # ratio of bugs to dogs at beginning of simulation
        a <- input$a        # bite rate
        b <- input$b        # transmission efficiency from bugs to dogs
        c <- input$c        # transmission efficiency from dogs to bugs
        g <- input$g        # daily probability of bug mortality, not from treatment
        n <- 45             # incubation period of the T.cruzi in vector host
        R <- input$R         # bug birthrate
        K <- input$K        # bug carrying capacity
        p <- input$p        #percentage of dead bugs the dogs consume
        k <- input$k        # probabilty of transmission  bug to dog, oral, oral transmission 
        d <- input$d        # day to implement treatment
        variable <- input#variable
        
        
        
        ###### Set initial conditions, x= DPT
            #- y: 5th stage resistant instars
            #- z: 5th stage suceptible instars
        x <- c(4, 30, 60, 90, 120, 210, 360)
        y <- c(0.99, 1.0, 1.0, 0.47, 0.49, 0, 0)  #% killed
        z <- c(1.0, 0.99, 0.99, 0.79, 0.7, 0.02, 0)
        fit2 <- nls(z ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, z))
        Asym<-summary(fit2)$parameters[1,1]
        xmid<-summary(fit2)$parameters[2,1]
        scal<-summary(fit2)$parameters[3,1]
        
        #Create empty dataframe to store values of the curve for bugs killed vs DPT
        times2 <- seq(0, 20000, by = 1)
        signal <- data.frame(times = times2, import = rep(0, length(times2)))
        
        #Function for creating the dataframe of bugs killed vs days post treatment depending on user input, d
        trt2.function <- function(fit2, trt.days, signal) {
            Asym<-summary(fit2)$parameters[1,1]
            xmid<-summary(fit2)$parameters[2,1]
            scal<-summary(fit2)$parameters[3,1]
            trt.start <- trt.days
            trt.end <- c(trt.start + 400)
            trt.segments <- unlist(Map(':', trt.start, trt.end))
            time.segments <- rep(0:401, times2=length(trt.start)) 
            signal$import[trt.segments] = (Asym / (1 + exp((xmid - times2[time.segments]) / scal)))  
            return(signal)
        }
        
        #Input number of days and times to implement treament  ***WANT TO TURN THIS INTO A USER DEFINE VECTOR
    #    trt.days <- textInput(variable)
        trt.days <-
        #        trt.days <- c(d) 
        day <-trt2.function(fit2, trt.days, signal)
        
        
        input <- approxfun(day, rule = 2)
        
        
        #Initial conditions for ODE
        y0 <- c(X = 0.01, Y = 0, m=m)
        parameters <- list(a=a, b=b, c=c, m=m, n=n, g=g, r= 1/(r*365), R=R, K=K, p=p, k=k)
        times <- seq(from=0, 365*nyrs+1, by=1 )
        
        ###### Specify system of differential equations 
        SetODEs <- function(times, y, parameters){
            X <- y[1]
            Y <- y[2]
            m <- y[3]
            with(
                as.list(c(parameters)), 
                {
                    z <- input(times)
                    dX <- ((m*a*b*Y)+(p*k*(a*m*z*Y)))*(1-X)-r*X
                    dY <- a*c*X*(exp(-g*n)-Y)-((g*(1-m/K)*Y)+(m*a*z*Y)) 
                    dm <- ((R*(1-m/K)*m )+(-m*a*z))
                    return(list(c(dX, dY, dm)))
                }
            )
        }
        
        ###### Run the ode solver function
        out <- as.data.frame(ode(y = y0, times = times, func = SetODEs, parms = parameters))
        RESULTS<-data.frame(out$X,out$Y)
        
        
        ###### Plot the output
        matplot(times/365, RESULTS, type = "l", xlab = "Time", ylab = "Proportion (X,Y)", main = "Model with day to implement treatment", lwd = 1, lty = 1, bty = "l", col = c("black","red"))
        legend("topright", c("Dogs", "Triatomines"), pch = 1, col = c("black","red"))
        
    })
    output$distPlot2 <- renderPlot({
        
        # Set parameter values
        #    ndays <- input$years
        #    nyrs  <- input$days
        #    doglifespan <-input$(1/(years*365)) 
        nyrs <- input$years
        r <- input$r        # dog lifespan in years
        m <- input$m        # ratio of bugs to dogs at beginning of simulation
        a <- input$a        # bite rate
        b <- input$b        # transmission efficiency from bugs to dogs
        c <- input$c        # transmission efficiency from dogs to bugs
        g <- input$g        # daily probability of bug mortality, not from treatment
        n <- 45             # incubation period of the T.cruzi in vector host
        R <- input$R         # bug birthrate
        K <- input$K        # bug carrying capacity
        p <- input$p        #percentage of dead bugs the dogs consume
        k <- input$k        # probabilty of transmission  bug to dog, oral, oral transmission 
        d <- input$d        # day to implement treatment
        
        
        
        ###### Set initial conditions
        x <- c(4, 30, 60, 90, 120, 210, 360)
        z <- c(1.0, 0.99, 0.99, 0.79, 0.7, 0.02, 0)
        fit2 <- nls(z ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, z))
        Asym<-summary(fit2)$parameters[1,1]
        xmid<-summary(fit2)$parameters[2,1]
        scal<-summary(fit2)$parameters[3,1]
        
        #Create empty dataframe to store values of the curve for bugs killed vs DPT
        times2 <- seq(0, 20000, by = 1)
        signal <- data.frame(times = times2, import = rep(0, length(times2)))
        
        #Function for creating the dataframe of bugs killed vs days post treatment depending on user input, d
        trt2.function <- function(fit2, trt.days, signal) {
            Asym<-summary(fit2)$parameters[1,1]
            xmid<-summary(fit2)$parameters[2,1]
            scal<-summary(fit2)$parameters[3,1]
            trt.start <- trt.days
            trt.end <- c(trt.start + 400)
            trt.segments <- unlist(Map(':', trt.start, trt.end))
            time.segments <- rep(0:401, times2=length(trt.start)) 
            signal$import[trt.segments] = (Asym / (1 + exp((xmid - times2[time.segments]) / scal)))  
            return(signal)
        }
        
        #Input number of days and times to implement treament  ***WANT TO TURN THIS INTO A USER DEFINE VECTOR
        trt.days <- c(d) 
        day <-trt2.function(fit2, trt.days, signal)
        
        
        input <- approxfun(day, rule = 2)
        
        
        #Initial conditions for ODE
        y0 <- c(X = 0.01, Y = 0, m=m)
        parameters <- list(a=a, b=b, c=c, m=m, n=n, g=g, r= 1/(r*365), R=R, K=K, p=p, k=k)
        times <- seq(from=0, 365*nyrs+1, by=1 )
        
        ###### Specify system of differential equations 
        SetODEs <- function(times, y, parameters){
            X <- y[1]
            Y <- y[2]
            m <- y[3]
            with(
                as.list(c(parameters)), 
                {
                    z <- input(times)
                    dX <- ((m*a*b*Y)+(p*k*(a*m*z*Y)))*(1-X)-r*X
                    dY <- a*c*X*(exp(-g*n)-Y)-((g*(1-m/K)*Y)+(m*a*z*Y)) 
                    dm <- ((R*(1-m/K)*m )+(-m*a*z))
                    return(list(c(dX, dY, dm)))
                }
            )
        }
        
        ###### Run the ode solver function
        out <- as.data.frame(ode(y = y0, times = times, func = SetODEs, parms = parameters))
        RESULTS<-data.frame(out$X,out$Y)
        RESULTS2<- data.frame(out$m, out$m*out$Y)
        
        
        ###### Plot the output
        matplot(times/365, RESULTS2, type = "l", xlab = "Time", ylab = "# Bugs", main = "Model with tx, bugs:dogs", lwd = 1, lty = 1, bty = "l", col = c("black","red"))
        legend("topright", c("Bugs:dogs", "Infected bugs:dogs"), pch = 1, col = c("blue","orange"))
        
    })
    
    
    
    
}