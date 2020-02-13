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

# Chagas fluralaner model, prior to treatment 
server <- function(input, output) {
    
    output$distPlot <- renderPlot({
        
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

        
        ###### Set initial conditions
        y0 <- c(X = 0.01, Y = 0, m=m)
        parameters <- list(a=a, b=b, c=c, m=m, n=n, g=g, r= 1/(r*365), R=R, K=K)
        times <- seq(from=0, 365*nyrs+1, by=1 )
        
        ###### Specify system of differential equations 
        SetODEs <- function(times, y, parameters){
            X <- y[1]
            Y <- y[2]
            m <- y[3]
            with(
                as.list(c(parameters)), 
                {
                    dX <- m*a*b*Y*(1-X)-r*X   
                    dY <- a*c*X*(exp(-g*n)-Y)-(g*Y)
                    dm <- (R*(1-m/K)*m )
                    return(list(c(dX, dY, dm)))
                }
            )
        }
        
        ###### Run the ode solver function
        out <- as.data.frame(ode(y = y0, times = times, func = SetODEs, parms = parameters))
        RESULTS<-data.frame(out$X,out$Y)
        
        
        ###### Plot the output
        matplot(times/365, RESULTS, type = "l", xlab = "Time", ylab = "Proportion (X,Y)", main = "Pre-Tx Model", lwd = 1, lty = 1, bty = "l", col = c("black","red"))
        legend(8000, .4, c("Dogs", "Triatomines"), pch = 1, col = c("black","red"))
        
    })
    
    
    
}
