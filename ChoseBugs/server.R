#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

server <- function(input, output) {
    
    output$distPlot1 <- renderPlot({
        
        # Set parameter values
        nyrs <- input$years
        r <- input$r        # dog lifespan in years
        m <- input$m        # ratio of bugs to dogs- starting off at carrying capacity
        a <- input$a        # bite rate
        b <- input$b        # transmission efficiency from bugs to dogs
        c <- input$c        # transmission efficiency from dogs to bugs
        g <- input$g        # daily probability of bug mortality, not from treatment
        n <- 45             # incubation period of the T.cruzi in vector host
        R <- input$R        # bug birthrate
        p <- input$p        # percentage of dead bugs the dogs consume
        d <- input$d        # day to implement treatment
        MM <- input$MM      # external bug population, unaffected by treatment
        k <- input$k        # transmission efficiency bugs-dogs through ingestion 
        K <- input$m        # carrying capacity
        
        
        #=============================================
        #   estimate parameter z, percent of bugs that die after feeding on treated dogs, from data in Laino et al. 2019 https://www.ncbi.nlm.nih.gov/pubmed/30981313
        #	 	x = days post fluralaner treatment
        #		y = percentage of bugs that died after feeding on treated dogs (5th stage pyrethroid resistant instars) NOT USED
        #		z = percentage of bugs that died after feeding on treated dogs (5th stage suceptible instars) USED
        #       trt.function2: function to make dataframe of % bugs killed vs. DPT based on tx start day and vector of tx days 
        #       force time depedent covariate, z, into the model
        #=============================================
        
        
        x <- c(4, 30, 60, 90, 120, 210, 360)
        y <- c(0.99, 1.0, 1.0, 0.47, 0.49, 0, 0)  #% killed
        z <- c(1.0, 0.99, 0.99, 0.79, 0.7, 0.02, 0) #% killed
        
        #fit2 <- nls(z ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, z))
        
        #Option to choice susceptible or resistant commented out because not loading in the app
        if (input$select == "fifth stage resistant nymphs"){

            fit2 <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y))
         }
        
        if(input$select == "fifth stage susceptible nymphs"){
            plot(x~z)
            fit2 <- nls(z ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, z))
        }
        
        
        times2 <- seq(0, 20000, by = 1)
        signal <- data.frame(times = times2, import = rep(0, length(times2)))
        
        vector_trt = c()
        if (input$number_treatments>0) {
            vector_trt = (1:input$number_treatments-1)*input$days_between_treatments}
        days_vector =  c()
        if (input$number_treatments>0) { 
            days_vector = input$start_day+vector_trt}
        
        
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
        
        day <-trt2.function(fit2, days_vector, signal)
        input <- approxfun(day, rule = 2)
        
        
        #=============================================
        #   Set initial conditions for ODE
        #=============================================
        y0 <- c(X = 0.01, Y = 0, m=m)
        parameters <- list(a=a, b=b, c=c, m=m, n=n, g=g, r= 1/(r*365), R=R, K=m, p=p, k=k, MM=MM)
        times <- seq(from=0, 365*nyrs+1, by=1 )
        
        SetODEs <- function(times, y, parameters){
            X <- y[1]
            Y <- y[2]
            m <- y[3]
            with(
                as.list(c(parameters)), 
                {
                    z <- input(times)
                    dX <- (((m+MM)*a*b*Y)+(p*k*(m*a*z*Y)))*(1-X)-r*X
                    dY <- a*c*X*(exp(-g*n)-Y)-((g*Y)+(m*a*z*Y))+(MM*a*c*X*(exp(-g*n)-Y))
                    #dY <- a*c*X*(exp(-g*n)-Y)-((g*(1-(m+MM)/K)*Y)+(m*a*z*Y))+(MM*a*c*X*(exp(-g*n)-Y)) 
                    dm <- ((R*(1-(m+MM)/K)*(m+MM))+(-m*a*z))
                    return(list(c(dX, dY, dm)))
                }
            )
        }
        
        ###### Run the ode solver function
        out <- as.data.frame(ode(y = y0, times = times, func = SetODEs, parms = parameters))
        RESULTS<-data.frame(out$X,out$Y)
        
        
        ###### Plot the output
        matplot(times/365, RESULTS, type = "l", xlab = "Time", ylab = "Proportion (X,Y)", main = "Model with density dependent death for bugs", lwd = 1, lty = 1, bty = "l", col = c("black","red"))
        legend("bottomright", c("Dogs", "Triatomines"), pch = 1, col = c("black","red"))
        
    })
    
    
    
    
    
}