#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # App title ----
    titlePanel("Chagas Model with Treatment"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Slider for the number of bins ----
            sliderInput(inputId = "years",
                        label = "Number of years to run simulation",
                        min = 1, max = 100, value = 40),
            sliderInput(inputId = "m",
                        label = "Ratio of  treated bugs:dogs",
                        min = 5, max = 80, value = 40),
            sliderInput(inputId = "MM",
                        label = "External Bug Pop- not effected by treatment",
                        min= 0, max= 20, value=5),
            sliderInput(inputId = "b",
                        label = "Transmission efficiency from bug to dog, vectorial",
                        min = 0.0005, max = 0.001, value = 0.00068),
            sliderInput(inputId = "c",
                        label = "Transmission efficiency from dog to bug",
                        min = 0.1, max = 0.5, value = 0.28),
            sliderInput(inputId = "g",
                        label = "Daily probabilty of vector mortality",
                        min = 0.001, max = 0.01, value = 0.005),
            sliderInput(inputId = "a",
                        label = "Bug biting frequency",
                        min = 1/21, max = 1/7, value = 1/14),
            sliderInput(inputId = "r",
                        label = "Lifespan of dog in years",
                        min = 1, max = 10, value = 3),
            sliderInput(inputId = "R",
                        label = "Bug birthrate",
                        min = 0.05, max = 0.12, value = 0.09),
            sliderInput(inputId = "p",
                        label = "Percentage of Dead bugs consumed by dog",
                        min = 0.02, max = 0.99, value = 0.8),
            sliderInput(inputId = "k",
                        label = "Transmission efficiency from bug to dog, oral transmission",
                        min = 0.05, max = 0.12, value = 0.1),
            sliderInput(inputId = "d",
                        label = "Day to implement treatment",
                        min = 8000, max = 19000, value = 10000),
            numericInput(inputId = "start_day",
                         label = "Day to implement treatment",
                         min = 10000, 10000),
            numericInput(inputId = "number_treatments",
                         label = "Number of treatments",
                         value = 1),
            numericInput(inputId = "days_between_treatments",
                         label = "Days between treatments",
                         value=90),

            
            
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            
            # Output: Proportions infected (Dogs and Bugs) ----
            plotOutput(outputId = "distPlot1"),
            #  plotOutput(outputId = "distPlot2"),  #Commented this out bc there's a bug in this graph- population crashes
            
        )
    )
)
