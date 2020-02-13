#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # App title ----
    titlePanel("Chagas Disease Bugs and Dogs Transmission Model, Prior to treatment"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Slider for the number of bins ----
            sliderInput(inputId = "years",
                        label = "Number of years to run simulation",
                        min = 1, max = 100, value = 20),
            sliderInput(inputId = "m",
                        label = "Ratio of bugs:dogs at beginning of the simulation",
                        min = 5, max = 100, value = 40),
            sliderInput(inputId = "b",
                        label = "Transmission efficiency from bug to dog",
                        min = 0.0005, max = 0.001, value = 0.00068),
            sliderInput(inputId = "c",
                        label = "Transmission efficiency from dog to bug",
                        min = 0.1, max = 0.5, value = 0.28),
            sliderInput(inputId = "g",
                        label = "Daily probabilty of vector mortality",
                        min = 0.001, max = 0.01, value = 0.005),
            sliderInput(inputId = "a",
                        label = "Bug bite rate",
                        min = 1/21, max = 1/7, value = 1/14),
            sliderInput(inputId = "r",
                        label = "Lifespan of dog in years",
                        min = 1, max = 10, value = 3),
            sliderInput(inputId = "R",
                        label = "Bug birthrate",
                        min = 0.05, max = 0.12, value = 0.09),
            sliderInput(inputId = "K",
                        label = "Bug Carrying Capacity",
                        min = 5, max = 80, value = 40),
            
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Proportions infected (Dogs and Bugs) ----
            plotOutput(outputId = "distPlot")
            
        )
    )
)
