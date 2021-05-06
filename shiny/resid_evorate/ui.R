# UI
library(shiny)
#library(dashboard)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput("properties", "Property:",choices = unique(prop$property),selected = c('essential_core','essential_dispensable'),selectize = T, multiple = T),
            selectInput("features", "Features:",choices = unique(feat$feature),selectize = T, multiple = T),
        ),

        # Show a plot of the generated distribution
        mainPanel(
            column(4,plotOutput("propfit")),
            column(4,plotOutput("propres")),
            column(4,plotOutput("sumres")),

            column(4,plotOutput("featfit")),
            column(4,plotOutput("featres")),
            column(4,plotOutput("corres"))
        )
    )
))



## Residual evolutionary rate
```{r shiny, echo = FALSE}
library(shiny)
library(plotly)


shinyApp(

    ),

    server = function(input, output) {


    },
)


```
