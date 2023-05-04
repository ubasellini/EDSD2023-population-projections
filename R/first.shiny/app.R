#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)

my.years <- unique(my.proj.long$year)
my.cols <- cividis(n+1)
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("First Shiny App of Population Projections"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("year",
                        "Year:",
                        min = min(my.years),
                        max = max(my.years),
                        value = min(my.years),
                        step=5)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotlyOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlotly({
        ggplot(data=my.proj.long,
               aes(x=AgeGroup,y=population,fill=yearF)) +
            geom_bar(data = subset(my.proj.long, year == input$year),
                     stat="identity",position = "dodge",color="black") +
            coord_flip() +
            ggtitle("Swedish female population") +
            theme_bw() + ylim(0,max(my.proj.long$population)) +
            scale_fill_manual(name="Year",values=my.cols[which(input$year==my.years)])
    })
}

# Run the application 
shinyApp(ui = ui, server = server)






