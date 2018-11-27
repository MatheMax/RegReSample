#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
# install_github('MatheMax/RegReSample')
library(RegReSample)

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Optimal Two-Stage Designs"),

   # Sidebar with a slider input for number of bins
   #sidebarLayout(
      fluidRow(
        column(3,
         sliderInput("alpha",
                     "Maximal Type One Error Rate",
                     min = .01,
                     max = .1,
                     step = .005,
                     value = .025),
         sliderInput("beta",
                     "Maximal Type Two Error Rate",
                     min = .05,
                     max = .4,
                     step = .01,
                     value = .2),
         sliderInput("beta.2",
                     "Minimal Conditional Power",
                     min = .05,
                     max = .95,
                     value = .5,
                     step = .01),
         sliderInput("delta.alt",
                     "Expected effect size",
                     min = .1,
                     max = 1,
                     value = .3,
                     step = .01)
        ),
        column(3,
         sliderInput("lambda1",
                     "Penalty for Expected Conditional Power",
                     min = 0.0,
                     max = 10,
                     value = 0.0,
                     step = 1),
         sliderInput("lambda2",
                     "Penalty for norm of derivative of n_2",
                     min = 0.0,
                     max = 5,
                     value = 0.0,
                     step = .1),
         sliderInput("lambda3",
                     "Penalty for derivative of Conditional Power",
                     min = 0.0,
                     max = 10,
                     value = 0.0,
                     step = 1)
         ),
        column(3,
         selectInput("prior",
                     h4("What is the prior distribution?"),
                     choices = list("normal", "point"),
                     selected = "point"),
         checkboxInput("pred.yn",
                       "Should the prior distribution be updated at interim?",
                       value = FALSE),
         sliderInput("delta.mcr",
                     "Minimal clinically relevant effect size",
                     min = .05,
                     max = .5,
                     value = .1,
                     step = .01),
         sliderInput("tau",
                     "Standard deviation of prior distribution",
                     min = .05,
                     max = .5,
                     value = .01,
                     step = .01)
        ),
        column(3,
         checkboxInput("n.max.yn",
                       "Should a maximal sample size be used?",
                       value = FALSE),
         sliderInput("n.max.v",
                     "Maximal sample size per stage and group",
                     min = 10,
                     max = 500,
                     value = 200,
                     step = 1),
         checkboxInput("n1.yn",
                       "Should the first-stage sample size be fixed?",
                       value = FALSE),
         sliderInput("n1.v",
                     "First-stage sample size per group",
                     min = 10,
                     max = 200,
                     value = 100,
                     step = 1),
         selectInput("distribution",
                     h4("What is the test statistic's distribution?"),
                      choices = list("normal", "t"),
                     selected = "normal")
      )
      ),

      hr(),

      submitButton("Compute!"),

     hr(),

     # mainPanel(
         plotOutput("design", width = "100%", height = "400px")
       #  )
   #)
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #reactive({
  output$design <- renderPlot({
    withProgress(message = "Computing design",{
    pow.pred <- ifelse(input$pred.yn == T, "predictive", "equal")
    n.max <- ifelse(input$n.max.yn == T, input$n.max.v, Inf)
    n1 <- ifelse(input$n1.yn == T, input$n1.v, NA)

    j <- RegReSample::opt_design(input$alpha,
                                 input$beta,
                                 1 - input$beta.2,
                                 c(input$lambda1, input$lambda2, input$lambda3),
                                 input$prior,
                                 pow.pred,
                                 input$delta.mcr,
                                 input$delta.alt,
                                 input$tau,
                                 n.max,
                                 input$distribution,
                                 n1
                                 )

    return(RegReSample::plot_design(j, input$distribution, input$delta.alt,
                                    input$prior, pow.pred, input$delta.mcr,
                                    input$beta, 1 - input$beta.2,
                                    input$tau))}
  )})
}

# Run the application
shinyApp(ui = ui, server = server)

