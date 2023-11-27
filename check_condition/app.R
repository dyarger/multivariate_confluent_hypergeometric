Rcpp::sourceCpp('../mch.cpp')
library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Visualization of validity conditions"),

    # Sidebar
    sidebarLayout(
      sidebarPanel(
        column(width = 6,
               sliderInput("nu1",
                           shiny::withMathJax("$$\\nu_{1}:$$"),
                           min = .1,
                           max = 8.1,
                           value = 2.5, step = .1),
               sliderInput("a1",
                           shiny::withMathJax("$$\\alpha_{1}:$$"),
                           min = .5,
                           max = 15.5,
                           value = 2.5, step = .1),
               sliderInput("beta1",
                           shiny::withMathJax("$$\\beta_{1}:$$"),
                           min = .1,
                           max = 8.1,
                           value = 2.5, step = .1)),
        column(width = 6,
               sliderInput("nu2",
                           shiny::withMathJax("$$\\nu_{2}:$$"),
                           min = .1,
                           max = 8.1,
                           value = 2.5, step = .1),
               sliderInput("a2",
                           shiny::withMathJax("$$\\alpha_{2}:$$"),
                           min = .5,
                           max = 15.5,
                           value = 2.5, step = .1),
               sliderInput("beta2",
                           shiny::withMathJax("$$\\beta_{2}:$$"),
                           min = .1,
                           max = 8.1,
                           value = 2.5, step = .1)),
        fluidRow(column(width = 6,
               radioButtons('take_nu', 
                            shiny::withMathJax("Force $$\\nu_{12} = (\\nu_1  + \\nu_2)/2:$$"),
                            #'Force nu12 = (nu1 + nu2)/2:', 
                            choices = c('Yes', 'No'),
                            selected = 'Yes')), 
        column(width = 6,
               conditionalPanel("input.take_nu == 'No'", sliderInput("nu12",
                                                                     shiny::withMathJax("$$\\nu_{12}:$$"),
                                                                        min = .1,
                                                                        max = 8.1,
                                                                        value = 2.5, step = .1)))),
        fluidRow(column(width = 6,               radioButtons('take_alpha', 
                                                              shiny::withMathJax("Force $$\\alpha_{12} = (\\alpha_1  + \\alpha_2)/2:$$"),
                                                     choices = c('Yes', 'No' = 'No'),
                                                     selected = 'Yes')),
        column(width = 6,
               conditionalPanel("input.take_alpha == 'No'", sliderInput("alpha12",
                                                                        shiny::withMathJax("$$\\alpha_{12}:$$"),
                           min = .5,
                           max = 15.5,
                           value = 2.5, step = .1)))),
        fluidRow(column(width = 6,               radioButtons('take_beta', 
                                                     'Force $$\\beta_{12}^2 = (\\beta_1^2  + \\beta_2^2)/2:$$', 
                                                     choices = c('Yes', 'No'),
                                                     selected = 'Yes')),
                 column(width = 6,conditionalPanel("input.take_beta == 'No'", sliderInput("beta12",
                                                                                          shiny::withMathJax("$$\\beta_{12}:$$"),
                                                                          min = .5,
                                                                          max = 15.5,
                                                                          value = 2.5, step = .1)))),
        sliderInput("d",
                    "d:",
                    min = 1,
                    max = 10,
                    value = 1, step = 1),
      ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot", height = '600px', width = '1000px')
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
      a1 <- input$a1 # alphas
      a2 <- input$a2
      if (input$take_alpha == 'Yes') {
        a12 <- (a1 + a2)/2  
      } else {
        a12 <- input$alpha12
      }
      beta1 <- input$beta1
      beta2 <- input$beta2
      if (input$take_beta == 'Yes') {
        beta12 <- sqrt((beta1^2 + beta2^2)/2) 
      } else {
        beta12 <- input$beta12
      }
      nu1 = input$nu1
      nu2 = input$nu2
      if (input$take_nu == 'Yes') {
        nu12 <- (nu1 + nu2)/2  
      } else {
        nu12 <- input$nu12
      }
      d <- input$d
      b_function <- base::beta
      
      x <- 10^(seq(-8, 8, length.out = 1500))
      spec <- mix <- rep(0, length(x))
      for (i in 1:length(x)) {
        spec[i] <- (beta12/sqrt(beta1*beta2))^d * 
          sqrt(b_function(nu1, a1) * b_function(nu2, a2))/(b_function(nu12, a12)) *# exact spectral density value
          gamma(nu12 + d/2)/(sqrt(gamma(nu1 + d/2) * gamma(nu2 + d/2))) *
          Ufun(nu12 + d/2, 1 - a12 + d/2, beta12^2*x[i]^2/2) /
          sqrt(Ufun(nu1 + d/2, 1 - a1 + d/2, beta1^2 * x[i]^2/2) *
                 Ufun(nu2 + d/2, 1 - a2 + d/2, beta2^2 * x[i]^2/2))  
        mix[i] <- sqrt(gamma(nu1)*gamma(nu2))/gamma(nu12) /  # exact form when using mixture representation
          sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * gamma(nu12 + d/2) * # depends on x/phi
          beta12^(2*a12)/(beta1^a1 * beta2^a2) * sqrt(gamma(a1) * gamma(a2))/gamma(a12) * 
          x[i]^(-2*a12 + a1 + a2) * 
          exp(-1/2/x[i]^2*(beta12^2 - beta1^2/2 - beta2^2/2))
      }
      prop3.3 <- 1/(b_function(nu12, a12)) * # bound presented in Prop 3.3
        sqrt(b_function(nu1, a1) * b_function(nu2, a2)) 
      thm3.1 <- gamma(nu12 + d/2) / sqrt(gamma(nu1 + d/2) * gamma(nu2 + d/2)) / # bound presented in Theorem 3.1
        gamma(nu12) * sqrt(gamma(nu1) * gamma(nu2)) /
        gamma(a12) * sqrt(gamma(a1) * gamma(a2)) *
        beta12^(2*a12)/beta1^(a1)/beta2^(a2)
      thm3.2 <- sqrt(gamma(nu1)*gamma(nu2))/gamma(nu12) * # bound presented in Theorem 3.2
        sqrt(gamma(a1)*gamma(a2))/gamma(a12) *
        (nu12)^(nu12 + d/2) / ((nu1)^(nu1/2 + d/4)*(nu2)^(nu2/2 + d/4)) * 
        exp(-nu12 + nu1/2 + nu2/2) * 
        beta12^(2*a12)/beta1^(a1)/beta2^(a2)
      # if (max(1/mix) > 10 & max(1/spec) < 10) {
      #   range_use <- range(c(0, (1/spec)[1/spec < 2], 1/prop3.3, 1/thm3.1, 1/thm3.2))
      # } else if (max(1/spec) > 10) {
      #   range_use <- range(c(0, 1/prop3.3, 1/thm3.1, 1/thm3.2))
      # } else {
      #   range_use <- range(c(0, 1/spec, 1/mix, 1/prop3.3, 1/thm3.1, 1/thm3.2))
      # }
      range_use <- range(c(0,  (1/mix)[1/mix < 2],(1/spec)[1/spec < 2], 1/prop3.3, 1/thm3.1, 1/thm3.2))
      
      if (beta1 == beta2 & beta12 == beta1 & a1 > d/2 & a2 > d/2 & 
          a12 >= (a1 + a2)/2 & nu12 == (nu1 + nu2)/2) {
        prop3.3_cond <- T
      } else {
        prop3.3_cond <- F
      }
      if (beta12^2 - beta1^2/2 - beta2^2/2 > -10^-5 &
          a12 == (a1 + a2)/2 & nu12 == (nu1 + nu2)/2) {
        thm3.1_cond <- T
      } else {
        thm3.1_cond <- F
      }
      if (beta12 >= (beta1 + beta2)/2 &
          a12 == (a1 + a2)/2 & nu12 >= (nu1 + nu2)/2) {
        thm3.2_cond <- T
      } else {
        thm3.2_cond <- F
      }
      if (a1 > d/2 & a2 > d/2) {
        spec_cond <- T
      } else {
        spec_cond <- F
      }
      if (nu12 >= (nu1 + nu2)/2) {
        mix_cond <- T
      } else {
        mix_cond <- F
      }
      levels_use <- c('Spectral density', 'Mixing density',
                      'Theorem 3.1', 'Theorem 3.2', 'Proposition 3.3')
      df_varying <- data.frame(x = log10(x), 
                               values = c(1/spec, 1/mix, rep(c(1/prop3.3, 1/thm3.1, 1/thm3.2), each = length(x))), 
                               type = factor(rep(c('Spectral density', 'Mixing density', 'Proposition 3.3', 'Theorem 3.1', 'Theorem 3.2'), each = length(x)),
                                             levels = levels_use),
                               'Conditions met' = rep(paste('Condition met:', c(spec_cond, mix_cond, prop3.3_cond, thm3.1_cond, thm3.2_cond)), each = length(x)))
      ggplot(mapping = aes(color = factor(type, levels = levels_use),
                linetype = factor(type, levels = levels_use))) + 
        geom_line(data = df_varying, aes(y = values, x = x), size = 1.5) +
        facet_wrap(~`Conditions.met`) +
        scale_y_continuous(limits = range_use, breaks = seq(0, 2.5, by = .25)) + 
        labs(y = 'Maximal correlation in p=2 case',
             x = expression(log(phi)), color = 'Type', linetype = 'Type') +
        theme_bw() +
        theme(legend.position = 'bottom', text = element_text(size = 20))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
