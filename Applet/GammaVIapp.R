library(shiny)
library(tidyverse)
library(reshape2)
library(cmdstanr)
library(r2symbols)
library(shinyjs)


## Gamma Poisson example for ADVI in cmdstanr
# EXAMPLE: Estimate of the average number of active users of a popular massively multiplier online role-playing game (mmorpg) 
# playing between the peak evening hours 7 pm and 10 pm.
y <- c(50, 47, 46, 52, 49, 55, 53, 48, 45, 51, 50, 53, 46, 47)
N <- length(y)
## Prior
a_pri <- 100
b_pri <- 2
# Closed form posterior is also gamma
a_pos <- a_pri + sum(y)
b_pos <- b_pri + N
# PLot prior VS posterior
x <- seq(from = 35, to = 65, by = 0.1)
distr <- data.frame(x = x,
                    prior = dgamma(x, shape = a_pri, rate = b_pri),
                    posterior = dgamma(x, shape = a_pos, rate = b_pos))


ELBO <- function(y, mu_q, sigma_q, a_pri, b_pri, MC_n = 2000){
  q_sample <- rlnorm(MC_n, meanlog = mu_q, sdlog = sigma_q)
  ELBO <- 0
  for (i in 1:MC_n){
    expectation <- sum(log(dpois(y, q_sample[i]))) + log(dgamma(q_sample[i], shape = a_pri, rate = b_pri)) - log(dlnorm(q_sample[i],  meanlog = mu_q, sdlog = sigma_q))
    ELBO <- ELBO + expectation/MC_n
  }

  ELBO
}

Relable <- function(dataset){
  replace(dataset$variable,dataset$variable=="posterior","True Posterior")
}

#### EVIDENCE Computaiton
factorial_y <- 0
for (i in 1:N) {
  factorial_y <- factorial_y + sum(log(1:y[i]))
}

log_evidence <- a_pri * log(b_pri) + lgamma(sum(y) + a_pri) - lgamma(a_pri) - factorial_y - (sum(y) + a_pri) * log(N+ b_pri)
#####


ui <- fluidPage(theme = bslib::bs_theme(bootswatch = "sandstone"),
  titlePanel(
    "Variational Inference with Gamma-Poisson Model for count data"
  ),
  fluidRow(
    column(4,
           "Variational approximation using log-normal variational family:",
           shiny::uiOutput("qm_text"),
           sliderInput("qm", "", min = 3.5, max = 4.2, value = 3.7, step = 0.01),
           shiny::uiOutput("qs_text"),
           sliderInput("qs", "" , min = 0.01, max = 0.1, value = 0.05),
           textOutput("ELBO_text"),
           verbatimTextOutput("ELBO"),
           textOutput("KL_text"),
           verbatimTextOutput("KL")),


  column(8,
         plotOutput("distPlot")
    )
  ),
  fluidRow(useShinyjs(),
    column(4,
           checkboxInput("elbo_check", "Fit a variational approximation"),
           textOutput("KL_text_fit"),
           shiny::uiOutput("vb_mu_text"),
           verbatimTextOutput("VB_fit_mu"),
           shiny::uiOutput("vb_sd_text"),
           verbatimTextOutput("VB_fit_sd"),
           textOutput("VB_fit_elbo"),
           verbatimTextOutput("VB_fit_elbo_text")
             ),


    column(8,
           plotOutput("elboPlot")
    )
  )
)

server <- function(input, output, session) {

  data <- list(y = y, N = N)
  poisson_model_cmd <- cmdstan_model(stan_file = "GammaPoissonLog.stan")

  distr_with_vb <- reactive(distr %>% mutate(vb = dlnorm(x, meanlog = input$qm, sdlog = input$qs)))
  ELBO_computation <- reactive(ELBO(y, input$qm, input$qs, a_pri, b_pri, MC_n = 2000))
  
  output$qm_text <- renderUI({
    s<-"\\( \\mu \\)"
    withMathJax(s)        
    
  })

  output$qs_text <- renderUI({
    s<-"\\( \\sigma \\)"
    withMathJax(s)        
    
  })
  output$distPlot <- renderPlot({
    
    if (input$elbo_check) {
      req(vbdraws())
      distr_with_vb_new <- distr_with_vb() %>% mutate(vb_fit = dlnorm(x, meanlog = vbdraws()$summary("lambda_real")$mean, sdlog = vbdraws()$summary("lambda_real")$sd))
      distr_melt <- melt(distr_with_vb_new, id = "x")
      levels(distr_melt$variable) <- c("Prior", "True posterior", "VI - Manual", "VI - ELBO maximization")
      ggplot(data = distr_melt, aes(x = x, y = value)) + 
        xlim(35, 65) + 
        ylim(0, 0.3) +
        geom_line( aes(linetype=variable, color = variable, size = variable)) +
        theme(text = element_text(size = 25), legend.position = "top", panel.background = element_rect(fill = "transparent", color = "lightgrey"), panel.grid.major = element_line(colour = "lightgrey"), legend.title=element_blank()) +
        ylab("PDF") +
        xlab(expression(theta)) +
        scale_linetype_manual(values=c("dotted","solid", "twodash", "longdash"))+ 
        scale_size_manual(values=c(1.5, 3, 1.5, 1.5))  +
        scale_color_manual(values = c("brown1", "chartreuse3", "blue3", "black"))
      
    } else {
      distr_melt <- melt(distr_with_vb(), id = "x")
      levels(distr_melt$variable) <- c("Prior", "True posterior", "VI - Manual")
    ggplot(data = distr_melt, aes(x = x, y = value, color = variable)) + 
      xlim(35, 65) +
      ylim(0, 0.3) +
      geom_line(aes(linetype=variable, size = variable)) + theme(text = element_text(size = 25), legend.position = "top", panel.background = element_rect(fill = "transparent", color = "lightgrey"), panel.grid.major = element_line(colour = "lightgrey"), legend.title=element_blank()) +
      ylab("PDF") +
      xlab(expression(theta)) +
      scale_linetype_manual(values=c("dotted", "solid", "twodash")) +
      scale_size_manual(values=c(1.5, 3, 1.5)) +
      scale_color_manual(values = c("brown1", "chartreuse3", "blue3"))
    }
  })
  output$ELBO_text <- renderText("ELBO value:")

  output$ELBO <- renderPrint({
    as.html(round(ELBO_computation(), digits = 3))
  })

  output$KL_text <- renderText("KL Divergence value:")

  output$KL <- renderPrint({
     as.html(round(log_evidence - ELBO_computation(), digits = 3))
  })


  observeEvent(input$elbo_check, {
    shinyjs::toggle("elboPlot")
    shinyjs::toggle("KL_text_fit")
    shinyjs::toggle("VB_fit_mu")
    shinyjs::toggle("VB_fit_sd")
    shinyjs::toggle("VB_fit_elbo")
    shinyjs::toggle("VB_fit_elbo_text")
    shinyjs::toggle("vb_sd_text")
    shinyjs::toggle("vb_mu_text")
  })

  vbdraws <- reactive({
    if (input$elbo_check) {
      poisson_model_cmd$variational(data = data, seed = 1, output_samples = 4000, eval_elbo = 1, grad_samples = 10, elbo_samples = 10, algorithm = "meanfield",
                                    output_dir = NULL, iter = 1000, save_latent_dynamics=TRUE) 
    }
  })


output$elboPlot <- renderPlot({
    req(vbdraws())
    vb_diag <- utils::read.csv(vbdraws()$latent_dynamics_files()[1], comment.char = "#")
    ELBO = data.frame(Iteration = vb_diag[,1], ELBO = vb_diag[,3])
    #file.remove(vbdraws()$latent_dynamics_files()[1])
    ggplot(data = ELBO, aes(x = Iteration, y = ELBO)) +
      geom_line(lwd=1.5) + theme(text = element_text(size = 25), panel.background = element_rect(fill = "transparent", color = "lightgrey"), panel.grid.major = element_line(colour = "lightgrey")) + xlim(0,200) + ylim(-200, 0)
  })

output$KL_text_fit <- renderText("Results of ELBO maximization via gradient ascent:")

output$VB_fit_mu <- renderPrint({
  req(vbdraws())
  as.html(paste0( as.html(round(vbdraws()$summary("lambda_real")[2], digits = 3))))
})

output$VB_fit_sd <- renderPrint({
  req(vbdraws())
  as.html(paste0( as.html(round(vbdraws()$summary("lambda_real")[4], digits = 3))))
})

output$VB_fit_elbo <- renderText({
  paste0("ELBO value:")
})

output$VB_fit_elbo_text <- renderPrint({
  req(vbdraws())
  vb_diag <- utils::read.csv(vbdraws()$latent_dynamics_files()[1], comment.char = "#")
  ELBO = vb_diag[,3]
  as.html(round(ELBO[length(ELBO)], digits = 3))
})

output$vb_mu_text <- renderUI({
  s<-"\\( \\mu \\)"
  withMathJax(s)        
  
})


output$vb_sd_text <- renderUI({
  s<-"\\( \\sigma \\)"
  withMathJax(s)        
  
})


}

shinyApp(ui, server)
