# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)

#setwd("...")

# Checking integrity of installation of cmdstanr
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
cmdstan_path()
cmdstan_version()

# Auxiliary packages
library(tidyverse)

## Get data
PSID <- read.csv("https://bit.ly/LaborVar")

## Input for stan model

data = list(N = dim(PSID)[1],
            x = PSID$FamilyIncome,
            y = PSID$Participation
)

## VI fit
Logistic_model_cmd <- cmdstan_model(stan_file = "Logistic.stan")
Logistic_model_cmd$print()  

vi_fit <- Logistic_model_cmd$variational(data = data,
                                    seed = 1,
                                    output_samples = 5000,
                                    eval_elbo = 1,
                                    grad_samples = 15,
                                    elbo_samples = 15,
                                    algorithm = "meanfield",
                                    output_dir = NULL,
                                    iter = 1000,
                                    adapt_iter = 20,
                                    save_latent_dynamics=TRUE,
                                    tol_rel_obj = 10^-4) 

## Plotting ELBO
vb_diag <- utils::read.csv(vi_fit$latent_dynamics_files()[1], comment.char = "#")
ELBO = data.frame(Iteration = vb_diag[,1], ELBO = vb_diag[,3])

ggplot(data = ELBO, aes(x = Iteration, y = ELBO)) + geom_line(lwd=1.5) + 
  theme(text = element_text(size = 20), panel.background = element_rect(fill = "transparent", color = "lightgrey"), panel.grid.major = element_line(colour = "lightgrey")) +
  xlim(0,110)

## Accessing parameters
vi_fit$summary("alpha") # 
vi_fit$summary("beta") #

## Posterior interval estimates for the probability of labor participation 

prob_interval <- function(x, post){
  lp <- post[, 1] + x * post[, 2]
  quantile(exp(lp) / (1 + exp(lp)),
           c(.05, .50, .95))
}

out <- sapply(seq(10, 70, by = 10),
              prob_interval, vi_fit$draws()[,3:4])

df_out <- data.frame(Income = seq(10, 70, by = 10),
                     Low = out[1, ],
                     M = out[2, ],
                     Hi = out[3,  ])


ggplot(df_out) +
  geom_line(aes(x = Income, y = M), lwd=1.5) +
  geom_segment(aes(x = Income, y = Low,
                   xend = Income, yend = Hi), size = 2) +
  ylab("Prob(Participate)") +
  ylim(0, 1) + 
  theme(text = element_text(size = 20), panel.background = element_rect(fill = "transparent", color = "lightgrey"),
        panel.grid.major = element_line(colour = "lightgrey"))

## MCMC fit

mcmc_fit <- Logistic_model_cmd$sample(
  data = data, 
  seed = 1, 
  chains = 1, 
  iter_warmup = 5000,
  iter_sampling = 5000
)

mcmc_fit$summary()

## VI/MCMC comparison

df <- data.frame(Method = c(rep('VI', 5000), rep('MCMC', 5000) ),
                 value = c(vi_fit$draws("alpha"),
                           mcmc_fit$draws("alpha")))

ggplot(df, aes(x=value, fill=Method)) +
  geom_histogram( color='#e9ecef', alpha=0.6, position='identity') + xlab("alpha")

df <- data.frame(Method = c(rep('VI', 5000), rep('MCMC', 5000) ),
                 value = c(vi_fit$draws("beta"),
                           mcmc_fit$draws("beta")))

ggplot(df, aes(x=value, fill=Method)) +
  geom_histogram( color='#e9ecef', alpha=0.6, position='identity') + xlab("beta")

## Speed comparison
# VI
PSID_rep <- PSID[rep(seq_len(nrow(PSID)), 50), ]

data = list(N = dim(PSID_rep)[1],
            x = PSID_rep$FamilyIncome,
            y = PSID_rep$Participation
)

vi_fit <- Logistic_model_cmd$variational(data = data,
                                         seed = 1,
                                         output_samples = 5000,
                                         eval_elbo = 1,
                                         grad_samples = 15,
                                         elbo_samples = 15,
                                         algorithm = "meanfield",
                                         output_dir = NULL,
                                         iter = 1000,
                                         adapt_iter = 20,
                                         save_latent_dynamics=TRUE,
                                         tol_rel_obj = 10^-4) 

# MCMC
mcmc_fit <- Logistic_model_cmd$sample(
  data = data, 
  seed = 1, 
  chains = 1, 
  iter_warmup = 5000,
  iter_sampling = 5000
)

