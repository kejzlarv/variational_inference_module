#setwd("C:/Users/Vojtech/Dropbox/Bayes VI ed paper/Stan")
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)

# Checking integrity of installation of cmdstanr
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
cmdstan_path()
cmdstan_version()

# Auxiliary packages
library(tm)
library(tidyverse)
library(tidytext)
library(topicmodels)

## Get data
data("AssociatedPress", package = "topicmodels")

## Removing rare words from the vocabulary
dtm = removeSparseTerms(AssociatedPress, 0.95)
dim(dtm)

## Input for stan model
N_TOPICS = 2

data = list(K = N_TOPICS,
            V = dim(dtm)[2],
            M = dim(dtm)[1],
            N = sum(dtm$v),
            w = rep(dtm$j,dtm$v),
            doc = rep(dtm$i,dtm$v),
            alpha = rep(50/N_TOPICS,N_TOPICS), #according to Griffiths and Steyvers(2004) 
            beta = rep(1,dim(dtm)[2])  #according to Griffiths and Steyvers(2004) 
)

### VB fit
LDA_model_cmd <- cmdstan_model(stan_file = "LDA.stan")
LDA_model_cmd$print()  

vb_fit <- LDA_model_cmd$variational(data = data,
                                    seed = 1,
                                    output_samples = 1000,
                                    eval_elbo = 1,
                                    grad_samples = 10,
                                    elbo_samples = 10,
                                    algorithm = "meanfield",
                                    output_dir = NULL,
                                    iter = 1000,
                                    adapt_iter = 20,
                                    save_latent_dynamics=TRUE,
                                    tol_rel_obj = 10^-4) 

# Plotting ELBO
vb_diag <- utils::read.csv(vb_fit$latent_dynamics_files()[1], comment.char = "#")
ELBO = data.frame(Iteration = vb_diag[,1], ELBO = vb_diag[,3])

ggplot(data = ELBO, aes(x = Iteration, y = ELBO)) + geom_line(lwd=1.5) + 
  theme(text = element_text(size = 20), panel.background = element_rect(fill = "transparent", color = "lightgrey"), panel.grid.major = element_line(colour = "lightgrey")) +
  xlim(0,110)

## Accessing parameters

vb_fit$summary("theta") # dim: M-by-K
vb_fit$summary("phi") # dim: K-by-V

## Word distribution per topic

V <- dim(dtm)[2]
odd_rows <- rep(c(1,0), times = V)
Topic1 <- vb_fit$summary("phi")[odd_rows == 1,]
Topic2 <- vb_fit$summary("phi")[odd_rows == 0,]

word_probs <-data.frame(Topic = c(rep("Topic 1", V), rep("Topic 2", V)),
                        Word = rep(dtm$dimnames$Terms,N_TOPICS),
                        Probability = c(Topic1$mean, Topic2$mean))

# Selecting top 10 words per topic
top_words <- word_probs %>% group_by(Topic) %>% top_n(10) %>% ungroup() %>% arrange(Topic, -Probability)

top_words %>%
  mutate(Word = reorder_within(Word, Probability, Topic)) %>%
  ggplot(aes(Probability, Word, fill = factor(Topic))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ Topic, scales = "free") +
  scale_y_reordered() + theme(text = element_text(size = 15)) + xlim(0,0.025) +
  xlab("Word distributions ( \u03d5 )")

# Word Cloud display
#install.packages("wordcloud")
library(wordcloud)


top_words <- word_probs %>% group_by(Topic) %>% top_n(20) %>% ungroup() %>% arrange(Topic, -Probability)

mycolors <- brewer.pal(8, "Dark2")
wordcloud(top_words %>% filter(Topic == "Topic 1") %>% .$Word ,
          top_words %>% filter(Topic == "Topic 1") %>% .$Probability,
          random.order = FALSE,
          color = mycolors)

mycolors <- brewer.pal(8, "Dark2")
wordcloud(top_words %>% filter(Topic == "Topic 2") %>% .$Word ,
          top_words %>% filter(Topic == "Topic 2") %>% .$Probability,
          random.order = FALSE,
          color = mycolors)

