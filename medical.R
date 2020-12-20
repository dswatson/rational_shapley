# Set working directory
setwd('~/Documents/rational_shapley')

# Load libraries
library(data.table)
library(tidyverse)
library(elasticnet)
library(glmnet)
library(caret)
library(shapr)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import data
data(diabetes)
x <- as.matrix(cbind(diabetes$x[, c(2, 1, 3:10)]))
y <- diabetes$y
df <- as.data.frame(cbind(x, y))

# Fit model
f <- train(y ~ ., data = df, method = 'glmnet', 
           trControl = trainControl(method = 'cv'),
           tuneLength = 100)
beta <- as.numeric(coef(f$finalModel, f$bestTune$lambda))
f <- lm(y ~ ., data = df)
f$coefficients <- beta

# Focus on poorest responders
x_i <- df %>%
  mutate(y_hat = predict(f, df)) %>%
  filter(y_hat >= quantile(y_hat, 0.9)) %>%
  select(-y, -y_hat) %>%
  as.matrix(.)

# Define reference distribution
explainer <- shapr(x, f)
phi_0 <- mean(y)

# Marginal Shapley values
msv <- explain(x_i, approach = 'empirical', type = 'independence',
               explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Conditional Shapley values
csv <- explain(x_i, approach = c('empirical', rep('gaussian', 9)),
               explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Interventional Shapley values
isv <- explain(x_i, approach = c('empirical', rep('causal', 9)), 
               explainer = explainer, prediction_zero = phi_0,
               ordering = list(c(1, 2), c(3, 4), c(5:10)))$dt %>%
  select(-none)

# Plot results
vals <- c('Marginal', 'Conditional', 'Interventional')
rbind(msv, csv, isv) %>%
  mutate(value = factor(rep(vals, each = 45), levels = vals)) %>%
  pivot_longer(cols = -value, names_to = 'feature', values_to = 'phi') %>%
  ggplot(aes(phi, feature, fill = phi)) + 
  geom_jitter(size = 2, width = 0, height = 0.1, color = 'black', pch = 21) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Shapley\nValue', option = 'B') +
  xlab('Shapley Value') + ylab('Feature') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ value)


### RATIONAL SHAPLEY ### 

# Define the subspace: same sex and age
subspace <- df %>% filter(sex < 0, age < quantile(age, 0.25), y < median(y))
phi_0 <- mean(subspace$y)
x_ref <- subspace %>% 
  select(-y) %>%
  as.matrix(.)
explainer <- shapr(x_ref, f)
phi_0 <- mean(subspace$y)

# Marginal Shapley values
msv <- explain(x_i, approach = 'empirical', type = 'independence',
               explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Conditional Shapley values
csv <- explain(x_i, approach = c('empirical', rep('gaussian', 9)),
               explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Interventional Shapley values
isv <- explain(x_i, approach = c('empirical', rep('causal', 9)), 
               explainer = explainer, prediction_zero = phi_0,
               ordering = list(c(1, 2), c(3, 4), c(5:10)))$dt %>%
  select(-none)

# Plot results
vals <- c('Marginal', 'Conditional', 'Interventional')
rbind(msv, csv, isv) %>%
  mutate(value = factor(rep(vals, each = 45), levels = vals)) %>%
  pivot_longer(cols = -value, names_to = 'feature', values_to = 'phi') %>%
  ggplot(aes(phi, feature, fill = phi)) + 
  geom_jitter(size = 2, width = 0, height = 0.1, color = 'black', pch = 21) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Shapley\nValue', option = 'B') +
  xlab('Shapley Value') + ylab('Feature') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ value)















dat <- rbind(msv, csv, isv) %>%
  mutate(value = rep(c('Marginal', 'Conditional', 'Interventional'), 
                     each = 45)) %>%
  mutate(value = factor(
    value, levels = c('Marginal', 'Conditional', 'Interventional'))
  ) %>%
  pivot_longer(cols = -value, names_to = 'feature', values_to = 'phi')
ggplot(dat, aes(phi, feature, fill = phi)) + 
  geom_jitter(size = 2, width = 0, height = 0.1, color = 'black', pch = 21) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_distiller('Shapley\nValue', type = 'div', palette = 'RdBu', 
                       limits = c(-1, 1) * max(abs(dat$phi))) +
  xlab('Shapley Value') + ylab('Feature') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ value)










