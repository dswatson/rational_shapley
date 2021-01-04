# Set working directory
setwd('~/Documents/rational_shapley')

# Load libraries
library(data.table)
library(tidyverse)
library(e1071)
library(shapr) 
library(ggsci)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import data
df <- fread('german.data')
colnames(df) <- c('chk_acct', 'duration', 'credit_his', 'purpose', 'amount', 
                  'savings2', 'present_emp', 'installment_rate', 'sex', 
                  'other_debtor', 'present_resid', 'property', 'age', 
                  'other_install', 'housing', 'n_credits', 'job2', 'n_people', 
                  'telephone', 'foreign', 'response')

# Recode
df[, gender := ifelse(sex %in% c('A91', 'A93', 'A94'), 1, 0)]
df[, marital := ifelse(sex %in% c('A93', 'A95'), 1, 0)]
df[savings2 == 'A65', savings := 0]
df[savings2 == 'A61', savings := 100]
df[savings2 == 'A62', savings := 500]
df[savings2 == 'A63', savings := 1000]
df[savings2 == 'A64', savings := 2000]
df[job2 == 'A171', job := 0]
df[job2 == 'A172', job := 1]
df[job2 == 'A173', job := 2]
df[job2 == 'A174', job := 3]

# Approximate normality
df[, amount := log(amount)]
df[, age := log(age)]
df[, savings := log(savings + 1) + rnorm(.N)]
df[, job := job + rnorm(.N, sd = 0.5)]
df[, duration := log(duration) + rnorm(.N, sd = 0.1)]
df[, y := ifelse(response == 2, 0, 1)]

# Restrict focus to eight features
df <- df %>% select(gender, marital, age, amount, duration, savings, job, y)

# Fit model
f <- svm(y ~ ., data = df, type = 'C-classification', probability = TRUE)
y_hat <- predict(f, df, probability = TRUE)
y_hat <- attr(y_hat, 'probabilities')[, 1]

# UCI says there's asymmetric costs
df[, y_hat := y_hat]
tau <- seq(0, 1, 0.01)
err_fn <- function(tau) {
  df[, y_cat := ifelse(y_hat >= tau, 1, 0)]
  fp <- df[y == 0, sum(y_cat == 1)]
  fn <- df[y == 1, sum(y_cat == 0)]
  err <- 5 * fp + fn
  return(err)
}
loss <- sapply(tau, err_fn)
(thresh <- tau[which.min(loss)])
tmp <- data.table(tau, loss)
ggplot(tmp, aes(tau, loss)) + 
  geom_line() + 
  geom_hline(yintercept = tmp[tau == thresh, loss], 
             linetype = 'dashed', color = 'red') +
  theme_bw()

# Custom functions for new model class
model_type.svm <- function(x) {
  if (x$type == 0) {
    'classification'
  } else if (x$type == 3) {
    'regression'
  }
}
predict_model.svm <- function(x, newdata) {
  type <- model_type(x)
  if (type == 'classification') {
    y_hat <- predict(x, newdata, probability = TRUE)
    y_hat <- attr(y_hat, 'probabilities')[, 1]
    y_hat <- qlogis(y_hat)
  } else {
    y_hat <- predict(x, newdata)
  }
  return(y_hat)
}

# Focus on borderline cases
x_i <- df %>%
  mutate(y_hat = y_hat) %>%
  filter(between(y_hat, 0.7, 0.74)) %>%
  select(-starts_with('y'))

# Define reference distribution
x <- select(df, -starts_with('y'))
explainer <- shapr(x, f, feature_labels = colnames(x))
phi_0 <- qlogis(0.7)

# Marginal Shapley values
msv <- explain(x_i, approach = 'empirical', type = 'independence',
               explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Conditional Shapley values
csv <- explain(x_i, approach = c(rep('empirical', 2), rep('gaussian', 5)),
               explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Interventional Shapley values
isv <- explain(x_i, approach = 'causal',
               explainer = explainer, prediction_zero = phi_0,
               ordering = list(c(1:3), c(4:7)))$dt %>%
  select(-none)

# Plot results
vals <- c('Marginal', 'Conditional', 'Interventional')
tmp <- rbind(msv, csv, isv) %>%
  mutate(reference = factor(rep(vals, each = nrow(msv)), levels = vals)) %>%
  pivot_longer(cols = -reference, names_to = 'feature', values_to = 'phi') %>%
  mutate(value = rep(as.numeric(t(x_i)), times = 3))%>%
  as.data.table(.)
tmp[!feature %in% c('gender', 'marital'), value := scale(value), by = feature]
fwrite(tmp, 'credit_classical.csv')
ggplot(tmp, aes(phi, feature, fill = value)) + 
  geom_jitter(size = 1.5, width = 0, height = 0.2, color = 'black', pch = 21,
              alpha = 0.75) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Feature\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ reference)
ggsave('credit_classical_shap.pdf', width = 10, height = 7)


### RATIONAL SHAPLEY ### 

# Define the subspace: top quartile of predictions
subspace <- df %>% 
  filter(y_hat >= quantile(y_hat, 0.75)) 
phi_0 <- qlogis(mean(subspace$y_hat))
x_ref <- subspace %>% 
  select(-starts_with('y'))
explainer <- shapr(x_ref, f, feature_labels = colnames(x))

# Marginal Shapley values
msv_r <- explain(x_i, approach = 'empirical', type = 'independence',
                 explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Conditional Shapley values
csv_r <- explain(x_i, approach = c(rep('empirical', 2), rep('gaussian', 5)),
                 explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Interventional Shapley values
isv_r <- explain(x_i, approach = 'causal',
                 explainer = explainer, prediction_zero = phi_0,
                 ordering = list(c(1:3), c(4:7)))$dt %>%
  select(-none)

# Plot results
tmp <- rbind(msv_r, csv_r, isv_r) %>%
  mutate(reference = factor(rep(vals, each = nrow(msv)), levels = vals)) %>%
  pivot_longer(cols = -reference, names_to = 'feature', values_to = 'phi') %>%
  mutate(value = rep(as.numeric(t(x_i)), times = 3)) %>%
  as.data.table(.)
tmp[!feature %in% c('gender', 'marital'), value := scale(value), by = feature]
fwrite(tmp, 'credit_rational.csv')
ggplot(tmp, aes(phi, feature, fill = value)) + 
  geom_jitter(size = 1.5, width = 0, height = 0.2, color = 'black', pch = 21,
              alpha = 0.75) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Feature\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ reference)
ggsave('credit_rational_shap.pdf', width = 10, height = 7)





