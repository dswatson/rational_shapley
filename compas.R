# Set working directory
setwd('~/Documents/rational_shapley')

# Load libraries
library(data.table)
library(tidyverse)
library(ranger)
library(fastshap)
library(treeshap)
library(shapr) 
library(ggsci)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import data
df <- read_csv('compas-scores-two-years.csv')

# From ProPublica GitHub page
df <- df %>%
  filter(days_b_screening_arrest <= 30) %>%
  filter(days_b_screening_arrest >= -30) %>%
  filter(is_recid != -1) %>%
  filter(c_charge_degree != 'O') %>%
  filter(score_text != 'N/A') %>%
  # My riffing 
  mutate(y = v_decile_score, priors = log(priors_count + 1), age = log(age)) %>%
  filter(race %in% c('African-American', 'Caucasian')) %>%
  as.data.table(.)
df[, race := ifelse(race == 'African-American', 1, 0)]
df[, sex := ifelse(sex == 'Male', 1, 0)]
df[, felony := ifelse(c_charge_degree == 'F', 1, 0)]
df <- df %>% select(sex, race, age, felony, priors, y)
n <- nrow(df)

# Preprocess
f <- ranger(y ~ ., data = df)

# Focus on high risk defendants
x <- df %>% select(-y)
x_i <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95)) %>%
  select(-starts_with('y'))

# Marginal Shapley values
p_wrap <- function(object, newdata) {
  predict(object, newdata)$predictions
}
msv <- fastshap::explain(f, X = x, nsim = 2000, pred_wrapper = p_wrap, 
                         newdata = x_i, adjust = TRUE) %>%
  as.data.table(.)

# Conditional Shapley values
csv_ref <- ranger.unify(f, x)
csv <- treeshap(csv_ref, x_i)$shaps

# Interventional Shapley values
explainer <- shapr(x, f)
phi_0 <- mean(df$y)
isv <- shapr::explain(x_i, approach = 'causal', 
                      explainer = explainer, prediction_zero = phi_0,
                      ordering = list(c(1:3), c(4:5)))$dt %>%
  select(-none)

# Plot results
vals <- c('Marginal', 'Conditional', 'Interventional')
tmp <- rbind(msv, csv, isv) %>%
  mutate(reference = factor(rep(vals, each = nrow(msv)), levels = vals)) %>%
  pivot_longer(cols = -reference, names_to = 'feature', values_to = 'phi') %>%
  mutate(value = rep(as.numeric(t(x_i)), times = 3)) %>%
  as.data.table(.)
tmp[feature %in% c('age', 'priors'), value := scale(value), by = feature]
fwrite(tmp, 'compas_classical.csv')
ggplot(tmp, aes(phi, feature, fill = value)) + 
  geom_jitter(size = 1.5, width = 0, height = 0.2, color = 'black', pch = 21,
              alpha = 0.75) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Feature\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature') +
  theme_bw() + 
  facet_wrap(~ reference)
ggsave('compas_classical_shap.pdf', width = 10, height = 7)


### RATIONAL SHAPLEY VALUES ###

# Define the subspace
subspace <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(age <= max(x_i$age), priors >= log(2), y_hat <= median(y_hat)) 
x_ref <- subspace %>% select(-starts_with('y'))

# Marginal Shapley values
msv_r <- fastshap::explain(f, X = x_ref, nsim = 2000, pred_wrapper = p_wrap, 
                           newdata = x_i, adjust = TRUE) %>%
  as.data.table(.)

# Conditional Shapley values (treeshap not working for some reason?)
approaches <- c(rep('empirical', 2), 'gaussian', 'empirical', 'gaussian')
csv_r <- shapr::explain(x_i, approach = approaches, 
                        explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Interventional Shapley values
isv_r <- shapr::explain(x_i, approach = 'causal', 
                        explainer = explainer, prediction_zero = phi_0,
                        ordering = list(c(1:3), c(4:5)))$dt %>%
  select(-none)

# Plot
tmp <- rbind(msv_r, csv_r, isv_r) %>%
  mutate(reference = factor(rep(vals, each = nrow(msv)), levels = vals)) %>%
  pivot_longer(cols = -reference, names_to = 'feature', values_to = 'phi') %>%
  mutate(value = rep(as.numeric(t(x_i)), times = 3)) %>%
  as.data.table(.)
tmp[feature %in% c('age', 'priors'), value := scale(value), by = feature]
fwrite(tmp, 'compas_rational.csv')
ggplot(tmp, aes(phi, feature, fill = value)) + 
  geom_jitter(size = 1.5, width = 0, height = 0.2, color = 'black', pch = 21,
              alpha = 0.75) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Feature\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature') +
  theme_bw() + 
  facet_wrap(~ reference)
ggsave('compas_rational_shap.pdf', width = 10, height = 7)


### Boxplots ###

# Classical
dat <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95)) %>%
  mutate(Marginal = msv$race, Conditional = csv$race, 
         Interventional = isv$race) %>%
  select(race, Marginal, Conditional, Interventional) %>%
  pivot_longer(cols = Marginal:Interventional, names_to = 'value', values_to = 'phi') %>%
  mutate(value = factor(value, levels = c('Marginal', 'Conditional', 'Interventional')),
         Race = ifelse(race == 1, 'African American', 'Caucasian')) %>%
  select(-race) %>%
  as.data.table(.)
ggplot(dat, aes(Race, phi, fill = Race)) + 
  geom_boxplot() + 
  geom_jitter(size = 1.5, alpha = 0.75, color = 'black', pch = 21) +
  scale_fill_d3() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  ylab('Shapley Value') +
  facet_wrap(~ value)
ggsave('compas_classical_boxplot.pdf', width = 10, height = 7)

# Rational
dat_r <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95)) %>%
  mutate(Marginal = msv_r$race, Conditional = csv_r$race, 
         Interventional = isv_r$race) %>%
  select(race, Marginal, Conditional, Interventional) %>%
  pivot_longer(cols = Marginal:Interventional, names_to = 'value', values_to = 'phi') %>%
  mutate(value = factor(value, levels = c('Marginal', 'Conditional', 'Interventional')),
         Race = ifelse(race == 1, 'African American', 'Caucasian')) %>% 
  select(-race) %>%
  as.data.table(.)
ggplot(dat_r, aes(Race, phi, fill = Race)) + 
  geom_boxplot() + 
  geom_jitter(size = 1.5, alpha = 0.75, color = 'black', pch = 21) +
  scale_fill_d3() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  ylab('Shapley Value') +
  facet_wrap(~ value)
ggsave('compas_rational_boxplot.pdf', width = 10, height = 7)




