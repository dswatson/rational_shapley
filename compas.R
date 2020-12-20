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
  mutate(y = v_decile_score) %>%
  filter(race %in% c('African-American', 'Hispanic', 'Caucasian')) %>%
  mutate(sex = factor(sex), race = relevel(as.factor(race), ref = 'Caucasian'), 
         crim = factor(c_charge_degree)) %>%
  select(sex, age, race, priors_count, crim, y)
n <- nrow(df)

# Preprocess
x <- model.matrix(~ ., data = select(df, -y))
x <- x[, -1]
y <- df$y
f <- ranger(x = x, y = y)

# Focus on high risk defendants
x_i <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95))
x_i <- model.matrix(~ ., data = select(x_i, -y, -y_hat))
x_i <- x_i[, -1]

# Marginal Shapley values
p_wrap <- function(object, newdata) {
  predict(object, newdata)$predictions
}
msv <- fastshap::explain(f, X = x, nsim = 2000, pred_wrapper = p_wrap, 
                         newdata = x_i, adjust = TRUE)
race <- msv[[3]] + msv[[4]]
msv <- msv %>%
  as.data.frame(.) %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

# Conditional Shapley values
csv_ref <- ranger.unify(f, x)
csv <- treeshap(csv_ref, x_i)
race <- csv$shaps[[3]] + csv$shaps[[4]]
csv <- csv$shaps %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

# Interventional Shapley values
explainer <- shapr(x, f)
phi_0 <- mean(y)
isv <- shapr::explain(x_i, approach = 'causal', n_samples = 1000,
                      explainer = explainer, prediction_zero = phi_0,
                      ordering = list(c(1:2), c(3:6)))$dt %>%
  select(-none)
race <- isv[[3]] + isv[[4]]
isv <- isv %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

# Plot results
vals <- c('Marginal', 'Conditional', 'Interventional')
rbind(msv, csv, isv) %>%
  rename(sex = sexMale, priors = priors_count, crim = crimM) %>%
  mutate(value = factor(rep(vals, each = 290), levels = vals)) %>%
  pivot_longer(cols = -value, names_to = 'feature', values_to = 'phi') %>%
  ggplot(aes(phi, feature, fill = phi)) + 
  geom_jitter(size = 2, alpha = 0.75, width = 0, height = 0.1, 
              color = 'black', pch = 21) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Shapley\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature', title = 'Classical Shapley Values') 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ value)
ggsave('compas_classical_shap.pdf', width = 10, height = 7)


### RATIONAL SHAPLEY VALUES ###

# Define the subspace
subspace <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(age <= 28, priors_count >= 1, y_hat <= median(y_hat)) 
x_ref <- model.matrix(~ ., data = select(subspace, -y, -y_hat))
x_ref <- x_ref[, -1]

# Marginal Shapley values
msv_r <- fastshap::explain(f, X = x_ref, nsim = 2000, pred_wrapper = p_wrap, 
                           newdata = x_i, adjust = TRUE)
race <- msv_r[[3]] + msv_r[[4]]
msv_r <- msv_r %>%
  as.data.frame(.) %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

# Treeshap isn't working with this subspace for some reason...switching to 
# shapr::approach = 'empirical'
explainer <- shapr(x_ref, f)
phi_0 <- mean(subspace$y)

# Conditional Shapley values
csv_r <- shapr::explain(x_i, approach = 'empirical', n_samples = 1000,
                        explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)
race <- csv_r[[3]] + csv_r[[4]]
csv_r <- csv_r %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

# Interventional Shapley values
isv_r <- shapr::explain(x_i, approach = 'causal', n_samples = 1000,
                        explainer = explainer, prediction_zero = phi_0,
                        ordering = list(c(1:2), c(3:6)))$dt %>%
  select(-none)
race <- isv_r[[3]] + isv_r[[4]]
isv_r <- isv_r %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

# Plot
rbind(msv_r, csv_r, isv_r) %>%
  rename(sex = sexMale, priors = priors_count, crim = crimM) %>%
  mutate(value = factor(rep(vals, each = 290), levels = vals)) %>%
  pivot_longer(cols = -value, names_to = 'feature', values_to = 'phi') %>%
  ggplot(aes(phi, feature, fill = phi)) + 
  geom_jitter(size = 2, alpha = 0.75, width = 0, height = 0.1, 
              color = 'black', pch = 21) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Shapley\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature', title = 'Rational Shapley Values') 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ value)
ggsave('compas_rational_shap.pdf', width = 10, height = 7)

### Boxplots ###

# Classical
dat <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95)) %>%
  mutate(Marginal = msv$race, Conditional = csv$race, 
         Interventional = isv$race) %>%
  select(race, Marginal, Conditional, Interventional) %>%
  rename(Race = race) %>%
  pivot_longer(cols = Marginal:Interventional, names_to = 'value', values_to = 'phi') %>%
  mutate(value = factor(value, levels = c('Marginal', 'Conditional', 'Interventional')),
         Race = factor(Race, levels = c('African-American', 'Caucasian', 'Hispanic')))
ggplot(dat, aes(Race, phi, color = Race)) + 
  geom_boxplot() + 
  geom_jitter(size = 1.5, alpha = 0.75) +
  scale_color_npg() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Classical Shapley Values', y = 'Shapley Value') +
  facet_wrap(~ value)
ggsave('compas_classical_boxplot.pdf', width = 10, height = 7)

# Rational
dat <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95)) %>%
  mutate(Marginal = msv_r$race, Conditional = csv_r$race, 
         Interventional = isv_r$race) %>%
  select(race, Marginal, Conditional, Interventional) %>%
  rename(Race = race) %>%
  pivot_longer(cols = Marginal:Interventional, names_to = 'value', values_to = 'phi') %>%
  mutate(value = factor(value, levels = c('Marginal', 'Conditional', 'Interventional')),
         Race = factor(Race, levels = c('African-American', 'Caucasian', 'Hispanic')))
ggplot(dat, aes(Race, phi, color = Race)) + 
  geom_boxplot() + 
  geom_jitter(size = 1.5, alpha = 0.75) +
  scale_color_npg() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Rational Shapley Values', y = 'Shapley Value') +
  facet_wrap(~ value)
ggsave('compas_rational_boxplot.pdf', width = 10, height = 7)




