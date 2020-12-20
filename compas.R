# Set working directory
setwd('~/Documents/rational_shapley')

# Load libraries
library(data.table)
library(tidyverse)
library(ranger)
library(shapr) 
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
x <- df %>% select(-y)
y <- df$y
f <- ranger(x = x, y = y, num.trees = 1000)

# Focus on high risk defendants
x_i <- x %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95))

# Prep
explainer <- shapr(x, f)
phi_0 <- mean(y)







# Marginal Shapley values
msv <- explain(select(x_i, -y_hat), approach = 'empirical', type = 'independence',
               explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Conditional Shapley values
csv <- explain(select(x_i, -y_hat), approach = 'ctree', explainer = explainer, 
               prediction_zero = phi_0)$dt %>%
  select(-none)

# Interventional Shapley values






# Conditional Shapley values
csv_ref <- ranger.unify(f, x)
csv <- treeshap(csv_ref, x_i)
race <- csv$shaps[[3]] + csv$shaps[[4]]
csv <- csv$shaps %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

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

# Plot results
rbind(msv, csv) %>%
  rename(sex = sexMale, priors = priors_count, crim = crimM) %>%
  mutate(value = rep(c('Marginal', 'Conditional'), each = nrow(msv))) %>%
  mutate(value = factor(value, levels = c('Marginal', 'Conditional'))) %>%
  pivot_longer(cols = -value, names_to = 'feature', values_to = 'phi') %>%
  ggplot(aes(phi, feature, fill = phi)) + 
  geom_jitter(size = 2, width = 0, height = 0.1, color = 'black', pch = 21) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Shapley\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature') + 
  theme_bw() + 
  facet_wrap(~ value)

# Boxplot of Shapley values by race
dat <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(y_hat >= quantile(y_hat, 0.95)) %>%
  mutate(Marginal = msv[, 5], Conditional = csv[, 5]) %>%
  select(race, Marginal, Conditional) %>%
  rename(Race = race) %>%
  pivot_longer(cols = Marginal:Conditional, names_to = 'value', values_to = 'phi') %>%
  mutate(value = factor(value, levels = c('Marginal', 'Conditional')),
         Race = factor(Race, levels = c('African-American', 'Caucasian', 'Hispanic')))
ggplot(dat, aes(Race, phi, color = Race)) + 
  geom_boxplot() + 
  geom_jitter(size = 1.5, alpha = 0.75) +
  scale_color_d3() +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme_bw() + 
  labs(y = 'Shapley Value') +
  facet_wrap(~ value)

### RATIONAL SHAPLEY ### 

# Define the subspace
subspace <- df %>%
  mutate(y_hat = f$predictions) %>%
  filter(age <= 28, priors_count >= 1, y_hat <= median(y_hat)) 
x_ref <- model.matrix(~ ., data = subspace)
x_ref <- x_ref[, 2:7]

# Conditional value function
csv_ref <- ranger.unify(f, x_ref)
csv_r <- treeshap(csv_ref, x_i)
race <- csv_r$shaps[[3]] + csv_r$shaps[[4]]
csv_r <- csv_r$shaps %>%
  select(-starts_with('race')) %>%
  mutate(race = race)

# Marginal value function
msv_r <- fastshap::explain(f, X = x_ref, nsim = 2000, pred_wrapper = p_wrap, 
                           newdata = x_i, adjust = TRUE)
race <- msv_r[[3]] + msv_r[[4]]
msv_r <- msv_r %>%
  as.data.frame(.) %>%
  select(-starts_with('race')) %>%
  mutate(race = race)








rbind(msv_r, csv_r) %>%
  rename(sex = sexMale, priors = priors_count, crim = crimM) %>%
  mutate(value = rep(c('Marginal', 'Conditional'), each = nrow(msv))) %>%
  mutate(value = factor(value, levels = c('Marginal', 'Conditional'))) %>%
  pivot_longer(cols = -value, names_to = 'feature', values_to = 'phi') %>%
  ggplot(aes(phi, feature, fill = phi)) + 
  geom_jitter(size = 2, width = 0, height = 0.1, color = 'black', pch = 21) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_viridis_c('Shapley\nValue', option = 'B') +
  labs(x = 'Shapley Value', y = 'Feature') + 
  theme_bw() + 
  facet_wrap(~ value)



# Scatter: admissible
tmp <- data.frame(
  Classical = abs(c(msv$age, msv$priors_count, csv$age, csv$priors_count)),
  Rational = abs(c(msv_r$age, msv_r$priors_count, csv_r$age, csv_r$priors_count)),
  Reference = rep(c('Marginal', 'Conditional'), each = 20),
  Feature = rep(c('Age', 'Priors'), each = 10)
)
tmp <- tmp %>% 
  mutate(Reference = factor(Reference, levels = c('Marginal', 'Conditional')))
ggplot(tmp, aes(Classical, Rational)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') + 
  facet_wrap(Reference ~ Feature, ncol = 2) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave('Conditioning.pdf', width = 6, height = 6)

# Binomial test
x <- ifelse(tmp$Classical - tmp$Rational > 0, 0, 1)
binom.test(x, length(x), alt = 'less')

# Scatter: inadmissible
# https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
tmp <- data.frame(
  Classical = c(msv$sex, msv$race, csv$sex, csv$race),
  Rational = c(msv_r$sex, msv_r$race, csv_r$sex, csv_r$race),
  Reference = rep(c('Marginal', 'Conditional'), each = 20),
  Feature = rep(c('Sex', 'Race'), each = 10)
)
tmp <- tmp %>% 
  mutate(Reference = factor(Reference, levels = c('Marginal', 'Conditional')))
ggplot(tmp, aes(Classical, Rational)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') + 
  facet_wrap_equal(Reference ~ Feature, ncol = 2, scales = 'free') + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave('Inadmissible.pdf', width = 6, height = 6)










