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
  select(-starts_with('y')) %>%
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
isv <- explain(x_i, approach = 'causal', 
               explainer = explainer, prediction_zero = phi_0,
               ordering = list(c(1:3), c(4:7), c(8:10)))$dt %>%
  select(-none)

# Tidy
vals <- c('Marginal', 'Conditional', 'Interventional')
res_c <- rbind(msv, csv, isv) %>%
  mutate(reference = factor(rep(vals, each = nrow(msv)), levels = vals)) %>%
  pivot_longer(cols = -reference, names_to = 'feature', values_to = 'phi') %>%
  mutate(value = rep(as.numeric(t(x_i)), times = 3), type = 'Classical')

### RATIONAL SHAPLEY ### 

# Define the subspace: same sex and age
subspace <- df %>% 
  filter(sex > 0, age > quantile(age, 0.25), y < median(y)) %>%
  mutate(sex = sex + rnorm(76, sd = 0.001)) # Crashes when a feature is totally invariant
x_ref <- subspace %>% 
  select(-y) %>%
  as.matrix(.)
explainer <- shapr(x_ref, f)
phi_0 <- mean(subspace$y)

# Marginal Shapley values
msv_r <- explain(x_i, approach = 'empirical', type = 'independence',
                 explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Conditional Shapley values
csv_r <- explain(x_i, approach = c('empirical', rep('gaussian', 9)),
                 explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Interventional Shapley values
isv_r <- explain(x_i, approach = 'causal', 
                 explainer = explainer, prediction_zero = phi_0,
                 ordering = list(c(1:3), c(4:7), c(8:10)))$dt %>%
  select(-none)

# Tidy, export
res_r <- rbind(msv_r, csv_r, isv_r) %>%
  mutate(reference = factor(rep(vals, each = nrow(msv)), levels = vals)) %>%
  pivot_longer(cols = -reference, names_to = 'feature', values_to = 'phi') %>%
  mutate(value = rep(as.numeric(t(x_i)), times = 3), type = 'Rational')
res <- rbind(res_c, res_r)
fwrite(res, 'medical_res.csv')

# Plot 
ggplot(res, aes(phi, feature, fill = value)) + 
  geom_jitter(size = 1.5, width = 0, height = 0.2, color = 'black', pch = 21,
              alpha = 0.75) + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed') +
  scale_fill_gradient2('Feature\nValue', low = 'blue', high = 'red') +
  labs(x = 'Shapley Value', y = 'Feature') +
  theme_bw() + 
  facet_grid(type ~ reference)
ggsave('medical_shap.pdf', width = 10, height = 7)


# The case of Bert & Ernie
x_i <- x_i %>% 
  as.data.frame(.) %>% 
  filter(bmi > 0.05, map > 0.05) %>% 
  mutate(idx = row_number()) %>% 
  filter(idx == 3) %>% 
  select(-idx) 

# Bert
b <- df %>%
  filter(bmi <= quantile(bmi, 0.2))
explainer <- shapr(select(b, -y), f)
phi_0 <- mean(b$y)
csv_b <- explain(x_i, approach = c('empirical', rep('gaussian', 9)),
                 explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)

# Ernie
e <- df %>% 
  filter(map <= quantile(map, 0.2))
explainer <- shapr(select(e, -y), f)
phi_0 <- mean(e$y)
csv_e <- explain(x_i, approach = c('empirical', rep('gaussian', 9)),
                 explainer = explainer, prediction_zero = phi_0)$dt %>%
  select(-none)




