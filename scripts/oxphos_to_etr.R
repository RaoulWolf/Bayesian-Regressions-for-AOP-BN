library(brms)
options(mc.cores = parallel::detectCores() - 1L)

data <- readRDS("data/random_data.rds")

oxphos_to_etr_data <- na.omit(data[, c("oxphos", "etr")])
oxphos_to_etr_data$oxphos <- ifelse(oxphos_to_etr_data$oxphos == 0L, 1L, oxphos_to_etr_data$oxphos)

## Model

oxphos_to_etr_formula <- brmsformula(formula = etr | trunc(lb = 0) ~ lower + (upper - lower) / (1 + exp(-slope * log(oxphos / ec50))),
                                     flist = list(slope ~ 1, lower ~ 1, upper ~ 1, ec50 ~ 1),
                                     family = gaussian(link = "identity"),
                                     nl = TRUE)

get_prior(formula = oxphos_to_etr_formula,
          data = oxphos_to_etr_data)

plot(etr ~ oxphos, data = oxphos_to_etr_data)

oxphos_to_etr_priors <- c(prior(prior = normal(6000, 1000), class = "b", nlpar = "ec50",  lb = 0),
                          prior(prior = normal(   1,    1), class = "b", nlpar = "lower", lb = 0),
                          prior(prior = normal(   5,    2), class = "b", nlpar = "slope", lb = 0),
                          prior(prior = normal(  20,    5), class = "b", nlpar = "upper", lb = 0))

oxphos_to_etr <- brm(formula = oxphos_to_etr_formula,
                     data = oxphos_to_etr_data,
                     prior = oxphos_to_etr_priors,
                     chains = 5,
                     iter = 4000,
                     file = "models/oxphos_to_etr")

summary(oxphos_to_etr)

plot(oxphos_to_etr)

pp_check(oxphos_to_etr, type = "dens_overlay")

pp_check(oxphos_to_etr, type = "ecdf_overlay")

conditional_effects(oxphos_to_etr)

fixef(oxphos_to_etr)

## Simulations

oxphos_to_etr_pred <- data.frame(oxphos = as.integer(seq(2/3 * min(data$oxphos, na.rm = TRUE), 1.5 * max(data$oxphos, na.rm = TRUE), length.out = 10000)))

oxphos_to_etr_pred$etr_pred <- as.vector(posterior_predict(oxphos_to_etr, newdata = oxphos_to_etr_pred, nsamples = 1))

head(oxphos_to_etr_pred)

plot(etr_pred ~ oxphos, data = oxphos_to_etr_pred, type = "l")
points(etr ~ oxphos, data = oxphos_to_etr_data, col = "red")

write.csv(oxphos_to_etr_pred, "simulations/oxphos_to_etr_pred.csv")
saveRDS(oxphos_to_etr_pred, "simulations/oxphos_to_etr_pred.rds")
