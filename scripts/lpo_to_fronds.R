library(brms)
options(mc.cores = parallel::detectCores() - 1L)

data <- readRDS("data/random_data.rds")

lpo_to_fronds_data <- na.omit(data[, c("lpo", "fronds_number")])
lpo_to_fronds_data$lpo <- ifelse(lpo_to_fronds_data$lpo == 0, .Machine$double.eps, lpo_to_fronds_data$lpo)

## Model

lpo_to_fronds_formula <- brmsformula(formula = fronds_number ~ lower + (upper - lower) / (1 + exp(-slope * log(lpo / ec50))),
                                     flist = list(slope ~ 1, lower ~ 1, upper ~ 1, ec50 ~ 1),
                                     family = poisson(link = "identity"),
                                     nl = TRUE)

get_prior(formula = lpo_to_fronds_formula,
          data = lpo_to_fronds_data)

plot(fronds_number ~ lpo, data = lpo_to_fronds_data)

lpo_to_fronds_priors <- c(prior(prior = normal(  3,  1), class = "b", nlpar = "ec50",  lb = 0),
                          prior(prior = normal( 20,  5), class = "b", nlpar = "lower", lb = 0),
                          prior(prior = normal(-10,  2), class = "b", nlpar = "slope", ub = 0),
                          prior(prior = normal(120, 10), class = "b", nlpar = "upper", lb = 0))

lpo_to_fronds <- brm(formula = lpo_to_fronds_formula,
                     data = lpo_to_fronds_data,
                     prior = lpo_to_fronds_priors,
                     chains = 5,
                     iter = 4000,
                     file = "models/lpo_to_fronds")

summary(lpo_to_fronds)

plot(lpo_to_fronds)

pp_check(lpo_to_fronds, type = "dens_overlay")

pp_check(lpo_to_fronds, type = "ecdf_overlay")

conditional_effects(lpo_to_fronds)

fixef(lpo_to_fronds)

## Simulations

lpo_to_fronds_pred <- data.frame(lpo = seq(2/3 * min(data$lpo, na.rm = TRUE), 1.5 * max(data$lpo, na.rm = TRUE), length.out = 10000))

lpo_to_fronds_pred$fronds_number_pred <- as.vector(posterior_predict(lpo_to_fronds, newdata = lpo_to_fronds_pred, nsamples = 1))

head(lpo_to_fronds_pred)

plot(fronds_number_pred ~ lpo, data = lpo_to_fronds_pred, type = "l")
points(fronds_number ~ lpo, data = lpo_to_fronds_data, col = "red")

write.csv(lpo_to_fronds_pred, "simulations/lpo_to_fronds_pred.csv")
saveRDS(lpo_to_fronds_pred, "simulations/lpo_to_fronds_pred.rds")
