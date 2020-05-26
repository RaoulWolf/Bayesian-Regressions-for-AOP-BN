library(brms)
options(mc.cores = parallel::detectCores() - 1L)

data <- readRDS("data/random_data.rds")

etr_to_fronds_data <- na.omit(data[, c("etr", "fronds_number")])
etr_to_fronds_data$etr <- ifelse(etr_to_fronds_data$etr == 0, .Machine$double.eps, etr_to_fronds_data$etr)

## Model

etr_to_fronds_formula <- brmsformula(formula = fronds_number ~ lower + (upper - lower) * (1 - exp(-etr / slope)),
                                     flist = list(slope ~ 1, lower ~ 1, upper ~ 1),
                                     family = poisson(link = "identity"),
                                     nl = TRUE)

get_prior(formula = etr_to_fronds_formula,
          data = etr_to_fronds_data)

plot(fronds_number ~ etr, data = etr_to_fronds_data)

etr_to_fronds_priors <- c(prior(prior = normal( 10, 10), class = "b", nlpar = "lower", lb = 0),
                          prior(prior = normal(  5,  5), class = "b", nlpar = "slope", lb = 0),
                          prior(prior = normal(110, 20), class = "b", nlpar = "upper", lb = 0))

etr_to_fronds <- brm(formula = etr_to_fronds_formula,
                     data = etr_to_fronds_data,
                     prior = etr_to_fronds_priors,
                     chains = 5,
                     iter = 4000,
                     file = "models/etr_to_fronds")

summary(etr_to_fronds)

plot(etr_to_fronds)

pp_check(etr_to_fronds, type = "dens_overlay")

pp_check(etr_to_fronds, type = "ecdf_overlay")

conditional_effects(etr_to_fronds)

fixef(etr_to_fronds)

## Simulations

etr_to_fronds_pred <- data.frame(etr = seq(2/3 * min(data$etr, na.rm = TRUE), 1.5 * max(data$etr, na.rm = TRUE), length.out = 10000))

etr_to_fronds_pred$fronds_number_pred <- as.vector(posterior_predict(etr_to_fronds, newdata = etr_to_fronds_pred, nsamples = 1))

head(etr_to_fronds_pred)

plot(fronds_number_pred ~ etr, data = etr_to_fronds_pred, type = "l")
points(fronds_number ~ etr, data = etr_to_fronds_data, col = "red")

write.csv(etr_to_fronds_pred, "simulations/etr_to_fronds_pred.csv")
saveRDS(etr_to_fronds_pred, "simulations/etr_to_fronds_pred.rds")
