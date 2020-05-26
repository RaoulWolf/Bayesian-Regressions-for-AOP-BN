library(brms)
options(mc.cores = parallel::detectCores() - 1L)

data <- readRDS("data/random_data.rds")

fv_fm_to_fronds_data <- na.omit(data[, c("fv_fm", "fronds_number")])
fv_fm_to_fronds_data$fv_fm <- ifelse(fv_fm_to_fronds_data$fv_fm == 0, .Machine$double.eps, fv_fm_to_fronds_data$fv_fm)

## Model

fv_fm_to_fronds_formula <- brmsformula(formula = fronds_number ~ lower + (upper - lower) * (1 - exp(-fv_fm / slope)),
                                       flist = list(slope ~ 1, lower ~ 1, upper ~ 1),
                                       family = poisson(link = "identity"),
                                       nl = TRUE)

get_prior(formula = fv_fm_to_fronds_formula,
          data = fv_fm_to_fronds_data)

plot(fronds_number ~ fv_fm, data = fv_fm_to_fronds_data)

fv_fm_to_fronds_priors <- c(prior(prior = normal( 10, 10), class = "b", nlpar = "lower", lb = 0),
                            prior(prior = normal(  1,  1), class = "b", nlpar = "slope", lb = 0),
                            prior(prior = normal(110, 20), class = "b", nlpar = "upper", lb = 0))

fv_fm_to_fronds <- brm(formula = fv_fm_to_fronds_formula,
                       data = fv_fm_to_fronds_data,
                       prior = fv_fm_to_fronds_priors,
                       chains = 5,
                       iter = 4000,
                       file = "models/fv_fm_to_fronds")

summary(fv_fm_to_fronds)

plot(fv_fm_to_fronds)

pp_check(fv_fm_to_fronds, type = "dens_overlay")

pp_check(fv_fm_to_fronds, type = "ecdf_overlay")

conditional_effects(fv_fm_to_fronds)

fixef(fv_fm_to_fronds)

## Simulations

fv_fm_to_fronds_pred <- data.frame(fv_fm = seq(2/3 * min(data$fv_fm, na.rm = TRUE), 1.5 * max(data$fv_fm, na.rm = TRUE), length.out = 10000))

fv_fm_to_fronds_pred$fronds_number_pred <- as.vector(posterior_predict(fv_fm_to_fronds, newdata = fv_fm_to_fronds_pred, nsamples = 1))

head(fv_fm_to_fronds_pred)

plot(fronds_number_pred ~ fv_fm, data = fv_fm_to_fronds_pred, type = "l")
points(fronds_number ~ fv_fm, data = fv_fm_to_fronds_data, col = "red")

write.csv(fv_fm_to_fronds_pred, "simulations/fv_fm_to_fronds_pred.csv")
saveRDS(fv_fm_to_fronds_pred, "simulations/fv_fm_to_fronds_pred.rds")
