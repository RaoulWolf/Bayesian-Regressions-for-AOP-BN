library(brms)
options(mc.cores = parallel::detectCores())

data <- readRDS("data/random_data.rds")

dcp_to_oxphos_data <- na.omit(data[, c("conc", "oxphos")])
dcp_to_oxphos_data$conc <- ifelse(dcp_to_oxphos_data$conc == 0, .Machine$double.eps, dcp_to_oxphos_data$conc)

## Model

dcp_to_oxphos_formula <- brmsformula(formula = oxphos ~ lower + (upper - lower) / (1 + exp(-slope * log(conc / ec50))),
                                     flist = list(slope ~ 1, lower ~ 1, upper ~ 1, ec50 ~ 1),
                                     family = poisson(link = "identity"),
                                     nl = TRUE)

get_prior(formula = dcp_to_oxphos_formula,
          data = dcp_to_oxphos_data)

plot(oxphos ~ conc, data = dcp_to_oxphos_data)

dcp_to_oxphos_priors <- c(prior(prior = normal(    1,    1), class = "b", nlpar = "ec50",  lb = 0),
                          prior(prior = normal( 3000, 1000), class = "b", nlpar = "lower", lb = 0),
                          prior(prior = normal(  -10,   10), class = "b", nlpar = "slope", ub = 0),
                          prior(prior = normal(10000, 1000), class = "b", nlpar = "upper", lb = 0))

dcp_to_oxphos <- brm(formula = dcp_to_oxphos_formula,
                     data = dcp_to_oxphos_data,
                     prior = dcp_to_oxphos_priors,
                     chains = 5,
                     iter = 4000,
                     file = "models/dcp_to_oxphos")

summary(dcp_to_oxphos)

plot(dcp_to_oxphos)

pp_check(dcp_to_oxphos, type = "dens_overlay")

pp_check(dcp_to_oxphos, type = "ecdf_overlay")

conditional_effects(dcp_to_oxphos)

fixef(dcp_to_oxphos)

## Simulations

dcp_to_oxphos_pred <- data.frame(conc = seq(2/3 * min(data$conc, na.rm = TRUE), 1.5 * max(data$conc, na.rm = TRUE), length.out = 10000))

dcp_to_oxphos_pred$oxphos_pred <- as.vector(posterior_predict(dcp_to_oxphos, newdata = dcp_to_oxphos_pred, nsamples = 1))

head(dcp_to_oxphos_pred)

plot(oxphos_pred ~ conc, data = dcp_to_oxphos_pred, type = "l")
points(oxphos ~ conc, data = dcp_to_oxphos_data, col = "red")

write.csv(dcp_to_oxphos_pred, "simulations/dcp_to_oxphos_pred.csv")
saveRDS(dcp_to_oxphos_pred, "simulations/dcp_to_oxphos_pred.rds")
