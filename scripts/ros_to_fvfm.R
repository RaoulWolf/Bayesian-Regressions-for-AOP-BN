library(brms)
options(mc.cores = parallel::detectCores() - 1L)

data <- readRDS("data/random_data.rds")

ros_to_fv_fm_data <- na.omit(data[, c("ros", "fv_fm")])
ros_to_fv_fm_data$ros <- ifelse(ros_to_fv_fm_data$ros == 0L, 1L, ros_to_fv_fm_data$ros)

## Model

ros_to_fv_fm_formula <- brmsformula(formula = fv_fm ~ lower + (upper - lower) / (1 + exp(-slope * log(ros / ec50))),
                                     flist = list(slope ~ 1, lower ~ 1, upper ~ 1, ec50 ~ 1),
                                     family = Beta(link = "identity"),
                                     nl = TRUE)

get_prior(formula = ros_to_fv_fm_formula,
          data = ros_to_fv_fm_data)

plot(fv_fm ~ ros, data = ros_to_fv_fm_data)

ros_to_fv_fm_priors <- c(prior(prior = normal(7000,   1000),   class = "b", nlpar = "ec50",  lb = 0),
                         prior(prior = normal(   0.1,    0.1), class = "b", nlpar = "lower", lb = 0, ub = 1),
                         prior(prior = normal(  -5,      5),   class = "b", nlpar = "slope", ub = 0),
                         prior(prior = normal(   0.7,    0.1), class = "b", nlpar = "upper", lb = 0, ub = 1))

ros_to_fv_fm <- brm(formula = ros_to_fv_fm_formula,
                     data = ros_to_fv_fm_data,
                     prior = ros_to_fv_fm_priors,
                     chains = 5,
                     iter = 4000,
                     file = "models/ros_to_fv_fm")

summary(ros_to_fv_fm)

plot(ros_to_fv_fm)

pp_check(ros_to_fv_fm, type = "dens_overlay")

pp_check(ros_to_fv_fm, type = "ecdf_overlay")

conditional_effects(ros_to_fv_fm)

fixef(ros_to_fv_fm)

## Simulations

ros_to_fv_fm_pred <- data.frame(ros = as.integer(seq(2/3 * min(data$ros, na.rm = TRUE), 1.5 * max(data$ros, na.rm = TRUE), length.out = 10000)))

ros_to_fv_fm_pred$fv_fm_pred <- as.vector(posterior_predict(ros_to_fv_fm, newdata = ros_to_fv_fm_pred, nsamples = 1))

head(ros_to_fv_fm_pred)

plot(fv_fm_pred ~ ros, data = ros_to_fv_fm_pred, type = "l")
points(fv_fm ~ ros, data = ros_to_fv_fm_data, col = "red")

write.csv(ros_to_fv_fm_pred, "simulations/ros_to_fv_fm_pred.csv")
saveRDS(ros_to_fv_fm_pred, "simulations/ros_to_fv_fm_pred.rds")
