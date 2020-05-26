library(brms)
options(mc.cores = parallel::detectCores() - 1L)

data <- readRDS("data/random_data.rds")

ros_to_lpo_data <- na.omit(data[, c("ros", "lpo")])
ros_to_lpo_data$ros <- ifelse(ros_to_lpo_data$ros == 0L, 1L, ros_to_lpo_data$ros)

## Model

ros_to_lpo_formula <- brmsformula(formula = lpo | trunc(lb = 0) ~ lower + (upper - lower) / (1 + exp(-slope * log(ros / ec50))),
                                  flist = list(slope ~ 1, lower ~ 1, upper ~ 1, ec50 ~ 1),
                                  family = gaussian(link = "identity"),
                                  nl = TRUE)

get_prior(formula = ros_to_lpo_formula,
          data = ros_to_lpo_data)

plot(lpo ~ ros, data = ros_to_lpo_data)

ros_to_lpo_priors <- c(prior(prior = normal(6000, 1000), class = "b", nlpar = "ec50",  lb = 0),
                       prior(prior = normal(   2,    1), class = "b", nlpar = "lower", lb = 0),
                       prior(prior = normal(  10,    5), class = "b", nlpar = "slope", lb = 0),
                       prior(prior = normal(   5,    1), class = "b", nlpar = "upper", lb = 0))

ros_to_lpo <- brm(formula = ros_to_lpo_formula,
                  data = ros_to_lpo_data,
                  prior = ros_to_lpo_priors,
                  chains = 5,
                  iter = 4000,
                  file = "models/ros_to_lpo")

summary(ros_to_lpo)

plot(ros_to_lpo)

pp_check(ros_to_lpo, type = "dens_overlay")

pp_check(ros_to_lpo, type = "ecdf_overlay")

conditional_effects(ros_to_lpo)

fixef(ros_to_lpo)

## Simulations

ros_to_lpo_pred <- data.frame(ros = as.integer(seq(2/3 * min(data$ros, na.rm = TRUE), 1.5 * max(data$ros, na.rm = TRUE), length.out = 10000)))

ros_to_lpo_pred$lpo_pred <- as.vector(posterior_predict(ros_to_lpo, newdata = ros_to_lpo_pred, nsamples = 1))

head(ros_to_lpo_pred)

plot(lpo_pred ~ ros, data = ros_to_lpo_pred, type = "l")
points(lpo ~ ros, data = ros_to_lpo_data, col = "red")

write.csv(ros_to_lpo_pred, "simulations/ros_to_lpo_pred.csv")
saveRDS(ros_to_lpo_pred, "simulations/ros_to_lpo_pred.rds")
