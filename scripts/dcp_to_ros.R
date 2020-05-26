library(brms)
options(mc.cores = parallel::detectCores() - 1L)

data <- readRDS("data/random_data.rds")

dcp_to_ros_data <- na.omit(data[, c("conc", "ros")])
dcp_to_ros_data$conc <- ifelse(dcp_to_ros_data$conc == 0, .Machine$double.eps, dcp_to_ros_data$conc)

## Model

dcp_to_ros_formula <- brmsformula(formula = ros ~ lower + (upper - lower) / (1 + exp(-slope * log(conc / ec50)))^shape,
                                     flist = list(slope ~ 1, lower ~ 1, upper ~ 1, ec50 ~ 1, shape ~ 1),
                                     family = poisson(link = "identity"),
                                     nl = TRUE)

get_prior(formula = dcp_to_ros_formula,
          data = dcp_to_ros_data)

plot(ros ~ conc, data = dcp_to_ros_data)

dcp_to_ros_priors <- c(prior(prior = normal(   2,    1),   class = "b", nlpar = "ec50",  lb = 0),
                       prior(prior = normal(3000,  500),   class = "b", nlpar = "lower", lb = 0),
                       prior(prior = normal(   0.5,  0.5), class = "b", nlpar = "shape", lb = 0),
                       prior(prior = normal(  20,    5),   class = "b", nlpar = "slope", lb = 0),
                       prior(prior = normal(8000, 1000),   class = "b", nlpar = "upper", lb = 0))

dcp_to_ros <- brm(formula = dcp_to_ros_formula,
                  data = dcp_to_ros_data,
                  prior = dcp_to_ros_priors,
                  chains = 5,
                  iter = 20000,
                  warmup = 18000,
                  control = list(adapt_delta = 0.99999999, max_treedepth = 20),
                  file = "models/dcp_to_ros")

summary(dcp_to_ros)

plot(dcp_to_ros)

pp_check(dcp_to_ros, type = "dens_overlay")

pp_check(dcp_to_ros, type = "ecdf_overlay")

conditional_effects(dcp_to_ros)

fixef(dcp_to_ros)

## Simulations

dcp_to_ros_pred <- data.frame(conc = seq(2/3 * min(data$conc, na.rm = TRUE), 1.5 * max(data$conc, na.rm = TRUE), length.out = 10000))

dcp_to_ros_pred$ros_pred <- as.vector(posterior_predict(dcp_to_ros, newdata = dcp_to_ros_pred, nsamples = 1))

head(dcp_to_ros_pred)

plot(ros_pred ~ conc, data = dcp_to_ros_pred, type = "l")
points(ros ~ conc, data = dcp_to_ros_data, col = "red")

write.csv(dcp_to_ros_pred, "simulations/dcp_to_ros_pred.csv")
saveRDS(dcp_to_ros_pred, "simulations/dcp_to_ros_pred.rds")

summary(drc::drm(formula = ros ~ conc, data = dcp_to_ros_data, fct = drc::LL.5()))
