# CARar forecasts for NHS trust network time series, forecasts are computed in a global-alpha GNAR form 
# where parameters are estimated using MCMC

library(GNAR)
library(CARBayesST)
library(igraph)


CAR_forecast_squared_error <- function(test_data, pred_data, S1) {
  set.seed(2023)
  pred_vector_form = vector_form_mvts(pred_data)
  CARar_model <- ST.CARar(formula = pred_vector_form ~ 1, family = "gaussian", data = data.frame(pred_vector_form),
                          W = as.matrix(S1), AR = 1, burnin = 2000, n.sample = 22000, thin = 10)
  mu_hat = CARar_model$summary.results[1, 1]
  alpha_hat = CARar_model$summary.results[5, 1]
  tt = nrow(pred_data)
  one_step_ahead = mu_hat + alpha_hat * (pred_data[tt, ] - mu_hat)
  one_step_ahead_nc = alpha_hat * pred_data[tt, ]
  out0 = round(sum((one_step_ahead - test_data)^2), digits = decimal_digits)
  out1 = round(sum((one_step_ahead_nc - test_data)^2), digits = decimal_digits)
  print(paste0("global-alpha: ", as.character(alpha_hat)))
  print(paste0("done forecast for ", as.character(tt + 1)))
  return (c(out0, out1))
}

vector_form_mvts <- function(data_mat) {
  time_steps = nrow(data_mat)
  ts_dim = ncol(data_mat)
  vector_form = rep(0, time_steps * ts_dim)
  for (i in 1:time_steps) {
    top = 1 + ts_dim * (i - 1)
    bottom = ts_dim + ts_dim * (i - 1)
    vector_form[top:bottom] = data_mat[i, ]
  }
  return (vector_form)
}

covid_data = logMVbedMVC.vts
k = 0
t_steps = 10
initial_start = 452 - t_steps
decimal_digits = 3
S1 = as.matrix(GNARtoigraph(NHSTrustMVCAug120.net))

CARar_squared_errors <- t(vapply(c(1:10), function(x) {CAR_forecast_squared_error(covid_data[(initial_start + x), ], 
                                                                                covid_data[1:(initial_start + 
                                                                                             (x - 1)), ], S1)}, c(0.0, 0.0)))
colnames(CARar_squared_errors) <- c("centred_spe", "non_centred_mse")

save(CARar_squared_errors, file = "/Users/danielsalnikov/Documents/PhD/my_stuff/papers/corbit_paper/R_scripts/R_scripts_revised/CARar_GNAR_comparison.RData")

CARar_execution_times = c(1313.7, 461.5, 462.6, 1161.5, 468.1, 473.4, 506.4, 1111, 754.6, 1109.1)
