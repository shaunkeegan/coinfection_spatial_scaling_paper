## Runs the brm model, given the 'model.data' dataset, and the specified 'var.names' explanatory variable,
##   Does this for each neighbourhood size in turn 
##   Returns estimate and lwr/upr CI of estimate of neighbourhood prevalence of effecting parasite

model.run <- function(model.data, var.names){
  
  estimates <- rep(NA, length(var.names))
  lowerCI <- rep(NA, length(var.names))
  upperCI <- rep(NA, length(var.names))
  
  model.10 <- brm(log(response.var) ~ Prev_10 + NoTotal_10 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                  model.data, iter = 4000,
                  chains = 10,
                  control = list(adapt_delta = 0.99),
                  prior = prior,
                  save_all_pars = TRUE)

  
  temp <- fixef(model.10, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[1] <- temp[var.names[1],"Estimate"]
  lowerCI[1] <- temp[var.names[1],"Q2.5"]
  upperCI[1] <- temp[var.names[1],"Q97.5"]
  
  
  model.14.1 <- brm(log(response.var) ~ Prev_14.1 + NoTotal_14.1 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                    model.data, iter = 4000,
                    chains = 10,
                    control = list(adapt_delta = 0.99),
                    prior = prior,
                    save_all_pars = TRUE)
  
  temp <- fixef(model.14.1, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[2] <- temp[var.names[2],"Estimate"]
  lowerCI[2] <- temp[var.names[2],"Q2.5"]
  upperCI[2] <- temp[var.names[2],"Q97.5"]
  
  model.20 <- brm(log(response.var) ~ Prev_20 + NoTotal_20 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                  model.data, iter = 4000,
                  chains = 10,
                  control = list(adapt_delta = 0.99),
                  prior = prior,
                  save_all_pars = TRUE)
  
  temp <- fixef(model.20, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[3] <- temp[var.names[3],"Estimate"]
  lowerCI[3] <- temp[var.names[3],"Q2.5"]
  upperCI[3] <- temp[var.names[3],"Q97.5"]
  
  model.22.4 <- brm(log(response.var) ~ Prev_22.4 + NoTotal_22.4 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                    model.data, iter = 4000,
                    chains = 10,
                    control = list(adapt_delta = 0.99),
                    prior = prior,
                    save_all_pars = TRUE)

  temp <- fixef(model.22.4, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[4] <- temp[var.names[4],"Estimate"]
  lowerCI[4] <- temp[var.names[4],"Q2.5"]
  upperCI[4] <- temp[var.names[4],"Q97.5"]
  
  model.28.3 <- brm(log(response.var) ~ Prev_28.3 + NoTotal_28.3 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                    model.data, iter = 4000,
                    chains = 10,
                    control = list(adapt_delta = 0.99),
                    prior = prior,
                    save_all_pars = TRUE)

  temp <- fixef(model.28.3, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[5] <- temp[var.names[5],"Estimate"]
  lowerCI[5] <- temp[var.names[5],"Q2.5"]
  upperCI[5] <- temp[var.names[5],"Q97.5"]
  
  
  model.30 <- brm(log(response.var) ~ Prev_30 + NoTotal_30 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                  model.data, iter = 4000,
                  chains = 10,
                  control = list(adapt_delta = 0.99),
                  prior = prior,
                  save_all_pars = TRUE)

  temp <- fixef(model.30, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[6] <- temp[var.names[6],"Estimate"]
  lowerCI[6] <- temp[var.names[6],"Q2.5"]
  upperCI[6] <- temp[var.names[6],"Q97.5"]
  
  model.31.6 <- brm(log(response.var) ~ Prev_31.6 + NoTotal_31.6 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                    model.data, iter = 4000,
                    chains = 10,
                    control = list(adapt_delta = 0.99),
                    prior = prior,
                    save_all_pars = TRUE)
  
  temp <- fixef(model.31.6, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[7] <- temp[var.names[7],"Estimate"]
  lowerCI[7] <- temp[var.names[7],"Q2.5"]
  upperCI[7] <- temp[var.names[7],"Q97.5"]
  
  model.36.1 <- brm(log(response.var) ~ Prev_36.1 + NoTotal_36.1 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                    model.data, iter = 4000,
                    chains = 10,
                    control = list(adapt_delta = 0.99),
                    prior = prior,
                    save_all_pars = TRUE)

  temp <- fixef(model.36.1, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[8] <- temp[var.names[8],"Estimate"]
  lowerCI[8] <- temp[var.names[8],"Q2.5"]
  upperCI[8] <- temp[var.names[8],"Q97.5"]
  
  model.40 <- brm(log(response.var) ~ Prev_40 + NoTotal_40 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                  model.data, iter = 4000,
                  chains = 10,
                  control = list(adapt_delta = 0.99),
                  prior = prior,
                  save_all_pars = TRUE)

  temp <- fixef(model.40, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[9] <- temp[var.names[9],"Estimate"]
  lowerCI[9] <- temp[var.names[9],"Q2.5"]
  upperCI[9] <- temp[var.names[9],"Q97.5"]
  
  model.50 <- brm(log(response.var) ~ Prev_50 + NoTotal_50 + HpolINF  + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
                  model.data, iter = 4000,
                  chains = 10,
                  control = list(adapt_delta = 0.99),
                  prior = prior,
                  save_all_pars = TRUE)

  temp <- fixef(model.50, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975))
  estimates[10] <- temp[var.names[10],"Estimate"]
  lowerCI[10] <- temp[var.names[10],"Q2.5"]
  upperCI[10] <- temp[var.names[10],"Q97.5"]
  
  output <- cbind(var.names, estimates, lowerCI, upperCI)
  
  return(output)
}