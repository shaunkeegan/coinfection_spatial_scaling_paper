##  Runs the brm model, given the 'model.data' dataset, and the specified 'var.names' neighbourhood explanatory variable,
##    Does this for each neighbourhood size in turn 
##    Returns estimates and lwr/upr CI for all terms in the model (inc intercept and Random Effect)
##  2 models: model.run.intensity for intensity (gaussian) or model.run.bernoulli for presence/absence (bernoulli) models



## ------------These 2 functions run the models -----------#
##  Called by functions below, passing relevant variables for each NH size

# 1) For intensity (EPG) response variables (gaussian models)
model.run.intensity.model <- function(model.data, tothosts, Eimprev){
  
  model.data$totalhosts <- tothosts
  model.data$EimNHprev <- Eimprev
    
  mod <- brm(log(response.var) ~ pred.var + totalhosts + EimNHprev + HpolINF + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
             model.data, iter = 4000,
             chains = 5,
             control = list(adapt_delta = 0.99),
             prior = prior,
             save_pars = save_pars(all=TRUE))
  
  #Compile output, comprising model fixed effects and random effects (Estimates and CIs)
  op <- rbind(fixef(mod, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975)),
              c(summary(mod)$random$ID$Estimate, summary(mod)$random$ID$Est.Error, summary(mod)$random$ID[1,"u-95% CI"], summary(mod)$random$ID[1,"l-95% CI"]))
  rownames(op)[10] <- "RanEff"
  
  return(op)
} #Close model.run.intensity.model


# 2) For presence/absence response variables (bernoulli models)
model.run.bernoulli.model <- function(model.data, tothosts, Eimprev){
  
  model.data$totalhosts <- tothosts
  model.data$EimNHprev <- Eimprev
  
  mod <- brm(response.var ~ pred.var + totalhosts + EimNHprev + HpolINF + Sex + Age + poly(as.Date(Capture.date), 2) + (1| ID), 
             model.data, family=bernoulli(),
             iter = 4000, chains = 5,
             control = list(adapt_delta = 0.99),
             prior = prior,
             save_pars = save_pars(all=TRUE))
  
  #Compile output, comprising model fixed effects and random effects (Estimates and CIs)
  op <- rbind(fixef(mod, summary = TRUE, robust = TRUE, probs = c(0.025, 0.975)),
              c(summary(mod)$random$ID$Estimate, summary(mod)$random$ID$Est.Error, summary(mod)$random$ID[1,"u-95% CI"], summary(mod)$random$ID[1,"l-95% CI"]))
  rownames(op)[10] <- "RanEff"
  
  return(op)
} #Close model.run.bernoulli.model


#-----------------------------------------------------------------





## ------------These 2 functions control the running of the models -----------#
##  Formats data structure with relevant variables for each NH size, to pass to above functions 

# 1) INTENSITY (gaussian) MODELS
model.run.intensity.control <- function(model.data, var.names){
  
  #Set up blank columns for estimates & CIs for neighbourhood parasite variable (variable of interest)
  estimates <- rep(NA, length(var.names))
  lowerCI <- rep(NA, length(var.names))
  upperCI <- rep(NA, length(var.names))
  
  #Set up blank cols for estimates & CIs for other variables in model (not of specific interest, but for reporting)
  estimates.tothosts <- rep(NA, length(var.names))
  lowerCI.tothosts <- rep(NA, length(var.names))
  upperCI.tothosts <- rep(NA, length(var.names))
  estimates.EimNHprev <- rep(NA, length(var.names))
  lowerCI.EimNHprev <- rep(NA, length(var.names))
  upperCI.EimNHprev <- rep(NA, length(var.names))
  estimates.HpolINF <- rep(NA, length(var.names))
  lowerCI.HpolINF <- rep(NA, length(var.names))
  upperCI.HpolINF <- rep(NA, length(var.names))
  estimates.Sex <- rep(NA, length(var.names))
  lowerCI.Sex <- rep(NA, length(var.names))
  upperCI.Sex <- rep(NA, length(var.names))
  estimates.Age <- rep(NA, length(var.names))
  lowerCI.Age <- rep(NA, length(var.names))
  upperCI.Age <- rep(NA, length(var.names))
  estimates.capdate1 <- rep(NA, length(var.names))
  lowerCI.capdate1 <- rep(NA, length(var.names))
  upperCI.capdate1 <- rep(NA, length(var.names))
  estimates.capdate2 <- rep(NA, length(var.names))
  lowerCI.capdate2 <- rep(NA, length(var.names))
  upperCI.capdate2 <- rep(NA, length(var.names))
  estimates.intercept <- rep(NA, length(var.names))
  lowerCI.intercept <- rep(NA, length(var.names))
  upperCI.intercept <- rep(NA, length(var.names))
  estimates.raneff <- rep(NA, length(var.names))
  lowerCI.raneff <- rep(NA, length(var.names))
  upperCI.raneff <- rep(NA, length(var.names))
  
  
  #drop all 0-infected hosts (analyse intensity, not abundance)
  model.data <- subset(model.data, response.var > 0) 
  
  
  ## Run the model for each NH size in turn, storing the estimates
  model.data$pred.var <- model.data[,var.names[1]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_10, model.data$EimPrev_10)
  
  estimates[1] <- temp["pred.var","Estimate"]
  lowerCI[1] <- temp["pred.var","Q2.5"]
  upperCI[1] <- temp["pred.var","Q97.5"]
  estimates.tothosts[1] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[1] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[1] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[1] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[1] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[1] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[1] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[1] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[1] <- temp["HpolINF","Q97.5"]
  estimates.Sex[1] <- temp["SexM","Estimate"]
  lowerCI.Sex[1] <- temp["SexM","Q2.5"]
  upperCI.Sex[1] <- temp["SexM","Q97.5"]
  estimates.Age[1] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[1] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[1] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[1] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[1] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[1] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[1] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[1] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[1] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[1] <- temp["Intercept","Estimate"]
  lowerCI.intercept[1] <- temp["Intercept","Q2.5"]
  upperCI.intercept[1] <- temp["Intercept","Q97.5"]
  estimates.raneff[1] <- temp["RanEff","Estimate"]
  lowerCI.raneff[1] <- temp["RanEff","Q2.5"]
  upperCI.raneff[1] <- temp["RanEff","Q97.5"]
  
  
  model.data$pred.var <- model.data[,var.names[2]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_14.1, model.data$EimPrev_14.1)
  
  estimates[2] <- temp["pred.var","Estimate"]
  lowerCI[2] <- temp["pred.var","Q2.5"]
  upperCI[2] <- temp["pred.var","Q97.5"]
  estimates.tothosts[2] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[2] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[2] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[2] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[2] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[2] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[2] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[2] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[2] <- temp["HpolINF","Q97.5"]
  estimates.Sex[2] <- temp["SexM","Estimate"]
  lowerCI.Sex[2] <- temp["SexM","Q2.5"]
  upperCI.Sex[2] <- temp["SexM","Q97.5"]
  estimates.Age[2] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[2] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[2] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[2] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[2] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[2] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[2] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[2] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[2] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[2] <- temp["Intercept","Estimate"]
  lowerCI.intercept[2] <- temp["Intercept","Q2.5"]
  upperCI.intercept[2] <- temp["Intercept","Q97.5"]
  estimates.raneff[2] <- temp["RanEff","Estimate"]
  lowerCI.raneff[2] <- temp["RanEff","Q2.5"]
  upperCI.raneff[2] <- temp["RanEff","Q97.5"]
  
  
  model.data$pred.var <- model.data[,var.names[3]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_20, model.data$EimPrev_20)
  
  estimates[3] <- temp["pred.var","Estimate"]
  lowerCI[3] <- temp["pred.var","Q2.5"]
  upperCI[3] <- temp["pred.var","Q97.5"]
  estimates.tothosts[3] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[3] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[3] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[3] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[3] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[3] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[3] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[3] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[3] <- temp["HpolINF","Q97.5"]
  estimates.Sex[3] <- temp["SexM","Estimate"]
  lowerCI.Sex[3] <- temp["SexM","Q2.5"]
  upperCI.Sex[3] <- temp["SexM","Q97.5"]
  estimates.Age[3] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[3] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[3] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[3] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[3] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[3] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[3] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[3] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[3] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[3] <- temp["Intercept","Estimate"]
  lowerCI.intercept[3] <- temp["Intercept","Q2.5"]
  upperCI.intercept[3] <- temp["Intercept","Q97.5"]
  estimates.raneff[3] <- temp["RanEff","Estimate"]
  lowerCI.raneff[3] <- temp["RanEff","Q2.5"]
  upperCI.raneff[3] <- temp["RanEff","Q97.5"]
  
  
  model.data$pred.var <- model.data[,var.names[4]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_22.4, model.data$EimPrev_22.4)
  
  estimates[4] <- temp["pred.var","Estimate"]
  lowerCI[4] <- temp["pred.var","Q2.5"]
  upperCI[4] <- temp["pred.var","Q97.5"]
  estimates.tothosts[4] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[4] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[4] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[4] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[4] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[4] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[4] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[4] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[4] <- temp["HpolINF","Q97.5"]
  estimates.Sex[4] <- temp["SexM","Estimate"]
  lowerCI.Sex[4] <- temp["SexM","Q2.5"]
  upperCI.Sex[4] <- temp["SexM","Q97.5"]
  estimates.Age[4] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[4] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[4] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[4] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[4] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[4] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[4] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[4] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[4] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[4] <- temp["Intercept","Estimate"]
  lowerCI.intercept[4] <- temp["Intercept","Q2.5"]
  upperCI.intercept[4] <- temp["Intercept","Q97.5"]
  estimates.raneff[4] <- temp["RanEff","Estimate"]
  lowerCI.raneff[4] <- temp["RanEff","Q2.5"]
  upperCI.raneff[4] <- temp["RanEff","Q97.5"]  
  
  
  
  model.data$pred.var <- model.data[,var.names[5]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_28.3, model.data$EimPrev_28.3)
  
  estimates[5] <- temp["pred.var","Estimate"]
  lowerCI[5] <- temp["pred.var","Q2.5"]
  upperCI[5] <- temp["pred.var","Q97.5"]
  estimates.tothosts[5] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[5] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[5] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[5] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[5] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[5] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[5] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[5] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[5] <- temp["HpolINF","Q97.5"]
  estimates.Sex[5] <- temp["SexM","Estimate"]
  lowerCI.Sex[5] <- temp["SexM","Q2.5"]
  upperCI.Sex[5] <- temp["SexM","Q97.5"]
  estimates.Age[5] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[5] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[5] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[5] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[5] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[5] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[5] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[5] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[5] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[5] <- temp["Intercept","Estimate"]
  lowerCI.intercept[5] <- temp["Intercept","Q2.5"]
  upperCI.intercept[5] <- temp["Intercept","Q97.5"]
  estimates.raneff[5] <- temp["RanEff","Estimate"]
  lowerCI.raneff[5] <- temp["RanEff","Q2.5"]
  upperCI.raneff[5] <- temp["RanEff","Q97.5"]  
  
  
  
  model.data$pred.var <- model.data[,var.names[6]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_30, model.data$EimPrev_30)
  
  estimates[6] <- temp["pred.var","Estimate"]
  lowerCI[6] <- temp["pred.var","Q2.5"]
  upperCI[6] <- temp["pred.var","Q97.5"]
  estimates.tothosts[6] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[6] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[6] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[6] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[6] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[6] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[6] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[6] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[6] <- temp["HpolINF","Q97.5"]
  estimates.Sex[6] <- temp["SexM","Estimate"]
  lowerCI.Sex[6] <- temp["SexM","Q2.5"]
  upperCI.Sex[6] <- temp["SexM","Q97.5"]
  estimates.Age[6] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[6] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[6] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[6] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[6] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[6] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[6] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[6] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[6] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[6] <- temp["Intercept","Estimate"]
  lowerCI.intercept[6] <- temp["Intercept","Q2.5"]
  upperCI.intercept[6] <- temp["Intercept","Q97.5"]
  estimates.raneff[6] <- temp["RanEff","Estimate"]
  lowerCI.raneff[6] <- temp["RanEff","Q2.5"]
  upperCI.raneff[6] <- temp["RanEff","Q97.5"]  
  
  
  
  model.data$pred.var <- model.data[,var.names[7]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_31.6, model.data$EimPrev_31.6)
  
  estimates[7] <- temp["pred.var","Estimate"]
  lowerCI[7] <- temp["pred.var","Q2.5"]
  upperCI[7] <- temp["pred.var","Q97.5"]
  estimates.tothosts[7] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[7] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[7] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[7] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[7] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[7] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[7] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[7] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[7] <- temp["HpolINF","Q97.5"]
  estimates.Sex[7] <- temp["SexM","Estimate"]
  lowerCI.Sex[7] <- temp["SexM","Q2.5"]
  upperCI.Sex[7] <- temp["SexM","Q97.5"]
  estimates.Age[7] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[7] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[7] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[7] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[7] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[7] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[7] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[7] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[7] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[7] <- temp["Intercept","Estimate"]
  lowerCI.intercept[7] <- temp["Intercept","Q2.5"]
  upperCI.intercept[7] <- temp["Intercept","Q97.5"]
  estimates.raneff[7] <- temp["RanEff","Estimate"]
  lowerCI.raneff[7] <- temp["RanEff","Q2.5"]
  upperCI.raneff[7] <- temp["RanEff","Q97.5"] 
  
  model.data$pred.var <- model.data[,var.names[8]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_36.1, model.data$EimPrev_36.1)
  
  estimates[8] <- temp["pred.var","Estimate"]
  lowerCI[8] <- temp["pred.var","Q2.5"]
  upperCI[8] <- temp["pred.var","Q97.5"]
  estimates.tothosts[8] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[8] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[8] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[8] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[8] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[8] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[8] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[8] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[8] <- temp["HpolINF","Q97.5"]
  estimates.Sex[8] <- temp["SexM","Estimate"]
  lowerCI.Sex[8] <- temp["SexM","Q2.5"]
  upperCI.Sex[8] <- temp["SexM","Q97.5"]
  estimates.Age[8] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[8] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[8] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[8] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[8] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[8] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[8] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[8] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[8] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[8] <- temp["Intercept","Estimate"]
  lowerCI.intercept[8] <- temp["Intercept","Q2.5"]
  upperCI.intercept[8] <- temp["Intercept","Q97.5"]
  estimates.raneff[8] <- temp["RanEff","Estimate"]
  lowerCI.raneff[8] <- temp["RanEff","Q2.5"]
  upperCI.raneff[8] <- temp["RanEff","Q97.5"]  
  
  model.data$pred.var <- model.data[,var.names[9]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_40, model.data$EimPrev_40)
  
  estimates[9] <- temp["pred.var","Estimate"]
  lowerCI[9] <- temp["pred.var","Q2.5"]
  upperCI[9] <- temp["pred.var","Q97.5"]
  estimates.tothosts[9] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[9] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[9] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[9] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[9] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[9] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[9] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[9] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[9] <- temp["HpolINF","Q97.5"]
  estimates.Sex[9] <- temp["SexM","Estimate"]
  lowerCI.Sex[9] <- temp["SexM","Q2.5"]
  upperCI.Sex[9] <- temp["SexM","Q97.5"]
  estimates.Age[9] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[9] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[9] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[9] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[9] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[9] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[9] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[9] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[9] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[9] <- temp["Intercept","Estimate"]
  lowerCI.intercept[9] <- temp["Intercept","Q2.5"]
  upperCI.intercept[9] <- temp["Intercept","Q97.5"]
  estimates.raneff[9] <- temp["RanEff","Estimate"]
  lowerCI.raneff[9] <- temp["RanEff","Q2.5"]
  upperCI.raneff[9] <- temp["RanEff","Q97.5"]  
  
  model.data$pred.var <- model.data[,var.names[10]]
  
  temp <- model.run.intensity.model(model.data, model.data$NemNoTotal_50, model.data$EimPrev_50)
  
  estimates[10] <- temp["pred.var","Estimate"]
  lowerCI[10] <- temp["pred.var","Q2.5"]
  upperCI[10] <- temp["pred.var","Q97.5"]
  estimates.tothosts[10] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[10] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[10] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[10] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[10] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[10] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[10] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[10] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[10] <- temp["HpolINF","Q97.5"]
  estimates.Sex[10] <- temp["SexM","Estimate"]
  lowerCI.Sex[10] <- temp["SexM","Q2.5"]
  upperCI.Sex[10] <- temp["SexM","Q97.5"]
  estimates.Age[10] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[10] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[10] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[10] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[10] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[10] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[10] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[10] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[10] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[10] <- temp["Intercept","Estimate"]
  lowerCI.intercept[10] <- temp["Intercept","Q2.5"]
  upperCI.intercept[10] <- temp["Intercept","Q97.5"]
  estimates.raneff[10] <- temp["RanEff","Estimate"]
  lowerCI.raneff[10] <- temp["RanEff","Q2.5"]
  upperCI.raneff[10] <- temp["RanEff","Q97.5"]
  
  
  output <- cbind(dist, estimates, lowerCI, upperCI,
                  estimates.tothosts, lowerCI.tothosts, upperCI.tothosts,
                  estimates.EimNHprev, lowerCI.EimNHprev, upperCI.EimNHprev,
                  estimates.HpolINF,lowerCI.HpolINF, upperCI.HpolINF,
                  estimates.Sex, lowerCI.Sex, upperCI.Sex,
                  estimates.Age, lowerCI.Age, upperCI.Age,
                  estimates.capdate1, lowerCI.capdate1, upperCI.capdate1,
                  estimates.capdate2, lowerCI.capdate2, upperCI.capdate2,
                  estimates.intercept, lowerCI.intercept, upperCI.intercept,
                  estimates.raneff, lowerCI.raneff, upperCI.raneff
  )
  
  return(output)
}






#---------------------------------------------------------------------#



# 2) PRESENCE/ABSENCE (bernoulli) MODEL
model.run.bernoulli.control <- function(model.data, var.names){
  
  #Set up blank columns for estimates & CIs for neighbourhood parasite variable (variable of interest)
  estimates <- rep(NA, length(var.names))
  lowerCI <- rep(NA, length(var.names))
  upperCI <- rep(NA, length(var.names))
  
  #Set up blank cols for estimates & CIs for other variables in model (not of specific interest, but for reporting)
  estimates.tothosts <- rep(NA, length(var.names))
  lowerCI.tothosts <- rep(NA, length(var.names))
  upperCI.tothosts <- rep(NA, length(var.names))
  estimates.EimNHprev <- rep(NA, length(var.names))
  lowerCI.EimNHprev <- rep(NA, length(var.names))
  upperCI.EimNHprev <- rep(NA, length(var.names))
  estimates.HpolINF <- rep(NA, length(var.names))
  lowerCI.HpolINF <- rep(NA, length(var.names))
  upperCI.HpolINF <- rep(NA, length(var.names))
  estimates.Sex <- rep(NA, length(var.names))
  lowerCI.Sex <- rep(NA, length(var.names))
  upperCI.Sex <- rep(NA, length(var.names))
  estimates.Age <- rep(NA, length(var.names))
  lowerCI.Age <- rep(NA, length(var.names))
  upperCI.Age <- rep(NA, length(var.names))
  estimates.capdate1 <- rep(NA, length(var.names))
  lowerCI.capdate1 <- rep(NA, length(var.names))
  upperCI.capdate1 <- rep(NA, length(var.names))
  estimates.capdate2 <- rep(NA, length(var.names))
  lowerCI.capdate2 <- rep(NA, length(var.names))
  upperCI.capdate2 <- rep(NA, length(var.names))
  estimates.intercept <- rep(NA, length(var.names))
  lowerCI.intercept <- rep(NA, length(var.names))
  upperCI.intercept <- rep(NA, length(var.names))
  estimates.raneff <- rep(NA, length(var.names))
  lowerCI.raneff <- rep(NA, length(var.names))
  upperCI.raneff <- rep(NA, length(var.names))
  
  
  #Run the model for each NH size in turn, storing the estimates
  model.data$pred.var <- model.data[,var.names[1]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_10, model.data$EimPrev_10)

  estimates[1] <- temp["pred.var","Estimate"]
  lowerCI[1] <- temp["pred.var","Q2.5"]
  upperCI[1] <- temp["pred.var","Q97.5"]
  estimates.tothosts[1] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[1] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[1] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[1] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[1] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[1] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[1] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[1] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[1] <- temp["HpolINF","Q97.5"]
  estimates.Sex[1] <- temp["SexM","Estimate"]
  lowerCI.Sex[1] <- temp["SexM","Q2.5"]
  upperCI.Sex[1] <- temp["SexM","Q97.5"]
  estimates.Age[1] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[1] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[1] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[1] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[1] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[1] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[1] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[1] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[1] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[1] <- temp["Intercept","Estimate"]
  lowerCI.intercept[1] <- temp["Intercept","Q2.5"]
  upperCI.intercept[1] <- temp["Intercept","Q97.5"]
  estimates.raneff[1] <- temp["RanEff","Estimate"]
  lowerCI.raneff[1] <- temp["RanEff","Q2.5"]
  upperCI.raneff[1] <- temp["RanEff","Q97.5"]

    
  model.data$pred.var <- model.data[,var.names[2]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_14.1, model.data$EimPrev_14.1)
  
  estimates[2] <- temp["pred.var","Estimate"]
  lowerCI[2] <- temp["pred.var","Q2.5"]
  upperCI[2] <- temp["pred.var","Q97.5"]
  estimates.tothosts[2] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[2] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[2] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[2] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[2] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[2] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[2] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[2] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[2] <- temp["HpolINF","Q97.5"]
  estimates.Sex[2] <- temp["SexM","Estimate"]
  lowerCI.Sex[2] <- temp["SexM","Q2.5"]
  upperCI.Sex[2] <- temp["SexM","Q97.5"]
  estimates.Age[2] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[2] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[2] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[2] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[2] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[2] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[2] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[2] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[2] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[2] <- temp["Intercept","Estimate"]
  lowerCI.intercept[2] <- temp["Intercept","Q2.5"]
  upperCI.intercept[2] <- temp["Intercept","Q97.5"]
  estimates.raneff[2] <- temp["RanEff","Estimate"]
  lowerCI.raneff[2] <- temp["RanEff","Q2.5"]
  upperCI.raneff[2] <- temp["RanEff","Q97.5"]
  
  
  model.data$pred.var <- model.data[,var.names[3]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_20, model.data$EimPrev_20)
  
  estimates[3] <- temp["pred.var","Estimate"]
  lowerCI[3] <- temp["pred.var","Q2.5"]
  upperCI[3] <- temp["pred.var","Q97.5"]
  estimates.tothosts[3] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[3] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[3] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[3] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[3] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[3] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[3] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[3] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[3] <- temp["HpolINF","Q97.5"]
  estimates.Sex[3] <- temp["SexM","Estimate"]
  lowerCI.Sex[3] <- temp["SexM","Q2.5"]
  upperCI.Sex[3] <- temp["SexM","Q97.5"]
  estimates.Age[3] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[3] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[3] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[3] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[3] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[3] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[3] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[3] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[3] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[3] <- temp["Intercept","Estimate"]
  lowerCI.intercept[3] <- temp["Intercept","Q2.5"]
  upperCI.intercept[3] <- temp["Intercept","Q97.5"]
  estimates.raneff[3] <- temp["RanEff","Estimate"]
  lowerCI.raneff[3] <- temp["RanEff","Q2.5"]
  upperCI.raneff[3] <- temp["RanEff","Q97.5"]
  
  
  model.data$pred.var <- model.data[,var.names[4]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_22.4, model.data$EimPrev_22.4)
  
  estimates[4] <- temp["pred.var","Estimate"]
  lowerCI[4] <- temp["pred.var","Q2.5"]
  upperCI[4] <- temp["pred.var","Q97.5"]
  estimates.tothosts[4] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[4] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[4] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[4] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[4] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[4] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[4] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[4] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[4] <- temp["HpolINF","Q97.5"]
  estimates.Sex[4] <- temp["SexM","Estimate"]
  lowerCI.Sex[4] <- temp["SexM","Q2.5"]
  upperCI.Sex[4] <- temp["SexM","Q97.5"]
  estimates.Age[4] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[4] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[4] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[4] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[4] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[4] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[4] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[4] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[4] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[4] <- temp["Intercept","Estimate"]
  lowerCI.intercept[4] <- temp["Intercept","Q2.5"]
  upperCI.intercept[4] <- temp["Intercept","Q97.5"]
  estimates.raneff[4] <- temp["RanEff","Estimate"]
  lowerCI.raneff[4] <- temp["RanEff","Q2.5"]
  upperCI.raneff[4] <- temp["RanEff","Q97.5"]  

  
    
  model.data$pred.var <- model.data[,var.names[5]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_28.3, model.data$EimPrev_28.3)
  
  estimates[5] <- temp["pred.var","Estimate"]
  lowerCI[5] <- temp["pred.var","Q2.5"]
  upperCI[5] <- temp["pred.var","Q97.5"]
  estimates.tothosts[5] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[5] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[5] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[5] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[5] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[5] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[5] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[5] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[5] <- temp["HpolINF","Q97.5"]
  estimates.Sex[5] <- temp["SexM","Estimate"]
  lowerCI.Sex[5] <- temp["SexM","Q2.5"]
  upperCI.Sex[5] <- temp["SexM","Q97.5"]
  estimates.Age[5] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[5] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[5] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[5] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[5] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[5] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[5] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[5] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[5] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[5] <- temp["Intercept","Estimate"]
  lowerCI.intercept[5] <- temp["Intercept","Q2.5"]
  upperCI.intercept[5] <- temp["Intercept","Q97.5"]
  estimates.raneff[5] <- temp["RanEff","Estimate"]
  lowerCI.raneff[5] <- temp["RanEff","Q2.5"]
  upperCI.raneff[5] <- temp["RanEff","Q97.5"]  

  
    
  model.data$pred.var <- model.data[,var.names[6]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_30, model.data$EimPrev_30)
  
  estimates[6] <- temp["pred.var","Estimate"]
  lowerCI[6] <- temp["pred.var","Q2.5"]
  upperCI[6] <- temp["pred.var","Q97.5"]
  estimates.tothosts[6] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[6] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[6] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[6] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[6] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[6] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[6] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[6] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[6] <- temp["HpolINF","Q97.5"]
  estimates.Sex[6] <- temp["SexM","Estimate"]
  lowerCI.Sex[6] <- temp["SexM","Q2.5"]
  upperCI.Sex[6] <- temp["SexM","Q97.5"]
  estimates.Age[6] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[6] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[6] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[6] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[6] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[6] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[6] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[6] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[6] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[6] <- temp["Intercept","Estimate"]
  lowerCI.intercept[6] <- temp["Intercept","Q2.5"]
  upperCI.intercept[6] <- temp["Intercept","Q97.5"]
  estimates.raneff[6] <- temp["RanEff","Estimate"]
  lowerCI.raneff[6] <- temp["RanEff","Q2.5"]
  upperCI.raneff[6] <- temp["RanEff","Q97.5"]  

  
    
  model.data$pred.var <- model.data[,var.names[7]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_31.6, model.data$EimPrev_31.6)
  
  estimates[7] <- temp["pred.var","Estimate"]
  lowerCI[7] <- temp["pred.var","Q2.5"]
  upperCI[7] <- temp["pred.var","Q97.5"]
  estimates.tothosts[7] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[7] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[7] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[7] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[7] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[7] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[7] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[7] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[7] <- temp["HpolINF","Q97.5"]
  estimates.Sex[7] <- temp["SexM","Estimate"]
  lowerCI.Sex[7] <- temp["SexM","Q2.5"]
  upperCI.Sex[7] <- temp["SexM","Q97.5"]
  estimates.Age[7] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[7] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[7] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[7] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[7] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[7] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[7] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[7] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[7] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[7] <- temp["Intercept","Estimate"]
  lowerCI.intercept[7] <- temp["Intercept","Q2.5"]
  upperCI.intercept[7] <- temp["Intercept","Q97.5"]
  estimates.raneff[7] <- temp["RanEff","Estimate"]
  lowerCI.raneff[7] <- temp["RanEff","Q2.5"]
  upperCI.raneff[7] <- temp["RanEff","Q97.5"] 
  
  model.data$pred.var <- model.data[,var.names[8]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_36.1, model.data$EimPrev_36.1)
  
  estimates[8] <- temp["pred.var","Estimate"]
  lowerCI[8] <- temp["pred.var","Q2.5"]
  upperCI[8] <- temp["pred.var","Q97.5"]
  estimates.tothosts[8] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[8] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[8] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[8] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[8] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[8] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[8] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[8] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[8] <- temp["HpolINF","Q97.5"]
  estimates.Sex[8] <- temp["SexM","Estimate"]
  lowerCI.Sex[8] <- temp["SexM","Q2.5"]
  upperCI.Sex[8] <- temp["SexM","Q97.5"]
  estimates.Age[8] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[8] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[8] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[8] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[8] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[8] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[8] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[8] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[8] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[8] <- temp["Intercept","Estimate"]
  lowerCI.intercept[8] <- temp["Intercept","Q2.5"]
  upperCI.intercept[8] <- temp["Intercept","Q97.5"]
  estimates.raneff[8] <- temp["RanEff","Estimate"]
  lowerCI.raneff[8] <- temp["RanEff","Q2.5"]
  upperCI.raneff[8] <- temp["RanEff","Q97.5"]  
  
  model.data$pred.var <- model.data[,var.names[9]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_40, model.data$EimPrev_40)
  
  estimates[9] <- temp["pred.var","Estimate"]
  lowerCI[9] <- temp["pred.var","Q2.5"]
  upperCI[9] <- temp["pred.var","Q97.5"]
  estimates.tothosts[9] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[9] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[9] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[9] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[9] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[9] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[9] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[9] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[9] <- temp["HpolINF","Q97.5"]
  estimates.Sex[9] <- temp["SexM","Estimate"]
  lowerCI.Sex[9] <- temp["SexM","Q2.5"]
  upperCI.Sex[9] <- temp["SexM","Q97.5"]
  estimates.Age[9] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[9] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[9] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[9] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[9] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[9] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[9] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[9] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[9] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[9] <- temp["Intercept","Estimate"]
  lowerCI.intercept[9] <- temp["Intercept","Q2.5"]
  upperCI.intercept[9] <- temp["Intercept","Q97.5"]
  estimates.raneff[9] <- temp["RanEff","Estimate"]
  lowerCI.raneff[9] <- temp["RanEff","Q2.5"]
  upperCI.raneff[9] <- temp["RanEff","Q97.5"]  
  
  model.data$pred.var <- model.data[,var.names[10]]
  
  temp <- model.run.bernoulli.model(model.data, model.data$NemNoTotal_50, model.data$EimPrev_50)
  
  estimates[10] <- temp["pred.var","Estimate"]
  lowerCI[10] <- temp["pred.var","Q2.5"]
  upperCI[10] <- temp["pred.var","Q97.5"]
  estimates.tothosts[10] <- temp["totalhosts","Estimate"]
  lowerCI.tothosts[10] <- temp["totalhosts","Q2.5"]
  upperCI.tothosts[10] <- temp["totalhosts","Q97.5"]
  estimates.EimNHprev[10] <- temp["EimNHprev","Estimate"]
  lowerCI.EimNHprev[10] <- temp["EimNHprev","Q2.5"]
  upperCI.EimNHprev[10] <- temp["EimNHprev","Q97.5"]
  estimates.HpolINF[10] <- temp["HpolINF","Estimate"]
  lowerCI.HpolINF[10] <- temp["HpolINF","Q2.5"]
  upperCI.HpolINF[10] <- temp["HpolINF","Q97.5"]
  estimates.Sex[10] <- temp["SexM","Estimate"]
  lowerCI.Sex[10] <- temp["SexM","Q2.5"]
  upperCI.Sex[10] <- temp["SexM","Q97.5"]
  estimates.Age[10] <- temp["AgeNotAdult","Estimate"]
  lowerCI.Age[10] <- temp["AgeNotAdult","Q2.5"]
  upperCI.Age[10] <- temp["AgeNotAdult","Q97.5"]
  estimates.capdate1[10] <- temp["polyas.DateCapture.date21","Estimate"]
  lowerCI.capdate1[10] <- temp["polyas.DateCapture.date21","Q2.5"]
  upperCI.capdate1[10] <- temp["polyas.DateCapture.date21","Q97.5"]
  estimates.capdate2[10] <- temp["polyas.DateCapture.date22","Estimate"]
  lowerCI.capdate2[10] <- temp["polyas.DateCapture.date22","Q2.5"]
  upperCI.capdate2[10] <- temp["polyas.DateCapture.date22","Q97.5"]
  estimates.intercept[10] <- temp["Intercept","Estimate"]
  lowerCI.intercept[10] <- temp["Intercept","Q2.5"]
  upperCI.intercept[10] <- temp["Intercept","Q97.5"]
  estimates.raneff[10] <- temp["RanEff","Estimate"]
  lowerCI.raneff[10] <- temp["RanEff","Q2.5"]
  upperCI.raneff[10] <- temp["RanEff","Q97.5"]
  
  
  output <- cbind(dist, estimates, lowerCI, upperCI,
                  estimates.tothosts, lowerCI.tothosts, upperCI.tothosts,
                  estimates.EimNHprev, lowerCI.EimNHprev, upperCI.EimNHprev,
                  estimates.HpolINF,lowerCI.HpolINF, upperCI.HpolINF,
                  estimates.Sex, lowerCI.Sex, upperCI.Sex,
                  estimates.Age, lowerCI.Age, upperCI.Age,
                  estimates.capdate1, lowerCI.capdate1, upperCI.capdate1,
                  estimates.capdate2, lowerCI.capdate2, upperCI.capdate2,
                  estimates.intercept, lowerCI.intercept, upperCI.intercept,
                  estimates.raneff, lowerCI.raneff, upperCI.raneff
  )
  
  return(output)
}
