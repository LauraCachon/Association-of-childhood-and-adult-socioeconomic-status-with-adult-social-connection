### THIS SCRIPT THE CODE FOR THE ANALYSES OF THE PAPER: "Association of childhood and adult socioeconomic status with adult social relationships" ###
#CONTENTS: 1) Data wrangling, imputation, 2) Regression models, 3) Causal Mediation Analysis, 4) Sensitivity analyses, 5) CMA with interaction ExpxMed, 6) Supplementary Material: empirical check of the positivity assumption

### 1. Preparing the dataset: Multiple Imputation and Re-Scaling the variables ###
#1.1. Model variables and auxiliary variables as factors
vars_to_convert <- c("edu80", "edu107", "edu111", "occ80", "occ107", "occ111", 
                     "income80", "income07", "income11", "toimi80", "toimi07", 
                     "toimi11", "sex")

new_compbas[vars_to_convert] <- lapply(new_compbas[vars_to_convert], as.factor)

#1.2. Multiple imputation by chained equations
library(mice)
init = mice(new_compbas, maxit=0) 
meth = init$method
predM = init$predictorMatrix

#id has no predictive value
predM[, c("id")]=0

meth[c("age", "loneliness", "social_network", "social_relations")]="pmm" #numerical
meth[c("income80", "edu107", "income07", "edu111", "income11")]="polr" #categorical ordered
#meth[c("sex", "edu80", "occ80", "toimi80")]="logreg" #binary #no need as I have selected complete cases for baseline variables
meth[c("occ107", "occ111","toimi07", "toimi11")]= "polyreg" #categorical unordered

set.seed(103)
imputed.bas = mice(new_compbas, method=meth, predictorMatrix=predM, m=5)
imputed_datasets <- complete(imputed.bas, action = "long", include = TRUE)

#1.3. Select those with complete data for the outcomes (not imputed outcomes)
new_compout<- new_compbas %>% filter(complete.cases(loneliness, social_network, social_relations)) #social_relations=MSPSS scores
compout<-new_compout %>% select(id, loneliness, social_network, social_relations)

imp_data<-left_join(compout, imputed_datasets, by="id")
imp_data<-imp_data %>% select(id, .imp, .id, age, sex, edu80, income80, edu107, income07, loneliness=loneliness.x, social_network=social_network.x, social_relations=social_relations.x)
imp_data<-imputed.bas

#1.4. Arranging the variables #
#education (primary and low secondary education as risk)
imp_data$edu107 <-recode(imp_data$edu107, `0` = 1, `1` = 1, `2` = 0, `3` = 0) #also upper secondary as at risk
summary(imp_data$edu107)

#income (lower 25% as risk group)
imp_data$income07 <-recode(imp_data$income07, `1` = 1, `2` = 1, `3` = 0, `4` = 0, `5` = 0, `6` = 0, `7` = 0, `8` = 0) #Q1=3

## original TILS scores imp_data ##
imp_data<-imp_data %>% mutate(lon_TILS=loneliness*3)
imp_data$lon_TILS<-recode(imp_data$lon_TILS, `3`= 3, `4`= 4, `4.5`=5, `5`= 5, `6`= 6, `7`=7, `8`=8, `9`=9)

#inverse mspss score (higher meaning lower ssup)
imp_data<-imp_data %>% mutate(sr=social_relations*12)
imp_data<-imp_data %>% mutate(inverse=sr*-1)
imp_data<-imp_data %>% mutate(mspss.inv=inverse + 72) #from 12 to 60, with higher scores meaning less social support

#scale 
imp_data <- imp_data %>% mutate(social_network_sc = scale(social_network))
imp_data <- imp_data %>% mutate(mspss.inv_sc = scale(mspss.inv))

#Separate the 5 imputed datasets (with original outcome values) each of them N=1,685
imputed_dataset_1 <- subset(imp_data, .imp == 1)
summary(imputed_dataset_1)

imputed_dataset_2 <- subset(imp_data, .imp == 2)
summary(imputed_dataset_2)

imputed_dataset_3 <- subset(imp_data, .imp == 3)
summary(imputed_dataset_3)

imputed_dataset_4 <- subset(imp_data, .imp == 4)
summary(imputed_dataset_4)

imputed_dataset_5 <- subset(imp_data, .imp == 5)
summary(imputed_dataset_5)

### 2. REGRESSION MODELS ###
#packages
library(MASS)

imp_reg<-imp_data %>% filter(.imp !=0)

#function for pooling, CIs and OR (ordinal and logistic)
alltherest<- function(x) {pooled_model <- pool(x)
summary_data <- summary(pooled_model)

confidence_level <- 0.95
summary_data$Lower_CI <- summary_data$estimate - qnorm(1 - (1 - confidence_level) / 2) * summary_data$std.error
summary_data$Upper_CI <- summary_data$estimate + qnorm(1 - (1 - confidence_level) / 2) * summary_data$std.error

# Calculate odds ratios for each parameter
summary_data$Odds_Ratio <- exp(summary_data$estimate)
confidence_level <- 0.95
summary_data$OR_Lower_CI <- exp(summary_data$Lower_CI)
summary_data$OR_Upper_CI <- exp(summary_data$Upper_CI)

# Print the updated data frame with odds ratios and confidence intervals
print(summary_data[, c("term", "estimate", "std.error", "Odds_Ratio", "Lower_CI", "Upper_CI", "OR_Lower_CI", "OR_Upper_CI", "p.value")])
}

#function for pooling and CIs (continuous variables)
CI<- function(x) {pooled_model <- pool(x)
summary_data <- summary(pooled_model)

confidence_level <- 0.95
summary_data$Lower_CI <- summary_data$estimate - qnorm(1 - (1 - confidence_level) / 2) * summary_data$std.error
summary_data$Upper_CI <- summary_data$estimate + qnorm(1 - (1 - confidence_level) / 2) * summary_data$std.error

# Print the updated data frame with odds ratios and confidence intervals
print(summary_data[, c("term", "estimate", "std.error", "Lower_CI", "Upper_CI", "p.value")])
}

#Example 1-Regression models: loneliness - education
#lon_edu exposure-outcome
# Assuming num_imputations is the correct number of imputations
num_imputations <- 5  # Change this to the actual number of imputations

fit_models <- lapply(1:num_imputations, function(i) {
  polr(lon_TILS ~ edu80 + age, Hess = TRUE, data = imp_reg %>% filter(.imp == i))
})

lon_edu_done<-alltherest(fit_models)

#edu-edu07 exposure-mediator
edu_edu07 <- lapply(1:num_imputations, function(i) {
  glm(formula = edu107 ~ edu80 + age, family = binomial(), data = imp_reg %>% filter(.imp == i))
})

edu_edu07_done<-alltherest(edu_edu07)

#lon_edu mediator-outcome 
lon_edu07 <- lapply(1:num_imputations, function(i) {
  polr(lon_TILS ~ edu107 + age + sex + income80, Hess = TRUE, data = imp_reg %>% filter(.imp == i))
})

lon_edu07_done<-alltherest(lon_edu07)

#Example 2-Regression models: perceived social support and income 
##exposure-outcome
sp_inc <- lapply(1:num_imputations, function(i) {
  lm(mspss.inv_sc ~ income80 + age +edu80, data = imp_reg %>% filter(.imp == i))
})

sp_inc_done<-CI(sp_inc)

#exposure-mediator 
inc_inc07 <- lapply(1:num_imputations, function(i) {
  glm(formula = income07 ~ income80 + age + edu80, family = binomial(), data = imp_reg %>% filter(.imp == i))
})

inc_inc07_done<-alltherest(inc_inc07)

#mediator-outcome
sp_inc07 <- lapply(1:num_imputations, function(i) {
  lm(mspss.inv_sc ~ income07 + age + sex + edu80 + edu107, data = imp_reg %>% filter(.imp == i))
})

sp_inc07<-CI(sp_inc07)


### 3. CAUSAL MEDIATION MODELS ###
library(CMAverse)
#FUNCTIONS
#Runs the CMA models in the 5 imputed datasets
#Ordinal outcome
edu_cmest_0 <- function(x, num_imputations = 5) {
  models <- list()
  
  for (i in 1:num_imputations) {
    # Create the dataset variable name (imputed_dataset_i)
    dataset_var <- paste0("imputed_dataset_", i)
    
    # Perform cmest for each imputation
    model <- cmest(
      data = get(dataset_var), 
      model = "msm", outcome = x, 
      exposure = "edu80", mediator = "edu107", 
      basec = c("age", "sex"), postc = "income80", 
      EMint = FALSE, ereg = "logistic", 
      yreg = "ordinal", mreg = list("logistic"),
      wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
      astar = 0, a = 1, mval = list(1), 
      estimation = "imputation", inference = "bootstrap", nboot=1000
    )
    
    if (!is.null(model)) {
      # Get the summary of the model
      model_summary <- summary(model)
      
      # Add the model summary to the list
      models[[paste0("m", i)]] <- model_summary
    }
  }
  
  return(models)
}

income_cmest_0 <- function(x, num_imputations = 5) {
  models <- list()
  
  for (i in 1:num_imputations) {
    # Create the dataset variable name (imputed_dataset_i)
    dataset_var <- paste0("imputed_dataset_", i)
    
    # Perform cmest for each imputation
    model <- cmest(
      data = get(dataset_var), 
      model = "msm", outcome = x, 
      exposure = "income80", mediator = "income07", basec = c("age", "sex", "edu80"), 
      postc = "edu107", EMint = F,
      ereg = "logistic", yreg = "ordinal", mreg = list("logistic"),
      wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
      astar = 0, a = 1, mval = list(1), 
      estimation = "imputation", inference = "bootstrap", nboot=1000)
    
    if (!is.null(model)) {
      # Get the summary of the model
      model_summary <- summary(model)
      
      # Add the model summary to the list
      models[[paste0("m", i)]] <- model_summary
    }
  }
  
  return(models)
}

## continuous outcomes #
edu_cmest_1 <- function(x, num_imputations = 5) {
  models <- list()
  
  for (i in 1:num_imputations) {
    # Create the dataset variable name (imputed_dataset_i)
    dataset_var <- paste0("imputed_dataset_", i)
    
    # Perform cmest for each imputation
    model <- cmest(
      data = get(dataset_var), 
      model = "msm", outcome = x, 
      exposure = "edu80", mediator = "edu107", basec = c("age", "sex"), 
      postc = "income80", EMint = F,
      ereg = "logistic", yreg = "linear", mreg = list("logistic"),
      wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
      astar = 0, a = 1, mval = list(1), 
      estimation = "imputation", inference = "bootstrap", nboot=1000)
    
    if (!is.null(model)) {
      # Get the summary of the model
      model_summary <- summary(model)
      
      # Add the model summary to the list
      models[[paste0("m", i)]] <- model_summary
    }
  }
  
  return(models)
}

income_cmest_1 <- function(x, num_imputations = 5) {
  models <- list()
  
  for (i in 1:num_imputations) {
    # Create the dataset variable name (imputed_dataset_i)
    dataset_var <- paste0("imputed_dataset_", i)
    
    # Perform cmest for each imputation
    model <- cmest(
      data = get(dataset_var), 
      model = "msm", outcome = x, 
      exposure = "income80", mediator = "income07", basec = c("age", "sex", "edu80"), 
      postc = "edu107", EMint = F,
      ereg = "logistic", yreg = "linear", mreg = list("logistic"),
      wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
      astar = 0, a = 1, mval = list(1), 
      estimation = "imputation", inference = "bootstrap", nboot=1000)
    
    if (!is.null(model)) {
      # Get the summary of the model
      model_summary <- summary(model)
      
      # Add the model summary to the list
      models[[paste0("m", i)]] <- model_summary
    }
  }
  
  return(models)
}


#This function extracts and pools the estimates when the outcome is ordinal

poolRR0<- function(x) {
  ##TOTAL EFFECT
  lonedu_TE<-filter(x, Parameter == "Rte") %>% select(Estimate)
  lonedu_TE <- lonedu_TE[["Estimate"]]
  lonedu_TEse<-filter(x, Parameter == "Rte") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_TEse <- lonedu_TEse[["variance"]]
  
  try<-pool.scalar(lonedu_TE, lonedu_TEse, n = 1,685)
  
  qbar <- try$qbar
  t <- try$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
  
  ##DIRECT EFFECT
  lonedu_DE<-filter(x, Parameter == "rRpnde") %>% select(Estimate)
  lonedu_DE <- lonedu_DE[["Estimate"]]
  lonedu_DEse<-filter(x, Parameter == "rRpnde") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_DEse <- lonedu_DEse[["variance"]]
  
  try1<-pool.scalar(lonedu_DE, lonedu_DEse, n = 1,685)
  
  qbar <- try1$qbar
  t <- try1$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
  
  ##INDIRECT EFFECT
  lonedu_IE<-filter(x, Parameter == "rRpnie") %>% select(Estimate)
  lonedu_IE <- lonedu_IE[["Estimate"]]
  lonedu_IEse<-filter(x, Parameter == "rRpnie") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_IEse <- lonedu_IEse[["variance"]]
  
  try2<- pool.scalar(lonedu_IE, lonedu_IEse, n = 1,685)
  
  qbar <- try2$qbar
  t <- try2$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
  
  ##PROPORTION MEDIATED
  lonedu_pm<-filter(x, Parameter == "pm") %>% select(Estimate)
  lonedu_pm <- lonedu_pm[["Estimate"]]
  lonedu_pmse<-filter(x, Parameter == "pm") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_pmse <- lonedu_pmse[["variance"]]
  
  try3<- pool.scalar(lonedu_pm, lonedu_pmse, n = 1,685)
  
  qbar <- try3$qbar
  t <- try3$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
}


#This function extracts and pools the estimates when the outcome is continuous

poolRR1<- function(x) {
  ##TOTAL EFFECT
  lonedu_TE<-filter(x, Parameter == "te") %>% select(Estimate)
  lonedu_TE <- lonedu_TE[["Estimate"]]
  lonedu_TEse<-filter(x, Parameter == "te") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_TEse <- lonedu_TEse[["variance"]]
  
  try<-pool.scalar(lonedu_TE, lonedu_TEse, n = 1,685)
  
  qbar <- try$qbar
  t <- try$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
  
  ##DIRECT EFFECT
  lonedu_DE<-filter(x, Parameter == "rpnde") %>% select(Estimate)
  lonedu_DE <- lonedu_DE[["Estimate"]]
  lonedu_DEse<-filter(x, Parameter == "rpnde") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_DEse <- lonedu_DEse[["variance"]]
  
  try1<-pool.scalar(lonedu_DE, lonedu_DEse, n = 1,685)
  
  qbar <- try1$qbar
  t <- try1$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
  
  ##INDIRECT EFFECT
  lonedu_IE<-filter(x, Parameter == "rpnie") %>% select(Estimate)
  lonedu_IE <- lonedu_IE[["Estimate"]]
  lonedu_IEse<-filter(x, Parameter == "rpnie") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_IEse <- lonedu_IEse[["variance"]]
  
  try2<- pool.scalar(lonedu_IE, lonedu_IEse, n = 1,685)
  
  qbar <- try2$qbar
  t <- try2$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
  
  ##PROPORTION MEDIATED
  lonedu_pm<-filter(x, Parameter == "pm") %>% select(Estimate)
  lonedu_pm <- lonedu_pm[["Estimate"]]
  lonedu_pmse<-filter(x, Parameter == "pm") %>% select(Std.error) %>%  mutate(variance=Std.error^2)
  lonedu_pmse <- lonedu_pmse[["variance"]]
  
  try3<- pool.scalar(lonedu_pm, lonedu_pmse, n = 1,685)
  
  qbar <- try3$qbar
  t <- try3$t
  
  # Critical value for a 95% confidence interval (two-tailed)
  critical_value <- qnorm(0.975)  # Since it's a z-distribution
  
  # Calculate margin of error
  margin_of_error <- critical_value * sqrt(t)
  
  # Calculate confidence interval bounds
  lower_bound <- qbar - margin_of_error
  upper_bound <- qbar + margin_of_error
  
  # Print confidence interval
  print(qbar)
  cat("Lower bound:", lower_bound, "\n")
  cat("Upper bound:", upper_bound, "\n")
}


#EXAMPLE: Loneliness x Income (Same procedure for all mediation models)
set.seed(5)

#Runs the models in the 5 imputed datasets 
lonel_inc<-income_cmest_0("lon_TILS")
print(lonel_inc)

# Extracts each model and adds a parameter column
loninc_dfs <- list()

for (i in 1:5) {
  loninc_dfs[[i]] <- lonel_inc[[i]][[15]]
  loninc_dfs[[i]]$Parameter <- rownames(loninc_dfs[[i]])
  rownames(loninc_dfs[[i]]) <- NULL
}

# Combine data frames vertically
all_loninc <- bind_rows(loninc_dfs, .id = "Imputation")

#pooled results
poolRR0(all_lonedu)

#### 4. Sensitivity analyses ###
#4.1. Multiple imputation without auxiliary variables

library(mice)
init = mice(data_MI, maxit=0) 
meth = init$method
predM = init$predictorMatrix

#id has not predictive value
predM[, c("id")]=0

meth[c("age", "loneliness", "social_network", "social_relations", "edu80")]="pmm" #numerical
meth[c("income80", "edu107", "income07")]="polr" #categorical ordered

set.seed(103)
imputedMI = mice(data_MI, method=meth, predictorMatrix=predM, m=5)

imputed_datasetsMI <- complete(imputedMI, action = "long", include = TRUE)

#From here, the same procedure as described in 1.3.

#4.2. Sensitivity analysis for unmeasured confounding (e-value) 

sens_analysis <- function(models_result) {
  sens_results <- list()
  
  for (i in seq_along(models_result)) {
    # Extract the model from the result
    model <- models_result[[i]]
    
    # Perform sensitivity analysis for unmeasured confounding
    sens_result <- cmsens(object = model, sens = "uc")
    
    # Add the sensitivity analysis result to the list
    sens_results[[paste0("sens", i)]] <- sens_result
  }
  
  return(sens_results)
}

# Example: Loneliness and education
lonedu_models<- edu_cmest_0sa("lon_TILS") #edu_cmest0sa is a function as the one used before but that storages in the final lists the models instead of the models' summaries 
lon_edu_sa<-sens_analysis(lonedu_models)
print(lon_edu_sa)

#4.3. Sensitivity analysis for the adult educational attainment measure
#educational attainment (primary and low secondary education as risk)
imp_data$edu107 <-recode(imp_data$edu107, `0` = 1, `1` = 0, `2` = 0, `3` = 0) #less that upper secondary school as at risk 

#Same CMA models as in the main analysis

### 5. CMA with interaction Exp x Med
#For models including the interaction, the procedure was the same, but the function sets EMint=T. Example:
edu_cmest_int0 <- function(x, num_imputations = 5) {
  models <- list()
  
  for (i in 1:num_imputations) {
    # Create the dataset variable name (imputed_dataset_i)
    dataset_var <- paste0("imputed_dataset_", i)
    
    # Perform cmest for each imputation
    model <- cmest(
      data = get(dataset_var), 
      model = "msm", outcome = x, 
      exposure = "edu80", mediator = "edu107", 
      basec = c("age", "sex"), postc = "income80", 
      EMint = T, ereg = "logistic", 
      yreg = "ordinal", mreg = list("logistic"),
      wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
      astar = 0, a = 1, mval = list(1), 
      estimation = "imputation", inference = "bootstrap", nboot=1000
    )
    
    if (!is.null(model)) {
      # Get the summary of the model
      model_summary <- summary(model)
      
      # Add the model summary to the list
      models[[paste0("m", i)]] <- model_summary
    }
  }
  
  return(models)
}

### 6. Supplementary Material: empirical check of the positivity assumption

#Following #https://simonejdemyr.com/r-tutorials/statistics/tutorial8.html 
#Example using imputed_dataset1

#Eductional attainment models 
# Fit a propensity score model for the exposure
ps_model_exposure1 <- glm(edu80 ~ age, data = imputed_dataset_1, family = binomial())

# Fit a propensity score model for the mediator
ps_model_mediator1 <- glm(edu107 ~ edu80 + age + income80 + sex, data = imputed_dataset_1, family = binomial())

# Create a new data frame that contains the PSs and the actual exposure/mediator status 
Ps_exp <- data.frame(ps = predict(ps_model_exposure1, type = "response"), exposure = ps_model_exposure1$model$edu80)
Ps_med <- data.frame(ps = predict(ps_model_mediator1, type = "response"), mediator = ps_model_mediator1$model$edu107)

head(Ps_exp, 40)
head(Ps_med, 40)

# Summary statistics for exposure propensity scores (rounded to 3 decimals)
exp_summary_0 <- summary(Ps_exp$ps[Ps_exp$exposure == 0])
exp_summary_1 <- summary(Ps_exp$ps[Ps_exp$exposure == 1])

print(round(exp_summary_0, 3))
print(round(exp_summary_1, 3))

# Summary statistics for mediator propensity scores (rounded to 3 decimals)
med_summary_0 <- summary(Ps_med$ps[Ps_med$mediator == 0])
med_summary_1 <- summary(Ps_med$ps[Ps_med$mediator == 1])

print(round(med_summary_0, 3))
print(round(med_summary_1, 3))

#Income models

# Fit a propensity score model for the exposure
ps_model_exposure2 <- glm(income80 ~ age  + edu80, data = imputed_dataset_1, family = binomial())

# Fit a propensity score model for the mediator
ps_model_mediator2 <- glm(income07 ~ income80 + age + edu80 + edu107 + sex, data = imputed_dataset_1, family = binomial())

# Create a new data frame that contains the PSs and the actual exposure/mediator status 
Ps_exp2 <- data.frame(ps = predict(ps_model_exposure2, type = "response"), exposure = ps_model_exposure2$model$income80)
Ps_med2 <- data.frame(ps = predict(ps_model_mediator2, type = "response"), mediator = ps_model_mediator2$model$income07)

head(Ps_exp2)
head(Ps_med2)

# Summary statistics for exposure propensity scores (rounded to 3 decimals)
exp2_summary_0 <- summary(Ps_exp2$ps[Ps_exp2$exposure == 0])
exp2_summary_1 <- summary(Ps_exp2$ps[Ps_exp2$exposure == 1])

print(round(exp2_summary_0, 3))
print(round(exp2_summary_1, 3))

# Summary statistics for mediator propensity scores (rounded to 3 decimals)
med2_summary_0 <- summary(Ps_med2$ps[Ps_med2$mediator == 0])
med2_summary_1 <- summary(Ps_med2$ps[Ps_med2$mediator == 1])

print(round(med2_summary_0, 3))
print(round(med2_summary_1, 3))
