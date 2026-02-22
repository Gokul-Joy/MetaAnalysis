
library(readxl)

my_data <- read_excel("GC_DFS_NMA.xlsx",na="NA")


colSums(is.na(my_data))
#===========================================

#install.packages("gemtc")
#install.packages("rjags")
library(gemtc)
library(rjags)

my_data$staging_covariate <- as.numeric(my_data$staging_covariate)

# Only log-transform non-reference arms, keep reference arms as NA
my_data$diff <- ifelse(is.na(my_data$std.err), NA, log(my_data$`HR Value`))
# Create the network object
# Note: 'std.err' and 'diff(LnHR)' are your column names from the image
# Create study-level covariate (one row per study)

study_covariates <- data.frame(
  study = c("ACTS-GC", "JACCRO GC-07", "ARTIST 2", "PRODIGY", 
            "RESONCE", "Zhao et al.", "POST", "Liu et al."),
  staging_covariate = c(0, 0, 0, 1, 1, 1, 0, 0)
)

# Rebuild network with covariate
network <- mtc.network(
  data.re = my_data,
  studies = study_covariates,
  description = "Gastric Cancer NMR"
)

plot(network)

# Verify covariate is in the network
network$studies


model <- mtc.model(network, 
                   type = "regression",
                   likelihood = "normal",
                   link = "identity",
                   regressor = list(
                     coefficient = "shared", 
                     variable = "staging_covariate",
                     control = "S_1"
                   ))

network$treatments


#Run Model B ( regression model)
results_regression <- mtc.run(model, n.adapt = 5000, n.iter = 20000, thin = 10)
summary(results_regression)


# Run Model A (standard NMA without covariate, for comparison):
model_standard <- mtc.model(network, 
                            type = "consistency",
                            likelihood = "normal",
                            link = "identity")

results_standard <- mtc.run(model_standard, n.adapt = 5000, n.iter = 20000, thin = 10)
summary(results_standard)

#Check convergence
plot(results_regression)       # trace plots — chains should mix well
gelman.plot(results_regression) # all values should be close to 1.0


#Compare DIC between models:
results_standard$deviance$DIC
results_regression$deviance$DIC


#Get adjusted treatment effects:
forest(relative.effect(results_regression, t1 = "S_1"))

#==============================================================================
#"How would all treatments perform in a purely pathological-staged population?"
#Adjusted Estimates at Covariate = 0 (Pathological Staging):
# Get effects adjusted to covariate = 0 (pathological staging population)
adj_effects <- relative.effect(results_regression, 
                               t1 = "S_1",
                               covariate = 0)
summary(adj_effects)
forest(adj_effects)



#===================
# Unadjusted (from standard NMA)
unadj <- relative.effect(results_standard, t1 = "S_1")

# Adjusted to covariate = 0
adj <- relative.effect(results_regression, t1 = "S_1", covariate = 0)

summary(unadj)
summary(adj)
#===================

# Extract β and heterogeneity (τ²) for your results table
# Beta coefficient
regression_summary <- summary(results_regression)
regression_summary$summaries$statistics["B",]  # mean, SD, SE

# Tau (sd.d) - between study heterogeneity
# Model A
regression_summary$summaries$statistics["sd.d",]

# Model B  
std_summary <- summary(results_standard)
std_summary$summaries$statistics["sd.d",]

# tau-squared
tau_regression <- regression_summary$summaries$statistics["sd.d","Mean"]^2
tau_standard <- std_summary$summaries$statistics["sd.d","Mean"]^2

cat("tau² Standard NMA:", tau_standard, "\n")
cat("tau² Regression NMA:", tau_regression, "\n")
