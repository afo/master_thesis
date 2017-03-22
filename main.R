###### INSTALL & LOAD PACKAGES ######

install.packages("psych", dependencies = TRUE)
install.packages("mcmcplots")
install.packages("R2jags")
install.packages("Boruta")
install.packages("jagsUI")
install.packages("matrixStats")
install.packages("mlbench")
install.packages("caret")

library(psych)
library(mcmcplots)
library(R2jags)
library(Boruta)
library(jagsUI)
library(matrixStats)
library(mlbench)
library(caret)

###### SET WORKING DIRECTORY & SOURCE FUNCTIONS ######

set.seed(71)
# setwd("~/thesis/R/JAGS/final")

source(paste(getwd(), '/Functions/loadData.R', sep=""))
source(paste(getwd(), '/Functions/simulateData.R', sep=""))
source(paste(getwd(), '/Functions/posteriorPredictiveChecks.R', sep=""))
source(paste(getwd(), '/Functions/itemCharacteristicCurves.R', sep=""))
source(paste(getwd(), '/Functions/itemInformationCurves.R', sep=""))
source(paste(getwd(), '/Functions/plotCor.R', sep=""))
source(paste(getwd(), '/Functions/featureRanking.R', sep=""))
source(paste(getwd(), '/Functions/variableReduction.R', sep=""))


###### MAIN ######
#
# Authors: Alexander Fred Ojala & Johan Eng Larsson
#
# Note to run this code contact the authors, as you'll need an anonymized data set.
#
# This is the main script for the analysis performed in the Master Thesis
# "Construction of the Berkeley Innovation Index: A higher order item response
# theory model approach".
#
# Thesis URL: https://lup.lub.lu.se/student-papers/search/publication/8895254
#
# All results obtained in the analysis is presented in Chapter 5: Results and 
# analysis with the the exceptions of the parts "Exploratory Analysis" and 
# "Prior Analysis". For these two parts the result is presented in Chapter 4 of the 
# thesis. All model parameter estimates and plots are presented in the appendix.
# 
# With the two exceptions presented above the code also follows the same structure
# as Chapter 5. This means that the code for the analysis is presented in the same 
# order as the result is presented in the thesis. All results are obtained in code
# section "Result and Analysis" and each new section of the code is marked
# ###### Section Name ######.
#
# This sections are named as the sections of Chapter 5 and the result found in that
# section of the thesis is obtained in the section of the code with same name. 
#
# Results presented in the thesis, either as tables or plots, are marked with a 
# comment on what figure/table in the thesis that presents the obtained results.


data.full <- readData() # read in data


#----- MCMC simulations for HO-GRM and HO-GPCM -----#

# MCMC simulations to obtain model parameter estimates for all model 
# parameters. 
#
# Two different model estimations are used: jags.parallell ('R2Jags')
# and jags ('jagsUI'). The resulting parameter estimations are identical
# but the graphical covergence check is better when running 'R2Jags' while
# 'jagsUI' facilitates the data handling.
#
# Model parameters:
#   theta:  second order latent traits, [N x 1]
#   gamtes: first order latent traits, [N x d]
#   eps:    error parameter, [N x d]
#   lambda: factor loading, [d x 1]
#   a:      item discrimination, [J x 1]
#   b:      item difficulty, [J x (K-1)]

# Initialize MCMC parameters----
Y <- data.full                  # Item response matrix

N <- nrow(Y)                    # Number of samples/subjects
J <- ncol(Y)                    # Number of items
d <- 6                          # Number of domains
K <- max(Y)                     # Number of possible item responses
test <- sort(rep(1:d, J/d))     # Vector defining which domain an item belongs to


# 'R2jags' runs----
data <- list("Y" = Y, "N" = N, "J" = J, "cat"=d, "K" = K,"test" = test)
param <- c("a", "b", "lambda", "gamtes", "theta", "eps")

mod.grm.graph <- jags.parallel(data = data, inits = NULL, parameters.to.save = param, 
                          model.file = paste(getwd(),"/JAGSmodels/hogrm.txt", sep=""),
                          n.chains = 3, n.iter = 3, n.burnin = 1)

mod.gpcm.graph <- jags.parallel(data = data, inits = NULL, parameters.to.save = param, 
                          model.file = paste(getwd(),"/JAGSmodels/hogpcm.txt", sep=""),
                          n.chains = 3, n.iter = 3, n.burnin = 1)

# 'jagsUI' runs----
data <- list("Y" = Y, "N" = N, "J" = J, "cat"=d, "K" = K,"test" = test)
param <- c("a", "b", "lambda", "gamtes", "theta", "eps", "Ynew")

mod.grm <- jags(data = data, inits = NULL, parameters.to.save = param, 
                     model.file = paste(getwd(),"/JAGSmodels/hogrm_ppcm.txt", sep=""), 
                     n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE, codaOnly = c("Ynew"))

mod.gpcm <- jags(data = data, inits = NULL, parameters.to.save = param, 
                      model.file = paste(getwd(),"/JAGSmodels/hogpcm_ppcm.txt", sep=""), 
                      n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE, codaOnly = c("Ynew"))


#------  Results & Analysis  ------#
# The rest of the script covers the reults presented in Chapter 5 of the thesis.
# All new sections of Chapter 5 is here marked: ###### Section Name ######


###### Convergence Checks ######

# Convergence check through a graphical convergence check, conducted
# by manually evaluating trace plot, acf, posterior distribution and
# running mean (all obtained through 'mcmcplot'), and the Gelman-Rubin
# diagnostic (Rhat).

# Graphical convergence check----
mcmc.grm <- as.mcmc(mod.grm.graph)
mcmc.gpcm <- as.mcmc(mod.gpcm.graph)


mcmcplot(mcmc.grm) # Figure 5.1: HO-GRM
mcmcplot(mcmc.gpcm) # Figure 5.2: HO-GPCM                   

# Gelman-Rubin Diagnostics
param <- names(mod.gpcm$Rhat)[1:6]
Rhats <- matrix(nrow=2, ncol=6)
colnames(Rhats) <- param
row.names(Rhats) <- c("HO-GRM", "HO-GPCM")

for (i in 1:6){
  Rhats[1,i] <- max(unlist(mod.grm$Rhat[param[i]]))
  Rhats[2,i] <- max(unlist(mod.gpcm$Rhat[param[i]]))
}


Rhats # Table 5.1: G-R statistics                                


###### Model Selection (HO-GRM vs. HO-GPCM) ######

# Goodness of fit assessment through DIC, variance, CV and
# Posterior Predictive Check.

# DIC and variance----
DIC.variance <- matrix(nrow=2, ncol=8)
row.names(DIC.variance) <- c("HO-GRM", "HO-GPCM")

DIC.variance[1,1] <- mod.grm$DIC
DIC.variance[2,1] <- mod.gpcm$DIC

DIC.variance[1,2] <- mean(mod.grm$sd$theta^2)
DIC.variance[2,2] <- mean(mod.gpcm$sd$theta^2)

DIC.variance[1,3:8] <- colMeans(mod.grm$sd$gamtes^2)
DIC.variance[2,3:8] <- colMeans(mod.gpcm$sd$gamtes^2)

DIC.variance  # Table 5.2: DIC and variance                         


# Coefficient of variation----
CV <- matrix(nrow=2, ncol=7)
row.names(CV) <- c("HO-GRM", "HO-GPCM")

CV[1,1] <- mean(mod.grm$sd$theta^2)/mean(abs(mod.grm$mean$theta))
CV[2,1] <- mean(mod.gpcm$sd$theta^2)/mean(abs(mod.gpcm$mean$theta))

CV[1,2:7] <- colMeans(mod.grm$sd$gamtes^2)/colMeans(abs(mod.grm$mean$gamtes))
CV[2,2:7] <- colMeans(mod.gpcm$sd$gamtes^2)/colMeans(abs(mod.gpcm$mean$gamtes))

CV # Table 5.3: Coefficients of variation                                


# Posterior predicitve checks----
pmean.grm <- meanPPCM(mod.grm,Y)   # Figure 5.3(a)
pstd.grm <- stdPPCM(mod.grm,Y)     # Figure 5.3(b)

pmean.gpcm <- meanPPCM(mod.gpcm,Y) # Figure 5.4(a)
pstd.gpcm <- stdPPCM(mod.gpcm,Y)   # Figure 5.4(b)

PPCM <- matrix(c(pmean.grm, pstd.grm, pmean.gpcm, pstd.gpcm), nrow=2, ncol=6, byrow=TRUE)
row.names(PPCM) <- c("HO-GRM", "HO-GPCM")

PPCM # Table 5.4: PPC p-values


### Model Selection: Simulated data ###

# Two data sets are simulated: one with HO-GRM model structure and
# one with HO-GPCM model structure. Thereafter are the model parameters
# from both models estimated using the HO-GRM model and the HO-GPCM model.
# In total four MCMC-simualtions are done.
# For each model the DIC and the correlation between model parameter
# estimates and 'true' parameter value (from simulated data) is computed.
#
# The result are stored in the matrix modComp.result

# Simulate data----
grm.data <- simulateData()
gpcm.data <- simulateData(grm=FALSE)

# GRM data MCMC runs----
Y <- grm.data$Y                 # Item response matrix

N <- nrow(Y)                    # Number of samples/subjects
J <- ncol(Y)                    # Number of items
d <- length(grm.data$lambda)    # Number of domains
K <- max(Y)                     # Number of possible item responses
test <- sort(rep(1:d, J/d))     # Vector defining which domain an item belongs to

data <- list("Y" = Y, "N" = N, "J" = J, "cat"=d, "K" = K,"test" = test)
param <- c("a", "b", "lambda", "gamtes", "theta", "eps")

mod.grm.grmData <- jags(data = data, inits = NULL, parameters.to.save = param, 
                     model.file = paste(getwd(),"/JAGSmodels/hogrm.txt", sep=""), 
                     n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE)

mod.gpcm.grmData <- jags(data = data, inits = NULL, parameters.to.save = param, 
                      model.file = paste(getwd(),"/JAGSmodels/hogpcm.txt", sep=""), 
                      n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE)

# GPCM data MCMC runs----
Y <- gpcm.data$Y                # Item response matrix

N <- nrow(Y)                    # Number of samples/subjects
J <- ncol(Y)                    # Number of items
d <- length(gpcm.data$lambda)   # Number of domains
K <- max(Y)                     # Number of possible item responses
test <- sort(rep(1:d, J/d))     # Vector defining which domain an item belongs to

data <- list("Y" = Y, "N" = N, "J" = J, "cat"=d, "K" = K,"test" = test)
param <- c("a", "b", "lambda", "gamtes", "theta", "eps")

mod.grm.gpcmData <- jags(data = data, inits = NULL, parameters.to.save = param, 
                    model.file = paste(getwd(),"/JAGSmodels/hogrm.txt", sep=""), 
                    n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE)

mod.gpcm.gpcmData <- jags(data = data, inits = NULL, parameters.to.save = param, 
                     model.file = paste(getwd(),"/JAGSmodels/hogpcm.txt", sep=""), 
                     n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE)

# GRM data model results----
grmData.result <- matrix(nrow = 2, ncol = 8)

grmData.result[1,1] <- mod.grm.grmData$DIC
grmData.result[2,1] <- mod.gpcm.grmData$DIC

grmData.result[1,2] <- cor(mod.grm.grmData$mean$theta, grm.data$theta)
grmData.result[2,2] <- cor(mod.gpcm.grmData$mean$theta, grm.data$theta)

grmData.result[1,3:8] <- diag(cor(mod.grm.grmData$mean$gamtes, grm.data$gamtes))
grmData.result[2,3:8] <- diag(cor(mod.gpcm.grmData$mean$gamtes, grm.data$gamtes))

# GPCM data model results----
gpcmData.result <- matrix(nrow = 2, ncol = 8)

gpcmData.result[1,1] <- mod.grm.gpcmData$DIC
gpcmData.result[2,1] <- mod.gpcm.gpcmData$DIC

gpcmData.result[1,2] <- cor(mod.grm.gpcmData$mean$theta, gpcm.data$theta)
gpcmData.result[2,2] <- cor(mod.gpcm.gpcmData$mean$theta, gpcm.data$theta)

gpcmData.result[1,3:8] <- diag(cor(mod.grm.gpcmData$mean$gamtes, gpcm.data$gamtes))
gpcmData.result[2,3:8] <- diag(cor(mod.gpcm.gpcmData$mean$gamtes, gpcm.data$gamtes))

# Full result matrix: DIC and correlation----
modComp.result <- rbind(grmData.result, gpcmData.result)
row.names(modComp.result) <- c("HO-GRM (GRM)", "HO-GPCM (GRM)", "HO-GRM (GPCM)", "HO-GPCM (GPCM)")

modComp.result  # Table 5.5: Model fit comparison


###### HO-GRM Model Results ######

# Plot Item Characteristic Curves----
plotICC(mod.grm)  # Figure 5.5(a)-(d) & Figure B.1-3

# Plot Item Information Curves----
plotIIC(mod.grm)  # Figure 5.6 & Figure B.4


###### HO-GRM Ability Estimates ######

# Figure 5.7 & Figure B.5

# Plot trait distributions----
hist(mod.grm$mean$theta, xlab = expression(paste("Innovation, ", theta^(2))), main=expression("Distribution of Innovation: MCMC output"))

titles =  c("Trust", "Resilience", "Diversity", "Belief", "Perfection", "Collaboration")
for(i in 1:6){
  hist(mod.grm$mean$gamtes[,i], xlab = substitute(paste(titles, ", ", theta[i]^(1)), list(titles=titles[i])), main=substitute(paste("Distribution of ", titles, ": MCMC output"), list(titles=titles[i])))
}

# Lambda estimates----

# Table 5.6: Factor loadings
lambda.estimate <- mod.grm$mean$lambda 
lambda.sd <- mod.grm$sd$lambda        

# Correlation matrix----

theta <- mod.grm$mean$theta
gamtes <- mod.grm$mean$gamtes
traits <- cbind(theta, gamtes)

cor(traits) # Table 5.7: Trait correlation


# Example: 5 sample subjects (random samples s = 95, 629, 730)----
which.max(mod.grm$mean$theta)      # s = 375
which.min(mod.grm$mean$theta)      # s = 205

sample = c(95, 205, 375, 629, 730)

sample.estimates <- cbind(mod.grm$mean$theta[sample], mod.grm$mean$gamtes[sample,])
sample.sd <- cbind(mod.grm$sd$theta[sample], mod.grm$sd$gamtes[sample,])

row.names(sample.estimates) <- sample
row.names(sample.sd) <- sample

# Table 5.8
data.full[sample,]        
sample.estimates
sample.sd

# Plot correlations----
plotCorMean(data.full, mod.grm$mean$theta)   # Figure 5.8(a)
plotCorVar(data.full, mod.grm$mean$theta)    # Figure 5.8(b) 



###### Variable Reduction ######

# Data preperation for variable reduction----
Y <- data.full
names(Y) <- 1:24

theta<-mod.grm$mean$theta
gamtes<-mod.grm$mean$gamtes
traits <- cbind(gamtes, theta)

### Feature Ranking ###

# Boruta----

# Train models
boruta.models <- lapply(1:7, function(i)runBoruta(Y,traits[,i]))

# Compute total importance per item and obtain feature ranking
boruta.importance <- rowSums(sapply(1:7, function(i)colMeans(boruta.models[[i]]$ImpHistory)))[1:24]  # Figure 5.9
feature.ranking.boruta <- as.numeric(names(sort(boruta.importance, decreasing = TRUE)))

feature.ranking.boruta   # Section 5.4, presented as a vector

# GBM----

# Prepare training scheme and train models
control <- trainControl(method="repeatedcv", number=3, repeats=10)
gbm.models <- lapply(1:7, function(i)train(traits[,i]~., data=Y, method="gbm", 
                                           preProcess="scale", trControl=control))

# Compute total importance per item and obtain feature ranking
gbm.importance <- lapply(1:7, function(i)varImp(gbm.models[[i]], scale=FALSE))
gbm.importance.total <- rowSums(sapply(1:7, function(i)gbm.importance[[i]]$importance$Overall))
feature.ranking.gbm <- sort(gbm.importance.total, decreasing = TRUE, index.return = TRUE)$ix

feature.ranking.gbm   # Section 5.4, presented as a vector

### Variable reduction analysis ###

n <- 500        # number of cross-validations in each linear regression
a <- 0.9        # data split in cross-validation

getGoodnessOfFitMatrix(Y, traits, feature.ranking.boruta, n, a)  # Table 5.9: Boruta, goodness of fit
getGoodnessOfFitMatrix(Y, traits, feature.ranking.gbm, n, a)     # Table 5.9: GBM, goodness of fit 

reduced.model <- sort(feature.ranking.gbm[1:17])  # Table 5.10: Reduced model


###### Final Algorithm ######

# Linear regression on full and reduced model----
n <- 1000
a <- 0.9

lm.full.models <- lapply(1:7, function(i)lmCrossValidationRun(as.matrix(Y), traits[,i], n, a))
lm.reduced.models <- lapply(1:7, function(i)lmCrossValidationRun(as.matrix(Y[,reduced.model]), traits[,i], n, a))

# Get regression parameters----
beta.full <- sapply(1:7, function(i)lm.full.models[[i]]$mod$coefficient)
beta.reduced <- sapply(1:7, function(i)lm.reduced.models[[i]]$mod$coefficient)  # Table 5.11: Regression parameters

# Compute latent trait estimates and RMSE from linear regression parameters----
lm.estimates.full.model <- sapply(1:7, function(i)getYhat(beta.full[,i], Y))
lm.estimates.reduced.model <- sapply(1:7, function(i)getYhat(beta.reduced[,i], Y[,reduced.model]))

RMSE.full <- sapply(1:7, function(i)getRMSE(lm.estimates.full.model[,i], traits[,i]))        # Table 5.12: RMSE_full
RMSE.reduced <- sapply(1:7, function(i)getRMSE(lm.estimates.reduced.model[,i], traits[,i]))  # Table 5.12. RMSE_red

# Get max and min values used for scaling----
lo.values <- sapply(1:7, function(i)getYhat(beta.full[,i], t(rep(1,dim(Y)[2]))))  # Table 5.13: X_full, min
hi.values <- sapply(1:7, function(i)getYhat(beta.full[,i], t(rep(5,dim(Y)[2]))))  # Table 5.13: X_full, max

lo.values.reduced <- sapply(1:7, function(i)getYhat(beta.reduced[,i], t(rep(1,dim(Y[,reduced.model])[2]))))  # Table 5.13: X_red, min
hi.values.reduced <- sapply(1:7, function(i)getYhat(beta.reduced[,i], t(rep(5,dim(Y[,reduced.model])[2]))))  # Table 5.13: X_red, max

# Cast final scores casted on I[1,10]----
final.scores.full <- sapply(1:7, function(i)oneToTenScale(lm.estimates.full.model[,i], lo.values[i], hi.values[i]))  
final.scores.reduced <- sapply(1:7, function(i)oneToTenScale(lm.estimates.reduced.model[,i], lo.values.reduced[i], hi.values.reduced[i]))

# Histograms of final scores----
titles =  c("Trust", "Resilience", "Diversity", "Belief", "Perfection", "Collaboration")

# Figure 5.10(a) & Figure B.6: Full model
for(i in 1:6){
  hist(final.scores.full[,i], 
       xlab = substitute(paste(titles, ", ", theta[i]^(1)), list(titles=titles[i])), 
       main=substitute(paste("Distribution of ", titles, ": Full Model"), list(titles=titles[i])), 
       xlim=c(1,10))
  axis(1,c(1:10))
}
hist(final.scores.full[,7], 
     xlab = expression(paste("Innovation, ", theta^(2))), 
     main=expression("Distribution of Innovation: Full Model"), 
     xlim=c(1,10))
axis(1,c(1:10))

# Figure 5.10(b) & Figure B.7: Reduced model
for(i in 1:6){
  hist(final.scores.reduced[,i], 
       xlab = substitute(paste(titles, ", ", theta[i]^(1)), list(titles=titles[i])), 
       main=substitute(paste("Distribution of ", titles, ": Reduced Model"), list(titles=titles[i])), 
       xlim=c(1,10))
  axis(1,c(1:10))
}
hist(final.scores.reduced[,7], 
     xlab = expression(paste("Innovation, ", theta^(2))),
     main=expression("Distribution of Innovation: Reduced Model"), 
     xlim=c(1,10))
axis(1,c(1:10))


###### OUTLIER ANALYSIS ######

# Leave-one-out outlier analysis conducted by running the MCMC simulation 
# on two reduced data sets and comparing the result with the results
# obtained from the full data set.
# In the first model the sample with the highest mean item response
# was removed, in the second the sample with the lowest.
#
# The models are compared using DIC and variancs. Result is stored in
# the matrix 'outlier.result'.

# Create reduced data sets
Y.hi <- Y[-which.max(rowMeans(Y)),]
Y.lo <- Y[-which.min(rowMeans(Y)),]

# Initialize model parameters (all but Y same for both models)
N <- nrow(Y.hi)                 # Number of samples/subjects
J <- ncol(Y.hi)                 # Number of items
d <- 6                          # Number of domains
K <- max(Y.hi)                  # Number of possible item responses
test <- sort(rep(1:d, J/d))     # Vector defining which domain an item belongs to

data.hi <- list("Y" = Y.hi, "N" = N, "J" = J, "cat"=d, "K" = K,"test" = test)
data.lo <- list("Y" = Y.lo, "N" = N, "J" = J, "cat"=d, "K" = K,"test" = test)
param <- c("a", "b", "lambda", "gamtes", "theta", "eps")

# Run MCMC
mod.outlier.hi <- jags(data = data.hi, inits = NULL, parameters.to.save = param, 
                    model.file = paste(getwd(),"/JAGSmodels/hogrm.txt", sep=""), 
                    n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE)

mod.outlier.lo <- jags(data = data.lo, inits = NULL, parameters.to.save = param, 
                     model.file = paste(getwd(),"/JAGSmodels/hogrm.txt", sep=""), 
                     n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE)

# Result: DIC and variance
outlier.result <- matrix(nrow=3, ncol=8)

outlier.result[1,1] <- mod.grm$DIC
outlier.result[2,1] <- mod.outlier.lo$DIC
outlier.result[3,1] <- mod.outlier.hi$DIC

outlier.result[1,2] <- mean(mod.grm$sd$theta^2)
outlier.result[2,2] <- mean(mod.outlier.lo$sd$theta^2)
outlier.result[3,2] <- mean(mod.outlier.hi$sd$theta^2)

outlier.result[1,3:8] <- colMeans(mod.grm$sd$gamtes^2)
outlier.result[2,3:8] <- colMeans(mod.outlier.lo$sd$gamtes^2)
outlier.result[3,3:8] <- colMeans(mod.outlier.hi$sd$gamtes^2)

row.names(outlier.result) <- c("Full", "Lo", "Hi")

outlier.result  # Table 5.14: Outlier analysis



##### STRUCTURE ANALYSIS ######

# Simulate data according to model structures given in model.spec. Run MCMC
# for each model and compute the correlation between estimated values (from 
# MCMC) and 'true' values given by simulated data set.
#
# Correlations are stored in the matrix 'structure.result'.


# Define model structure specifications
N <- c(rep(1000,8), 3000)                   # subjects
J <- c(18, 18, 24, 24, 24, 36, 36, 36, 24)  # items
d <- c(3, 6, 3, 6, 8, 3, 6, 9, 6)           # domains

model.spec <- cbind(N,J,d)

# Simulate data, run MCMC, compute correlations.

mod.structure <- vector(mode = "list", length = length(N))
structure.result <- matrix(nrow = 9, ncol = 10)

for(i in 1:dim(model.spec)[1]){
  spec <- model.spec[i,]
  sim.data <- simulateData(n = spec["N"], m = spec["J"], d = spec["d"])
  test <- sort(rep(1:spec["d"], spec["J"]/spec["d"]))
  
  data <- list("Y" = sim.data$Y, "N" = spec["N"], "J" = spec["J"], "cat"= spec["d"], 
               "K" = 5,"test" = test)
  param <- c("a", "b", "lambda", "gamtes", "theta", "eps")
  
  mod.structure[[i]] <- jags(data = data, inits = NULL, parameters.to.save = param, 
                              model.file = paste(getwd(),"/JAGSmodels/hogrm.txt", sep=""), 
                              n.chains = 3, n.iter = 3, n.burnin = 1, parallel = TRUE)
  structure.result[i,1] <- cor(mod.structure[[i]]$mean$theta, sim.data$theta)
  structure.result[i,2:(spec["d"]+1)] <- diag(cor(mod.structure[[i]]$mean$gamtes, sim.data$gamtes))
}

  # Result: correlations
colnames(structure.result) <- c("2nd order", "1st (1)", "1st (2)", "1st (3)",
                                 "1st (4)", "1st (5)", "1st (6)", "1st (7)", "1st (8)", "1st (9)")

structure.result  # Table 5.15: Structure analysis

#----- Exploratory analysis -----#

# Exploratory analysis through Exploratory Factor Analysis, EFA ('omega')
# and Principal Component Analysis, PCA ('scree').
#
# Results are in thesis presented in Chapter 4.5.1: Exploratory Analysis

scree(data.full, factors=FALSE)  # Figure 4.5: PCA
omega(data.full, 6, sl=FALSE)     # Figure 4.6: EFA

#------ Prior Analysis ------#

# Simulate one HO-GRM data set. Run 12 different MCMC models where in each
# model one model parameter prior has been adjusted. For each model compute
# the correlation between estimated values (from MCMC) and 'true' values 
# given by simulated data set.
#
# Correlations are stored in the matrix 'prior.result'.
#
# Result is in thesis presented in Chapter 4.3.3: Prior selection

# Simulate data set
sim.data <- simulateData(seed=10)

# Initialize model parameters
Y <- sim.data$Y

N <- nrow(Y)                    # Number of samples/subjects
J <- ncol(Y)                    # Number of items
d <- 6                          # Number of domains
K <- max(Y)                     # Number of possible item responses
test <- sort(rep(1:d, J/d))     # Vector defining which domain an item belongs to

data <- list("Y" = Y, "N" = N, "J" = J, "cat"=d, "K" = K,"test" = test)
param <- c("a", "b", "lambda", "gamtes", "theta", "eps")

# Define path names to JAGS models
path.name <- paste(getwd(), "/JAGSmodels/hogrm",sep = "")

file.paths <- c(paste(path.name, ".txt", sep=""), paste(path.name, "_asdlo.txt", sep=""),
               paste(path.name, "_asdhi.txt", sep=""), paste(path.name, "_amuhi.txt", sep=""),
               paste(path.name, "_bsdlo.txt", sep=""), paste(path.name, "_bsdhi.txt", sep=""),
               paste(path.name, "_bmulo.txt", sep=""), paste(path.name, "_bmuhi.txt", sep=""),
               paste(path.name, "_lsdlo.txt", sep=""), paste(path.name, "_lmuhi.txt", sep=""),
               paste(path.name, "_tsdlo.txt", sep=""), paste(path.name, "_tsdhi.txt", sep=""))

# Run MCMC
mod.prior <- lapply(1:length(mod.names), 
                    function(i)jags(data = data, inits = NULL, parameters.to.save = param, 
                                    model.file = file.paths[i], n.chains = 3, n.iter = 3, 
                                    n.burnin = 1, parallel = TRUE))

# Result: correlations
prior.result.second.order <- sapply(1:length(mod.names), 
                                    function(i)cor(mod.prior[[i]]$mean$theta, sim.data$theta))
prior.result.first.order <- sapply(1:length(mod.names), 
                                   function(i)diag(cor(mod.prior[[i]]$mean$gamtes, sim.data$gamtes)))

prior.result <- t(rbind(prior.result.second.order, prior.result.first.order))
colnames(prior.result) <- c("2nd order", "1st (1)", "1st (2)", "1st (3)", 
                            "1st (4)", "1st (5)", "1st (6)")


prior.result  # Table 4.2: Prior analysis



