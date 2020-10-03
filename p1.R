rm(list = ls())
gc()

# Load packages

if (!require(rstudioapi)) {
  install.packages("rstudioapi")
  library(rstudioapi)
}
if (!require(colorspace)) {
  install.packages("colorspace")
  library(colorspace)
}
if (!require(recommenderlab)) {
  install.packages("recommenderlab")
  library(recommenderlab)
}

# Set working directory

current.path = getActiveDocumentContext()$path 
setwd(dirname(current.path))

# Setting

M = 1000        # Number of consumers
N = 100         # Number of products
Nc = 2          # Number of products in the choice set
P = 10          # Number of purchases recorded in ratings matrix in cold start
Q = 10          # Number of purcheses recorded in ratings matrix by RS
Nr = 1          # Number of products recommended
Kc = 2          # Number of latent factors in consumers preferences
Krs.max = 4     # Maximum number of latent factors used by RS
Etheta = c(0, 0)      # Parameters of the distribution of preferences factors
Vtheta = c(2, 0.5)    # NB: Kc = 2!
Ebeta = c(0, 0)
Vbeta = c(2, 0.5)

# Set seed

set.seed(1)

# Draw population of consumers

Theta = matrix(data = NA, nrow = M, ncol = Kc)
Beta = matrix(data = NA, nrow = Kc, ncol = N)
for (k in 1:Kc) {
  Theta[, k] = rnorm(n = M, mean = Etheta[k], sd = sqrt(Vtheta[k]))
  Beta[k, ] = rnorm(n = N, mean = Ebeta[k], sd = sqrt(Vbeta[k]))
}
Eu = Theta %*% Beta         # "Expected" utility matrix

S = 1000000
epsilon = -log(-log(runif(n = S)))
stdBeta1 = rnorm(n = S, sd = sqrt(Vbeta[1]))
stdBeta2 = rnorm(n = S, sd = sqrt(Vbeta[2]))

CP = matrix(data = NA, nrow = M, ncol = N)        # Choice probabilities matrix
U0 = matrix(data = NA, nrow = M, ncol = N)        # Utility matrix (random draw)
RA = matrix(data = NA, nrow = M, ncol = N)        # Absolute rating matrix (based on U0)
EA = matrix(data = NA, nrow = M, ncol = N)        # Absolute rating matrix (based on Eu)

for (m in 1:M) {
  CP[m, ] = exp(Eu[m, ])/sum(exp(Eu[m, ]))
  
  Usim = Theta[m, 1]*stdBeta1+Theta[m, 2]*stdBeta2+epsilon
  qUsim = quantile(Usim, probs = c(0.2, 0.4, 0.6, 0.8))
  for (n in 1:N) {
    U0[m, n] = Eu[m, n]-log(-log(runif(1)))
    if (U0[m, n] <= qUsim[1]) {
      RA[m, n] = 1
    } else if (U0[m, n] < qUsim[2]) {
      RA[m, n] = 2
    } else if (U0[m, n] < qUsim[3]) {
      RA[m, n] = 3
    } else if (U0[m, n] < qUsim[4]) {
      RA[m, n] = 4
    } else {
      RA[m, n] = 5
    }
  }
  
  Usim = Theta[m, 1]*stdBeta1+Theta[m, 2]*stdBeta2
  qUsim = quantile(Usim, probs = c(0.2, 0.4, 0.6, 0.8))
  for (n in 1:N) {
    if (Eu[m, n] <= qUsim[1]) {
      EA[m, n] = 1
    } else if (Eu[m, n] < qUsim[2]) {
      EA[m, n] = 2
    } else if (Eu[m, n] < qUsim[3]) {
      EA[m, n] = 3
    } else if (Eu[m, n] < qUsim[4]) {
      EA[m, n] = 4
    } else {
      EA[m, n] = 5
    }
  }
    
}

# Generate ratings matrix

iProd = matrix(data = NA, nrow = M, ncol = P)
RVSparse = matrix(data = NA, nrow = M, ncol = N)        # Rating matrix come dice VD
RI = matrix(data = 0, nrow = M, ncol = N)               # Implicit rating matrix

# Cold start: Each consumer makes P choices

for (m in 1:M) {
  for (p in 1:P) {
    if (p == 1) {
      prods = sample(x = 1:N, size = Nc, replace = FALSE)
    } else {
      prods[1] = iProd[m, p-1]
      prods[2:Nc] = sample(x = (1:N)[-iProd[m, p-1]], size = Nc-1, replace = FALSE)
    }
    Um = Eu[m, prods]-log(-log(runif(n = Nc)))
    iProd[m, p] = prods[which.max(Um)]
    if (is.na(RVSparse[m, iProd[m, p]])) {
      RVSparse[m, iProd[m, p]] = 1
    } else {
      RVSparse[m, iProd[m, p]] = RVSparse[m, iProd[m, p]]+1
    }
    RI[m, iProd[m, p]] = 1
  }
  RVSparse[m, ] = RVSparse[m, ]/P
}

# Sparse matrices

EuSparse = matrix(data = NA, nrow = M, ncol = N)
CPSparse = matrix(data = NA, nrow = M, ncol = N)
U0Sparse = matrix(data = NA, nrow = M, ncol = N)
RASparse = matrix(data = NA, nrow = M, ncol = N)
EASparse = matrix(data = NA, nrow = M, ncol = N)
for (m in 1:M) {
  for (n in 1:N) {
    if (!is.na(RVSparse[m, n])) {
      EuSparse[m, n] = Eu[m, n]
      CPSparse[m, n] = CP[m, n]
      U0Sparse[m, n] = U0[m, n]
      RASparse[m, n] = RA[m, n]
      EASparse[m, n] = EA[m, n]
    }
  }
}

# Save parameters and ratings matrices

save(M, N, Nc, P, Q, Nr, Kc, Krs.max, Etheta, Vtheta, Ebeta, Vbeta, 
     Theta, Beta, Eu, U0, CP, RA, EA, RI,
     iProd, EuSparse, U0Sparse, CPSparse, RASparse, EASparse, RVSparse,
     file = "ParsRatings.Rdata")

# Regularized SVD decompositions

source("p1Functions.R")

# Eu matrix (expected utility, full M x N matrix)

R = Eu
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "Eu.Rdata")

# EuSparse matrix (sparse expected utility, M x N matrix)

R = EuSparse
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "EuSparse.Rdata")

# CP matrix (choice probabilities, full M x N matrix)

R = CP
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "CP.Rdata")

# CPSparse matrix (sparse choice probabilities, M x N matrix)

R = CPSparse
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "CPSparse.Rdata")

# U0 matrix (random utilities, full M x N matrix)

R = U0
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "U0.Rdata")

# U0Sparse matrix (sparse random utilities, M x N matrix)

R = U0Sparse
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "U0Sparse.Rdata")

# EA matrix (absolute ratings based on Eu, full M x N matrix)

R = EA
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "EA.Rdata")

# EASparse matrix (absolute ratings based on Eu, sparse M x N matrix)

R = EASparse
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "EASparse.Rdata")

# RA matrix (absolute ratings based on U0, full M x N matrix)

R = RA
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "RA.Rdata")

# RASparse matrix (absolute ratings based on U0, M x N matrix)

R = RASparse
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "RASparse.Rdata")

# RVSparse matrix (observed choice frequency, M x N matrix)

R = RVSparse
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "RVSparse.Rdata")

# RI matrix (implicit ratings full M x N matrix)

R = RI
X.RS = FitRS(Krs.max, R)
RMSEPlot(Krs.max, X.RS, R)

X = Eu
RankCorrelationPlots(Krs.max, X.RS, X, "spearman")
RankCorrelationPlots(Krs.max, X.RS, X, "kendall")
RankBoxPlots(Krs.max, X.RS, X)
RankStackedPlots(Krs.max, X.RS, X, 5)

save(X.RS, R, X, file = "RI.Rdata")

