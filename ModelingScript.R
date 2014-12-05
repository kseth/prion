library("rJava")
.jinit(".")
computationClass <- J("ComputationUtils")
source("PrionModeling.R")

endTime <- 100
numiterations <- 10000
lambda <- 4400
delta_m <- 5
beta <- 0.3
delta_p <- 0.04
b <- 0.001
polymerThreshold <- 6
polymerLengths0_size <- c(6)
polymerLengths0_count <- c(1)
monomers0 <- 0

sampledValues <- prionMultipleIterations(endTime, numiterations, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

par(mfrow = c(2, 2))
hist(sampledValues[, 1])
hist(sampledValues[, 2])
hist(sampledValues[, 3])
hist(sampledValues[, 3]/sum(sampledValues[, c(1, 3)]))