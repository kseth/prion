library("rJava")
.jinit(".")
computationClass <- J("ComputationUtils")
source("PrionModeling.R")
source("PlotFunctions.R")

lambda <- 4400
delta_m <- 5
beta <- 0.3
delta_p <- 0.04
b <- 0.001
polymerThreshold <- 6
polymerLengths0_size <- c(6)
polymerLengths0_count <- c(1)
monomers0 <- lambda/delta_m

endTime <- 10
numiterations <- 1000
multipleIterationResults <- prionMultipleIterations(endTime, numiterations, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

dev.new()
par(mfrow = c(2, 2))
plotHistogramMultipleIterations(multipleIterationResults)

timeSteps <- 1:1500/10
numiterations <- 100
timeStepResults <- prionTimeStepsMultipleIterations(timeSteps, 100, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

bHigh <- 0.01
timeStepResultsBHigh <- prionTimeStepsMultipleIterations(timeSteps, 100, lambda, delta_m, beta, delta_p, bHigh, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

bLow <- 0.0001
timeStepResultsBLow <- prionTimeStepsMultipleIterations(timeSteps, 100, lambda, delta_m, beta, delta_p, bLow, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

## Figure 2
dev.new()
par(mfrow = c(3, 4))
plotTimeStepsMultipleIterations(timeSteps, timeStepResults)
plotTimeStepsMultipleIterations(timeSteps, timeStepResultsBHigh)
plotTimeStepsMultipleIterations(timeSteps, timeStepResultsBLow)

endTime <- 150
numiterations <- 1000
bThreshold <- 0.00001
multipleIterationResultsBThreshold <- prionMultipleIterations(endTime, numiterations, lambda, delta_m, beta, delta_p, bThreshold, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

## Figure 3
dev.new()
par(mfrow = c(3, 1))
plotHistogramMultipleIterations(multipleIterationResultsBThreshold)

timeSteps <- 1:1500/10
numiterations <- 100
delta_pHigh <- 0.4
timeStepResultsPHigh <- prionTimeStepsMultipleIterations(timeSteps, 100, lambda, delta_m, beta, delta_pHigh, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

betaLow <- 0.03
timeStepResultsBetaLow <- prionTimeStepsMultipleIterations(timeSteps, 100, lambda, delta_m, betaLow, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass)

## Figure 4
dev.new()
par(mfrow = c(3, 4))
plotTimeStepsMultipleIterations(timeSteps, timeStepResultsPHigh)
plotTimeStepsMultipleIterations(timeSteps, timeStepResultsBetaLow)