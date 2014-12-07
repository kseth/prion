plotTimeStepsMultipleIterations <- function(timeSteps, timeStepResults) {
	plot(timeSteps, timeStepResults$mean[, 1], ylim = c(0, max(timeStepResults$mean[, 1]+timeStepResults$sd[, 1])), type = "l", xlab = "Time (days)", ylab = "Free Monomers")
	lines(timeSteps, timeStepResults$mean[, 1]+timeStepResults$sd[, 1], lty = 2)
	lines(timeSteps, timeStepResults$mean[, 1]-timeStepResults$sd[, 1], lty = 2)
	plot(timeSteps, timeStepResults$mean[, 2], ylim = c(0, max(timeStepResults$mean[, 2]+timeStepResults$sd[, 2])), type = "l", xlab = "Time (days)", ylab = "Prion Polymers")
	lines(timeSteps, timeStepResults$mean[, 2]+timeStepResults$sd[, 2], lty = 2)
	lines(timeSteps, timeStepResults$mean[, 2]-timeStepResults$sd[, 2], lty = 2)
	plot(timeSteps, timeStepResults$mean[, 3], ylim = c(0, max(timeStepResults$mean[, 3]+timeStepResults$sd[, 3])), type = "l", xlab = "Time (days)", ylab = "Bound Monomers")
	lines(timeSteps, timeStepResults$mean[, 3]+timeStepResults$sd[, 3], lty = 2)
	lines(timeSteps, timeStepResults$mean[, 3]-timeStepResults$sd[, 3], lty = 2)
	plot(timeSteps, timeStepResults$mean[, 4], ylim = c(0, 1), type = "l", xlab = "Time (days)", ylab = "Fraction Bound Monomers")
	lines(timeSteps, timeStepResults$mean[, 4]+timeStepResults$sd[, 4], lty = 2)
	lines(timeSteps, timeStepResults$mean[, 4]-timeStepResults$sd[, 4], lty = 2)
}

plotHistogramMultipleIterations <- function(multipleIterationResults) {
	hist(multipleIterationResults[, 1], xlab = "Free Monomers", main = "")
	hist(multipleIterationResults[, 3], xlab = "Bound Monomers", main = "")
	hist(multipleIterationResults[, 3]/apply(multipleIterationResults[, c(1, 3)], 1, sum), xlab = "Fraction Bound Monomers", main = "")
}