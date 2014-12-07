prionMultipleIterations <- function(endTime, numIterations, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass) {
	
	polymerDistribution0_innerhm <- convertIntegerVectorsHM(polymerLengths0_size, polymerLengths0_count, computationClass)
	sampledValues <- computationClass$prionGillespieInfectiousNonInfectious(as.numeric(endTime), as.integer(numIterations), as.numeric(lambda), as.numeric(delta_m), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))


	return(.jevalArray(sampledValues, simplify=TRUE))
}

prionTimeSteps <- function(sampleTimes, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass) {

	sampleTimes <- sort(sampleTimes)
	polymerDistribution0_innerhm <- convertIntegerVectorsHM(polymerLengths0_size, polymerLengths0_count, computationClass)
	sampledValues <- computationClass$prionGillespieInfectiousNoninfectious(.jarray(as.numeric(sampleTimes)), as.numeric(lambda), as.numeric(delta_m), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))

	return(.jevalArray(sampledValues, simplify=TRUE))
}

prionTimeStepsMultipleIterations <- function(sampleTimes, numIterations, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass) {

	sampleTimes <- sort(sampleTimes)
	polymerDistribution0_innerhm <- convertIntegerVectorsHM(polymerLengths0_size, polymerLengths0_count, computationClass)
	
	sampledValues <- computationClass$prionGillespieInfectiousNoninfectious(.jarray(as.numeric(sampleTimes)), as.numeric(lambda), as.numeric(delta_m), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))

	currentTotal <- as.matrix(.jevalArray(sampledValues, simplify=TRUE))
	currentTotal <- cbind(currentTotal, currentTotal[, 3]/apply(currentTotal[, c(1, 3)], 1, sum))
	currentTotalSquared <- currentTotal * currentTotal

	if(numIterations > 1) {	
		for(iter in 2:numIterations) {
			sampledValues <- computationClass$prionGillespieInfectiousNoninfectious(.jarray(as.numeric(sampleTimes)), as.numeric(lambda), as.numeric(delta_m), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))

			eval <- as.matrix(.jevalArray(sampledValues, simplify=TRUE))
			eval <- cbind(eval, eval[, 3]/apply(eval[, c(1, 3)], 1, sum))
			currentTotal <- currentTotal + eval
			currentTotalSquared <- currentTotalSquared + eval * eval
		}
	}

	overallMean <- currentTotal / numIterations
	overallMeanSquare <- currentTotalSquared / numIterations
	overallVariance <- overallMeanSquare - overallMean * overallMean

	return(list(mean=overallMean, var=overallVariance, sd = sqrt(overallVariance)))
}

prionPolymerDistribution <- function(endTime, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass) {

	polymerDistribution0_innerhm <- convertIntegerVectorsHM(polymerLengths0_size, polymerLengths0_count, computationClass)
	polymerDistributionFinal_innerhm <- computationClass$prionGillespiePolymerLengths(as.numeric(endTime), as.numeric(lambda), as.numeric(delta_m), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))

	return(convertIntegerHMVectors(polymerDistributionFinal_innerhm))
}

convertIntegerVectorsHM <- function(keys, values, computationClass) {
	if(length(keys) != length(values))
		stop("keys and values must have same length")

	innerhm <- computationClass$createEmptyInnerHashMap()
	for(i in 1:length(keys)){
		innerhm$put(as.integer(keys[i]), as.integer(values[i]))
	}

	return(innerhm)
}

convertIntegerHMVectors <- function(innerhm) {
	hm <- J(innerhm, "getHashMap")
	entrySet <- J(hm, "entrySet")
	iterator <- J(entrySet, "iterator")
	size <- J(entrySet, "size")
	keys <- rep(0, size)
	values <- rep(0, size)
	for(index in 1:size) {
		entry <- J(iterator, "next")
		keys[index] <- J(entry, "getKey")
		values[index] <- J(entry, "getValue")
	}

	return(list(keys=keys, values=values))
}