prion.multiple.iterations <- function(endTime, numIterations, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass) {
	
	polymerDistribution0_innerhm <- convert.integer.vectors.hm(polymerLengths0_size, polymerLengths0_count)
	sampledValues <- computationClass$prionGillespieInfectiousNonInfectious(as.numeric(endTime), as.integer(numIterations), as.numeric(lambda), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))


	return(.jevalArray(sampledValues, simplify=TRUE))
}

prion.time.steps <- function(sampleTimes, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass) {

	sampleTimes <- sort(sampleTimes)
	polymerDistribution0_innerhm <- convert.integer.vectors.hm(polymerLengths0_size, polymerLengths0_count)
	sampledValues <- computationClass$prionGillespieInfectiousNonInfectious(.jarray(as.numeric(sampleTimes)), as.numeric(lambda), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))

	return(.jevalArray(sampledValues, simplify=TRUE))
}

prion.polymer.distribution <- function(endTime, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0, computationClass) {

	polymerDistribution0_innerhm <- convert.integer.vectors.hm(polymerLengths0_size, polymerLengths0_count)
	polymerDistributionFinal_innerhm <- computationClass$prionGillespiePolymerLengths(as.numeric(endTime), as.numeric(lambda), as.numeric(beta), as.numeric(delta_p), as.numeric(b), as.integer(polymerThreshold), polymerDistribution0_innerhm, as.integer(monomers0))

	return(convert.integer.hm.vectors(polymerDistributionFinal_innerhm))
}

convert.integer.vectors.hm <- function(keys, values, computationClass) {
	if(length(keys) != length(values))
		stop("keys and values must have same length")

	innerhm <- computationClass$createEmptyInnerHashMap()
	for(i in 1:length(keys)){
		innerhm$put(as.integer(keys[i]), as.integer(values[i]))
	}

	return(innerhm)
}

convert.integer.hm.vectors <- function(innerhm) {
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

library("rJava")
.jinit(".")
computationClass <- J("ComputationUtils")

keys <- c(1, 2, 3, 4)
values <- c(10, 20, 30, 40)
innerhm <- convert.integer.vectors.hm(keys, values, computationClass)
returnedSet <- convert.integer.hm.vectors(innerhm)
print(returnedSet)
