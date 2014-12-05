## model parameters
## beta - vector of size M containing redispersion rates for each age levels
## a - vector of size M containing age levels
## alpha - rate of infection in attack zone
## c - participation rate
## p - prevalence of infection
## mu - service rate
single.server <- function(beta, a = 0:(length(beta)-1), alpha, c, p, mu, javaClass) {
    treatedStats <- .jcall(javaClass,"[D","singleServer",.jarray(as.numeric(a)),.jarray(as.numeric(beta)), as.numeric(alpha), as.numeric(c), as.numeric(p), as.numeric(mu))
    return(treatedStats)
}

prion.multiple.iterations <- function() {

}

prion.time.steps <- function() {

}

prion.polymer.distribution <- function(endTime, lambda, delta_m, beta, delta_p, b, polymerThreshold, polymerLengths0_size, polymerLengths0_count, monomers0) {

	
}

convert.integer.vectors.hm <- function(keys, values) {
	if(length(keys) != length(values))
		stop("keys and values must have same length")

	hm <- .jnew("java/util/HashMap<Integer, Integer>")
	for(i in 1:length(keys)){
		.jcall(hm, "put", as.integer(keys[i]), as.integer(values[i]))
	}

	return(hm)
}

convert.integer.hm.vectors <- function(hm) {
	entrySet <- .jcall(hm, "entrySet")
	iterator <-.jcall(entrySet,"iterator")
	size <- .jcall(entrySet, "size")
	keys <- rep(0, size)
	values <- rep(0, size)
	for(index in 1:size) {
		entry <- .jcall(iterator, "next")
		keys[index] <- .jcall(entry, "getKey")
		values[index] <- .jcall(entry, "getValue")
	}

	return(list(keys=keys, values=values))
}

library("rJava")
.jinit(".")
computationClass <- .jnew("ComputationUtils")