library(dgof)


### Pr(tree of size n has root degree <= r)
root_degree_cdf <- function(n, r) {
	if (r == 0) { return(0) }
	total = 0
	for (k in 1:r) {
		numerators = (n+1):(n-k+1)
		denominators = (2*n):(2*n-k)
		part = k*prod(numerators / denominators)
		total = total + part
	}
	return(total)
}


### Pr(tree of size n has at most L leaves)
num_leaves_cdf <- function(n, l) {
	if (l == 0) { return(0) }
	
	numerators = (n+1):1
	denominators = c(n, (2*n):(n+1))
	coef = prod(numerators / denominators)
	if (l == 1) { return(coef*n) }

	total = n ### replaces k=1 in loop below
	for (k in 2:l) {
		numerators_1 = n:(n-k+1) 
		denominators_1 = k:1
		numerators_2 = n:(n-k+2)
		denominators_2 = (k-1):1
		total = total + prod(numerators_1 / denominators_1)*prod(numerators_2 / denominators_2)
		cat(sprintf("total: %f\n", total))
	}
	cat(sprintf("coef: %.20f\n", coef))
	cat(sprintf("coef == 0: %s\n", coef == 0))
	return(coef*total)
}


vectorize_function <- function(f) {
	return(function(n, vals) {sapply(vals, function(val) {f(n, val)})})
}


fixed_n_function <- function(n, f) {
	return(function(vals) {f(n, vals)})
}


testing = F
if (!testing) {
	stats = c("root_degree", "num_leaves")
	functions = c(root_degree_cdf, num_leaves_cdf)
	stats_and_functions = mapply(list, stats, functions, SIMPLIFY=F)

	n = 1000
	mixingTime = 10000000
	sampleInterval = 1075

	numSamples = 930
	for (stat_and_function in stats_and_functions) {
		stat = stat_and_function[[1]]
		func = stat_and_function[[2]]

		filename = sprintf(
			"data/by_sample/%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
			stat,
			n,
			mixingTime,
			sampleInterval,
			numSamples)

		actual = scan(filename)

		fn_vectorized = vectorize_function(func)
		fn_fixed_n = fixed_n_function(n, fn_vectorized)
		print(dgof::ks.test(actual, "fn_fixed_n"))
	}
}


### FUNCTIONS FOR TESTING ###

### Calculates nth Catalan number (will overflow for high values of n)
Cn <- function(n) {
	return( 1/(n+1) * choose(2*n, n))
}

### Calculates number of trees of size n with root degree
### <= r. May overflow for high values of n
root_degree <- function(n, r) {
	total = 0
	for (k in 1:r) {
		total = total + (k/n)*choose(2*n-1-k, n-1)
	}
	return(total)
}

### Pr(tree of size n has root degree <= r), with no
### attempt to avoid overflow
root_degree_cdf_naive <- function(n, r) {
	return(root_degree(n, r) / Cn(n))
}

### Calculates number of trees of size n with l leaves or fewer.
### May overflow for high values of n
num_leaves <- function(n, l) {
	total = 0
	for (k in 1:l) {
		total = total + (1/n)*choose(n, k)*choose(n, k-1)
	}
	return(total)
}

### Pr(tree of size n has l leaves or fewer), with no
### attempt to avoid overflow
num_leaves_cdf_naive <- function(n, l) {
	return(num_leaves(n, l) / Cn(n))
}

if (testing) {
	n = 1000
	r = 117
	cat(sprintf("root deg prob: %f\n", root_degree_cdf(n, r)))
	cat(sprintf("root deg prob: %f (naive)\n", root_degree_cdf_naive(n, r)))
	
	cat(sprintf("num leaves prob: %f\n", num_leaves_cdf(n, r)))
	cat(sprintf("num leaves prob: %f (naive)\n", num_leaves_cdf_naive(n, n)))

	cat("\nwarnings:\n")
	warnings()
}
