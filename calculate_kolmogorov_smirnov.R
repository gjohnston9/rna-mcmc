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
		# cat(sprintf("total: %f\n", total))
	}
	# cat(sprintf("coef: %.20f\n", coef))
	# cat(sprintf("coef == 0: %s\n", coef == 0))
	return(coef*total)
}


### Pr(tree of size on has height <= h)
height_cdf <- function(n, h) {
	total = 0
	for (k in 1:(((n+1)/h) + 2)) {
		numerators_1 = n:(n + 2 - k - k*h)
		denominators_1 = (n - 1 + k + k*h):(n + 1)

		numerators_2 = n:(n - k*h)
		denominators_2 = (n + k*h + 1):(n + 1)

		numerators_3 = n:(n - k - k*h)
		denominators_3 = (n + 1 + k + k*h):(n+1)

		total = total + prod(numerators_1 / denominators_1)
		total = total + 2*prod(numerators_2 / denominators_2)
		total = total + prod(numerators_3 / denominators_3)
	}

	total = (n + 1) * total
	return(total)
}


testing = F
if (!testing) {
	stats = c("root_degree", "num_leaves", "height")

	n = 1000
	mixingTime = 10000000
	sampleInterval = 1075
	numSamples = 930

	for (stat in stats) {
		print(sprintf("calculating KS for %s with n=%d, using %d samples", stat, n, numSamples))

		actual_values_filename = sprintf(
			"data/by_sample/%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
			stat,
			n,
			mixingTime,
			sampleInterval,
			numSamples)

		expected_values_filename = sprintf(
			"expectations/%s_n=%d_cdf.txt",
			stat,
			n)

		actual_values = scan(actual_values_filename)
		expected_values = scan(expected_values_filename)

		# expected_values_function <- function(indices) {cat(sprintf("%s\n", indices))}
		expected_values_function <- function(indices) { return(expected_values[indices]) }
		print(dgof::ks.test(actual_values, "expected_values_function"))
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


### Calculates number of trees of size n of size >= h
### May overflow
min_height <- function(n, h) {
    total = 0
    for (k in 1:((ceiling(n+1)/h) + 1)) {
        total = total + choose(2*n, n + 1 - k*h) - 2*choose(2*n, n - k*h) + choose(2*n, n - 1 - k*h)
    }
    return(total)
}


### Pr(tree of size on has height <= h), with no
### attempt to avoid overflow
height_cdf_naive <- function(n, h) {
	return(1 - (min_height(n, h) / Cn(n)))
}


if (testing) {
	n = 5
	r = 5
	cat(sprintf("root deg prob: %f\n", root_degree_cdf(n, r)))
	cat(sprintf("root deg prob: %f (naive)\n", root_degree_cdf_naive(n, r)))
	
	cat(sprintf("num leaves prob: %f\n", num_leaves_cdf(n, r)))
	cat(sprintf("num leaves prob: %f (naive)\n", num_leaves_cdf_naive(n, r)))

	cat(sprintf("height prob: %f\n", height_cdf(n, r)))
	cat(sprintf("height prob: %f (naive)\n", height_cdf_naive(n, r)))	

	cat(sprintf("# of trees of size %d with at least %d leaves: %d\n", n, r, min_height(n, r)))
	cat(sprintf("Catalan number for n=%d: %d\n", n, Cn(n)))

	cat("\nwarnings:\n")
	warnings()
}
