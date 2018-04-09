library(dgof)

n = 1000

stats = c("cd_averages", "avg_branching", "cd_sums", "root_degree", "num_leaves", "height")
samples = c("2-", "2-", "1-", "1-", "1-", "1-") # distinguishes between 1- and 2-sample KS test (1-sample is discrete, 2-sample is continuous)
dirs = c("by_sample", "by_sample", "by_frequency", "by_sample", "by_sample", "by_sample")
n_vals = c(2*n-1, n, 2*n-1, n, n, n)

mixingTime = 10000000
sampleInterval = 10000
numSamples = 10000

### parameters for uniform data to load for 2-sample test
uniform_mixingTime = 100000
uniform_sampleInterval = 1000
uniform_numSamples = 10000

for (j in 1:length(stats)) {
	stat = stats[j]
	n_val = n_vals[j]
	dir = dirs[j]
	sample = samples[j]

	if (sample == "1-") {
		print("performing 1-sample discrete KS test")
		expected_values_filename = sprintf(
			"expectations/%s_n=%d_cdf.txt",
			stat,
			n)

		expected_values = scan(expected_values_filename, quiet=TRUE)
		x = 1:n_val
		y = c(0, expected_values)
		f = stepfun(x, y)
	} else if (sample == "2-") {
		print("performing 2-sample continuous KS test")
		uniform_values_filename = sprintf(
			"data/%s/%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
			dir,
			stat,
			n,
			uniform_mixingTime,
			uniform_sampleInterval,
			uniform_numSamples)
		uniform_values = scan(uniform_values_filename, quiet=TRUE)
	} else {
		stop("unrecognized KS sample value")
	}
	
	print(sprintf("calculating KS for %s with n=%d, using %d samples", stat, n, numSamples))
	ks_statistic_vals = c()
	ks_p_vals = c()
	
	actual_values_filename = sprintf(
		"data/%s/%s_n=%d_dist=nntm_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
		dir,
		stat,
		n,
		mixingTime,
		sampleInterval,
		numSamples)

	actual_values = scan(actual_values_filename, quiet=TRUE)

	# cat(actual_values[1:5])
	# cat("\n")
	# cat(uniform_values[1:5])
	# cat("\n")

	if (sample == "1-") {
		ks_results = dgof::ks.test(actual_values, f)
	} else {
		ks_results = stats::ks.test(actual_values, uniform_values)
	}

	ks_statistic_vals = append(ks_statistic_vals, ks_results["statistic"][[1]][[1]])
	ks_p_value = ks_results["p.value"][[1]][[1]]
	ks_p_vals = append(ks_p_vals, ks_p_value)

	cat(sprintf(
	 	"KS statistic:\n%f\n\nKS p-value:\n%f\n\n",
	 	mean(ks_statistic_vals),
	 	mean(ks_p_vals)))
	
}

warnings()
