library(dgof)

# stats = c("root_degree", "num_leaves", "height")
stats = c("root_degree")
limits = c(1000, 10000, 100000, 500000, 1000000)

root_degree = c(31, 47, 431, 1886, 3624)
num_leaves = c(2, 9, 97, 478, 958)
height = c(2, 12, 99, 529, 1062)
sampleSizes = data.frame(root_degree, num_leaves, height)

num_runs = 4

n = 1000
# mixingTime = 10000000
# sampleInterval = 1064
# numSamples = 940

mixingTime = 0
sampleInterval = 1
numSamples = 1000000

for (stat in stats) {
	expected_values_filename = sprintf(
		"expectations/%s_n=%d_cdf.txt",
		stat,
		n)

	expected_values = scan(expected_values_filename)
	x = 1:n
	y = c(0, expected_values)
	f = stepfun(x, y)
	
	for (i in 1:length(limits)) {
		limit = limits[i]
		sampleSize = sampleSizes[i, stat]

		print(sprintf("calculating KS for %s with n=%d, using %d samples", stat, n, limit))
		ks_statistic_vals = c()
		ks_p_vals = c()
		cvm_statistic_vals = c()
		cvm_p_vals = c()
		for (run in 1:num_runs) {
			actual_values_filename = sprintf(
				"data/by_sample/run%d_%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
				run,
				stat,
				n,
				mixingTime,
				sampleInterval,
				numSamples)

			actual_values = scan(actual_values_filename, nmax=limit, quiet=TRUE)
			interval = round(limit / sampleSize)
			### get every kth element, ending up with $sampleSize samples
			actual_values = actual_values[seq(1, length(actual_values), interval)]
			ks_results = dgof::ks.test(actual_values, f)
			cvm_results = dgof::cvm.test(actual_values, f)

			ks_statistic_vals = append(ks_statistic_vals, ks_results["statistic"][[1]][[1]])
			cvm_statistic_vals = append(cvm_statistic_vals, cvm_results["statistic"][[1]][[1]])
			ks_p_value = ks_results["p.value"][[1]][[1]]
			cvm_p_value = cvm_results["p.value"][[1]][[1]]
			ks_p_vals = append(ks_p_vals, ks_p_value)
			cvm_p_vals = append(cvm_p_vals, cvm_p_value)
		}

		cat(sprintf(
		 	"KS statistic:\n%f\n\nKS p-value:\n%f\n\n",
		 	mean(ks_statistic_vals),
		 	mean(ks_p_vals)))

		cat(sprintf(
		 	"CVM statistic:\n%f\n\nCVM p-value:\n%f\n\n",
		 	mean(cvm_statistic_vals),
		 	mean(cvm_p_vals)))
	}
}

warnings()