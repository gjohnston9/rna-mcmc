library(dgof)

stats = c("root_degree", "num_leaves", "height")
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
		statistic_vals = c()
		p_vals = c()
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
			# cat(sprintf("length:%d\n", length(actual_values)))
			actual_values = actual_values[seq(1, length(actual_values), interval)]
			# cat(sprintf("sampleSize:%d\nnew length:%d\n", sampleSize, length(actual_values)))
			results = dgof::ks.test(actual_values, f)

			# print(results)
			# browser()

			statistic_vals = append(statistic_vals, results["statistic"][[1]][[1]])
			p_value = results["p.value"][[1]][[1]]
			# print(p_value)
			p_vals = append(p_vals, p_value)
		}

		# cat(sprintf(
		# 	"statistic:\n%f\n\np-value:\n%f\n\n",
		# 	mean(statistic_vals),
		# 	mean(p_vals)))

		# cat(sprintf(
		# 	"p-value:\n%f\n",
		# 	mean(p_vals)))

		cat(sprintf(
			"statistic:\n\t%.3f\n",
			mean(statistic_vals)))
	}
}