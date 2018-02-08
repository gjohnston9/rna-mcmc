library(dgof)

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
