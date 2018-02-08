library(coda)

stats = c("num_leaves", "height", "root_degree")
for (stat in stats)
{
	n = 1000
	mixingTime = 10000000

	### first run
	# sampleInterval = 1
	# numSamples = 1000000

	### second run based on first effective sample size
	# sampleInterval = 354
	# numSamples = 2820

	### let's try that again
	# sampleInterval = 1000
	# numSamples = 1000

	### let's try that again, with more accurate numbers
	sampleInterval = 1075
	numSamples = 930

	filename = sprintf(
		"data/by_sample/%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
		stat,
		n,
		mixingTime,
		sampleInterval,
		numSamples)

	t = read.table(filename)
	colnames(t)[1] <- stat
	print(effectiveSize(t))
}
