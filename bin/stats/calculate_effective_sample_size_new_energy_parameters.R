library(coda)

stats = c("root_degree", "num_leaves", "height", "ladder_distance", "cd_averages", "avg_branching")
limits = seq(from=5e6, to=40e6, by=5e6)

n = 1000
mixingTime = 0
sampleInterval = 1
numSamples = 40e6

num_runs = 4

stopifnot(all(limits <= numSamples))

for (limit in limits) {
	for (stat in stats)	{
		cat(sprintf("\ncalculating effective sample size for %s with %d samples\n", stat, limit))
		vals = c()
		for (run in 1:num_runs)
		{
			filename = sprintf(
				"data/by_sample/run%d_%s_n=%d_dist=nntm_mixingTime=%d_sampleInterval=%d_numSamples=%d_-0.9_-1.8_-1.7_-8.8.txt",
				run,
				stat,
				n,
				mixingTime,
				sampleInterval,
				numSamples)

			my_data = scan(filename, nmax=limit, quiet=TRUE)
			stopifnot(length(my_data) != limit)

			my_ess = effectiveSize(my_data) 
			vals = append(vals, my_ess)
		}

		cat(sprintf("mean: %.0f\nstandard deviation: %.0f\n", mean(vals), sd(vals)))
	}
}
