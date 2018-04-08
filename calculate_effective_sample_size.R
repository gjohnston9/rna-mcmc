library(coda)

# stats = c("root_degree", "num_leaves", "height")
stats = c("root_degree")

limits = c(1000, 10000, 100000, 500000, 1000000)
num_runs = 4

n = 1000
mixingTime = 0
sampleInterval = 1
numSamples = 1000000

stopifnot(all(limits <= numSamples))

for (stat in stats)
{
	for (limit in limits)
	{
		cat(sprintf("\ncalculating effective sample size for %s with %d samples\n", stat, limit))
		vals = c()
		for (run in 1:num_runs)
		{
			filename = sprintf(
				"data/by_sample/run%d_%s_n=%d_dist=nntm_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
				run,
				stat,
				n,
				mixingTime,
				sampleInterval,
				numSamples)

			my_data = scan(filename, nmax=limit)
			my_ess = effectiveSize(my_data) 
			vals = append(vals, my_ess)
		}

		cat(sprintf("mean=%.0f\nstd_dev=%.0f\n", mean(vals), sd(vals)))
	}
}
