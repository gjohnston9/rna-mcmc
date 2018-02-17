library(coda)

m1 = as.mcmc(v1)
m2 = as.mcmc(v2)

mcmc_list = mcmc.list(m1, m2)

stats = c("root_degree", "num_leaves", "height")

n = 1000
mixingTime = 0
sampleInterval = 0
numSamples = 30000000

for (stat in stats) {
	print(sprintf(
		"calculating gelman convergence for %s with n=%d, using %d samples",
		stat,
		n,
		numSamples))

	chains = vector("list", 10)
	for (run in 1:10) {
		actual_values_filename = sprintf(
			"data/by_sample/run%d_%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
			run,
			stat,
			n,
			mixingTime,
			sampleInterval,
			numSamples)

		chains[[run]] = as.mcmc(scan(actual_values_filename))
	}

	chains_list = do.call(mcmc.list, chains)
	print(gelman.diag(chains_list))
}
