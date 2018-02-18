library(coda)

stats = c("root_degree", "num_leaves", "height")

n = 1000
mixingTime = 0
sampleInterval = 1
numSamples = 10000000

pdf('gelman.pdf')

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

		chains[[run]] = as.mcmc(scan(actual_values_filename, nmax=1000000))
	}

	chains_list = do.call(mcmc.list, chains)
	print(gelman.diag(chains_list, transform=TRUE))
	print(gelman.plot(chains_list, transform=TRUE))
}
