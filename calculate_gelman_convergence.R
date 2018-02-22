library(coda)

stats = c("root_degree", "num_leaves", "height")
limits = c(1000, 10000, 100000, 500000, 1000000)
plot_at = 10000
transform = FALSE

n = 1000
mixingTime = 0
sampleInterval = 1
numSamples = 1000000

num_runs = 4

pdf('gelman.pdf')

for (limit in limits) {
	for (stat in stats) {
		print(sprintf(
			"calculating gelman convergence for %s with n=%d, using %d samples",
			stat,
			n,
			limit))

		chains = vector("list", num_runs)
		for (run in 1:num_runs) {
			actual_values_filename = sprintf(
				"data/by_sample/run%d_%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
				run,
				stat,
				n,
				mixingTime,
				sampleInterval,
				numSamples)

			chains[[run]] = as.mcmc(scan(actual_values_filename, nmax=limit, quiet=TRUE))
		}

		chains_list = do.call(mcmc.list, chains)
		print(gelman.diag(chains_list, transform=transform))

		# if (limit == plot_at) {
		# 	my_plot = gelman.plot(chains_list, transform=transform)
		# 	title(stat)
		# 	print(my_plot)
		# }
	}
}