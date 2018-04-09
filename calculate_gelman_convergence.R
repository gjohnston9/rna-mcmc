library(coda)

stats = c("root_degree", "num_leaves", "height")
limits = c(1000, 10000, 100000, 500000, 1000000, 10000000)

### If "plot_at" matches a value in "limits", a plot of number of iterations vs. shrink factor will be produced,
### 	with number of iterations going from 0 to "plot_at".
### Set to -1 or any number not in "limits" to disable plotting.
### For a value of 10,000,000, plotting takes 20+ minutes to finish
plot_at = -1

### If not plotting anything, make sure to comment out this line.
### Otherwise the plot will be overwritten with an empty pdf.
# pdf('gelman.pdf')

n = 1000
mixingTime = 0
sampleInterval = 1
numSamples = 10000000

num_runs = 4

stopifnot(all(limits <= numSamples))

for (limit in limits) {
	for (stat in stats) {
		cat(sprintf(
			"calculating gelman convergence for %s with n=%d, using %d samples\n",
			stat,
			n,
			limit))

		chains = vector("list", num_runs)
		for (run in 1:num_runs) {
			actual_values_filename = sprintf(
				"data/by_sample/run%d_%s_n=%d_dist=nntm_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
				run,
				stat,
				n,
				mixingTime,
				sampleInterval,
				numSamples)

			chains[[run]] = as.mcmc(scan(actual_values_filename, nmax=limit, quiet=TRUE))
		}

		chains_list = do.call(mcmc.list, chains)
		print(gelman.diag(chains_list, transform=TRUE))

		if (limit == plot_at) {
			cat(sprintf("creating plot for limit %d\n", limit))
			my_plot = gelman.plot(chains_list, transform=TRUE)
			title(stat)
			print(my_plot)
		}
	}
}
