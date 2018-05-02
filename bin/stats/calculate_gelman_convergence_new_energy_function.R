library(coda)

stats = c("root_degree", "num_leaves", "height")
limits = c(seq(from=5e5, to=2e6, by=5e5), seq(from=3e6, to=10e6, by=1e6), seq(15e6, 40e6, by=5e6))

# 500,000
# 1,000,000
# 1,500,000
# 2,000,000

# 3,000,000
# 4,000,000
# 5,000,000
# 6,000,000
# 7,000,000
# 8,000,000
# 9,000,000
# 10,000,000

# 15,000,000
# 20,000,000
# 25,000,000
# 30,000,000
# 35,000,000
# 40,000,000

### If "plot_at" matches a value in "limits", a plot of number of iterations vs. shrink factor will be produced,
### 	with number of iterations going from 0 to "plot_at".
### Set to -1 or any number not in "limits" to disable plotting.
### For a value of 10,000,000, plotting takes 20+ minutes to finish
plot_at = 40e6

### If not plotting anything, make sure to comment out this line.
### Otherwise the plot will be overwritten with an empty pdf.
pdf('results/gelman_new_energy_parameters.pdf')

n = 1000
mixingTime = 0
sampleInterval = 1
numSamples = 40000000

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
				"data/by_sample/run%d_%s_n=%d_dist=nntm_mixingTime=%d_sampleInterval=%d_numSamples=%d_-0.9_-1.8_-1.7_-8.8.txt",
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
