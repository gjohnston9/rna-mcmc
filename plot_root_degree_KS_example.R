library(dgof)

stats = c("root_degree", "num_leaves", "height")
limits = c(1000, 10000, 100000, 500000, 1000000)

root_degree = c(31, 47, 431, 1886, 3624)
num_leaves = c(2, 9, 97, 478, 958)
height = c(2, 12, 99, 529, 1062)
sampleSizes = data.frame(root_degree, num_leaves, height)

plot_at_stat = "root_degree"
plot_at_limit = 100000

num_runs = 4
n = 1000

mixingTime = 0
sampleInterval = 1
numSamples = 1000000

pdf("results/root_degree_KS_example.pdf")

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
			actual_values = actual_values[seq(1, length(actual_values), interval)]
			results = dgof::ks.test(actual_values, f)

			if ((run == 1) && (stat == plot_at_stat) && (limit == plot_at_limit)) {
				my_ecdf = ecdf(actual_values)
				plot(f, xval=0:10, xlim=c(0,10), col="blue", main="experimental and actual CDFs for root degree (n=1000)")
				plot(my_ecdf, verticals=T, add=T, col="red")
				x = 2.8
				y0 = min(f(x), my_ecdf(x))
				y1 = max(f(x), my_ecdf(x))
				arrows(x, y0, x, y1, length=0.05, code=3)
				text(x + 0.05, (y0+y1)/2, labels=sprintf("value of KS statistic: %.3f", y1-y0), pos=4, cex=1)
				legend(6, 0.3, c("experimental", "actual"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))
				quit()
			}
		}
	}
}