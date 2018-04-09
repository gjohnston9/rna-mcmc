library(dgof)

stat = "root_degree"
sample_size = 3624
plot_at_limit = 100000

num_runs = 4
n = 1000

mixingTime = 0
sampleInterval = 1
numSamples = 1000000

pdf("results/root_degree_KS_example.pdf")

expected_values_filename = sprintf(
	"expectations/%s_n=%d_cdf.txt",
	stat,
	n)

expected_values = scan(expected_values_filename)
x = 1:n
y = c(0, expected_values)
f = stepfun(x, y)

print(sprintf("calculating KS for %s with n=%d, using %d samples", stat, n, plot_at_limit))
statistic_vals = c()
p_vals = c()
actual_values_filename = sprintf(
	"data/by_sample/run1_%s_n=%d_dist=uniform_mixingTime=%d_sampleInterval=%d_numSamples=%d.txt",
	stat,
	n,
	mixingTime,
	sampleInterval,
	numSamples)

actual_values = scan(actual_values_filename, nmax=plot_at_limit, quiet=TRUE)
interval = round(plot_at_limit / sample_size)
### get every kth element, ending up with $sample_size samples
actual_values = actual_values[seq(1, length(actual_values), interval)]
results = dgof::ks.test(actual_values, f)

my_ecdf = ecdf(actual_values)
plot(f, xval=0:10, xlim=c(0,10), col="blue", main="experimental and actual CDFs for root degree (n=1000)")
plot(my_ecdf, verticals=T, add=T, col="red")
x = 2.8
y0 = min(f(x), my_ecdf(x))
y1 = max(f(x), my_ecdf(x))
arrows(x, y0, x, y1, length=0.05, code=3)
text(x + 0.05, (y0+y1)/2, labels=sprintf("value of KS statistic: %.3f", y1-y0), pos=4, cex=1)
legend(6, 0.3, c("experimental", "actual"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))
