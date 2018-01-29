library(coda)

stats = c("num_leaves", "height", "root_degree")
for (stat in stats)
{
	filename = sprintf("data/by_sample/%s_n=1000_dist=uniform_mixingTime=10000000_sampleInterval=1_numSamples=1000000.txt", stat)
	t = read.table(filename)
	colnames(t)[1] <- stat
	print(effectiveSize(t))
}
