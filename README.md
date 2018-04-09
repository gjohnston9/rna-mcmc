# Markov-Chain Monte Carlo
## Dependencies
- Python
	- 2.6+
	- NumPy to run `mcmc.py`
	- PyPy is not required, but in my experience it has led to a speedup of 5-10x over CPython (the default implementation of Python)
	- Matplotlib and NumPy to run `create_histograms.py`
- Sage
	- 7.6+, but earlier versions should work too
- R
	- The [coda](https://cran.r-project.org/web/packages/coda/index.html) package for `calculate_effective_sample_size.R` and `calculate_gelman_convergence.R`.
	- The [dgof](https://cran.r-project.org/web/packages/dgof/index.html) package for `calculate_kolmogorov_smirnov.R` and `calculate_kolmogorov_smirnov_averages.R`
## Running the chain
- `mcmc.py` requires four positional arguments: *n*, initial mixing time, gap size (number of moves to make before recording each sample), and number of samples. It also requires either the `--uniform` or `--nntm` flag to specify which distribution to use when making each move.
- Each time 10% of the total number of samples has been collected, the program will print a message and append per-sample data for each characteristic to files in the `data/by_sample` directory. This is done to reduce memory usage of the program since after writing these values, they can be discarded.
- After the program finishes, data about frequencies of characteristics (how many times each characteristic had a certain value) is written to the `data/by_frequency` directory.
## Creating plots
- To plot data output by the chain, the first step is to run `write_plot_data.sage`. This has slightly different behavior depending on whether the `--uniform` or `--nntm` flag is used.
	- With either flag, the first command line argument should be your value of *n* (1000 in most cases).
	- With the `--uniform` flag, you should provide initial mixing time, sample interval, and number of samples as command line arguments, in that order (this order is the same for any script in the project that requires command line arguments). For each characteristic that has an available CDF (height, number of leaves, root degree, and contact distance sums), a text file with two columns (experimental and expected values) will be created in the `data/processed_plot_data` directory.
	- With the `--nntm` flag, you should first provide uniform initial mixing time, sample interval, and number of samples as command line arguments. Then you should provide these three values used for the run under the nntm distribution. This will create files to the `data/processed_plot_data` directory, comparing experimental values under the nntm distribution to expected values under the uniform distribution. In addition, it will create two files comparing experimental values under the nntm distribution for ladder distance and contact distance averages with the experimental values for these characteristics under the uniform distribution.
	- Note: The comparison of nntm experimental data to uniform experimental data when the `--nntm` flag is used is the reason that both uniform and nntm parameters are required as command line arguments in that case.
- After the appropriate data is in the `data/processed_plot_data` directory, use the `plot_results.sage` and `create_histograms.py` scripts to create plots.
	- `plot_results.sage` creates scatterplots for characteristics that are integer-valued: ladder distance, height, contact distance sums, number of leaves, and root degree.
	- `create_histograms.py` creates histograms for characteristics that are non-negative and real-valued: average branching and contact distance averages.
	- Both of these scripts save plots to the `plots` directory.
	- `plot_results.sage` requires as command line arguments *n*, initial mixing time, gap size, and number of samples, followed by one of `--uniform` or `--nntm`.
	- `create_histograms.py` requires as command line arguments *n*, uniform initial mixing time, gap size, and number of samples, and nntm initial mixing time, gap size, and number of samples.

## Calculating statistics
- Detailed descriptions of the uses of these scripts is in the **experiments** section. This section describes instructions for use of these scripts and a short description of their functions.

| Script name  | Inputs | Function | Output |
| :-----------:  | :------: | :-----: | :-----: |
| `calculate_effective_sample_size.R`  |  At the top of the file: values for *n*, initial mixing time, gap size, and number of samples, as well as which characteristics to examine, and which initial parts of the chain to calculate gelman convergence for ("limits"). | For each given characteristic, calculates the effective sample size (an estimate for the number of independent samples) for each chain. | Prints the mean and standard deviation of effective sample size, for each characteristic and number of samples. |
| `calculate_gelman_convergence.R`  |  At the top of the file: values for *n*, initial mixing time, gap size, and number of samples, as well as which characteristics to examine, which initial parts of the chain to calculate gelman convergence for ("limits"), and optionally a limit for which a plot of number of iterations vs. shrink factor will be produced. | For each given characteristic, examines first portions of each of four chains (the first 1,000 samples, the first 10,000 samples, etc.) and uses Gelman and Rubin's convergence diagnostic to calculate a shrink factor. This diagnostic compares the variance between the chains to the variance within each chain, and produces a "shrink factor." A shrink factor of 1 indicates that between variance and within variance are equal, whereas a larger shrink factor indicates that there is still a noticeable difference in the two variances.  | Prints estimate of shrink factor and its upper confidence limit for each characteristic and limit value. |
| `calculate_kolmogorov_smirnov_averages.R` | At the top of the file: values for *n*, initial mixing time, gap size, number of samples, and distribution (nntm or uniform). Also, a list of which characteristics to examine, and a list of which initial parts of the chain to calculate KS for ("limits"). | For each given characteristic, examines first portions of each of four chains (the first 1,000 samples, the first 10,000 samples, etc.) and calculates KS for that portion of each chain, comparing it to the expected values under the uniform distribution. | For each characteristic and limit, prints the mean KS statistic and mean p-value. A p-value below our significance level of 0.05 indicates that the two underlying distributions definitely differ. A p-value of 0 can be expected for the nntm distribution, and a p-value much larger than 0.05 can be expected when the uniform distribution is used for experimental values. |
| `calculate_kolmogorov_smirnov.R` | At the top of the file: values for *n*, initial mixing time, gap size, and number of samples. Also, a list of which characteristics to examine, a list indicating whether a 1-sample discrete or 2-sample continuous KS test should be performed for a characteristic, a list indicating which directory (`data/by_sample` or `data/by_frequency`) data for a characteristic can be found in, and a list indicating the maximum value that each characteristic can have. | **Note**: this script is used for a different experiment than the three above scripts, and requires data from only one uniform and one nntm MCMC run, whereas the above scripts require data from multiple runs under the uniform or nntm distribution.<br/><br/> For each given characteristic, used the Kolmogorov-Smirnov test to test whether values for characteristics obtained under the nntm distribution could have been generated by the uniform distribution. <br/><br/> Produces some warnings about ties for average branching, even though after checking the inputs, there don't seem to be any ties. | Prints KS statistic and p-value for each characteristic. A p-value below our significance level of 0.05 indicates that the two underlying distributions definitely differ. |

## Other scripts
- Detailed descriptions of the uses of these scripts is in the **experiments** section. This section describes instructions for use of these scripts and a short description of their functions.

| Script name  | Inputs | Function | Output |
| :-----------:  | :------: | :-----: | :-----: |
| `calculate_cdfs.sage`  | Via command line: value for *n*.   | For each of the characteristics of contact distance sums, number of leaves, root degree, and height, this script writes data to a file in the `expectations` directory. The *n*th line in one of these files gives the probability that under the **uniform** distribution, a sample will have a value less than or equal to *n* for the file's characteristic. | Creates text files named, for example, `num_leaves_n=1000_cdf.txt` in the `expectations` directory. |
| `generate_polyhedron_vertices.py` | At the top of the file: value for *n*. | For the given value of *n*, creates Dyck word representations of four vertices of the bounding polyhedra, depicted in Figure 3 in page 765 of Hower and Heitsch: Parametric analysis of RNA branching configurations (2011). | Creates text files `v1.txt`, `v2.txt`, `v3.txt`, and `v4.txt` in the `start_words` directory. |
| `write_random_Dyck_words.sage` | Via command line: value for *n*, and number of random words to generate. | Uses Sage's `DyckWord` `random_element()` method to generate random Dyck words of order *n*. | Creates text files named, for example, `random1_n=1000.txt` in the `start_words` directory. |
| `multiple_runs.sh` | none | Bash script to run `mcmc.py` with the uniform distribution four times, and with the nntm distribution four times. For each distribution, starts once from each of the four vertices written in `generate_polyhedron_vertices.py`. | `mcmc.py` produces output in the `data/by_frequency` and `data/by_sample` directories. All data files generated from this script start with `run1`, `run2`, `run3`, or `run4`. |

## Experiments

#### Determining  adequate parameters (gap size and initial mixing time) for uniform and nntm distributions
- `generate_polyhedron_vertices.py`
	- Generate starting points for the `multiple_runs` script.
- `multiple_runs.sh`
	- Run `mcmc.py` with the uniform distribution four times, and with the nntm distribution four times.
- `calculate_gelman_convegence.R`
	- For each "limit" (initial portion of the chain), we see how across-chain variance compares to within-chain variance, giving us an idea of how close the chains are to converging to the same distribution. A value below 1.1 tells us that that limit is an adequate initial mixing time.
- `calculate_effective_sample_size.R`
	- For a given limit (number of iterations), for each characteristic we divide number of iterations by effective sample size. Taking the maximum across characteristics tells us what gap size we can use with an initial mixing time equal to the given limit.
- (optionally) `calculate_cdfs.sage`, then `calculate_kolmogorov_smirnov_averages.R`
	- Compares the values obtained under either the uniform or nntm in `multiple_runs`to the expected values under the uniform distribution. This tells us whether we can conclude that the distributions for any of these characteristics significantly differ under the uniform vs. the nntm distributions. A p-value below 0.05 means we can conclude this, but a p-value greater than or equal to 0.05 means we cannot make any conclusion.

#### Using Kolmogorov-Smirnov to find a difference between values for characteristics under uniform and nntm distributions
- `mcmc.py` under both uniform and nntm distributions
- `calculate_cdfs.sage` so that we can compare our nntm values with expected values under the uniform distributions, for characteristics that have a CDF available.
- `calculate_kolmogorov_smirnov.R` to determine whether we can conclude that the distributions for any of these characteristics significantly differ under the uniform vs. the nntm distributions. A p-value below 0.05 means we can conclude this, but a p-value greater than or equal to 0.05 means we cannot make any conclusion.

#### Visually comparing values for characteristics under uniform and nntm distributions
- `mcmc.py` under both uniform and nntm distributions
- `calculate_cdfs.sage` so that we can compare our nntm values with expected values under the uniform distributions, for characteristics that have a CDF available.
- `plot_results.sage` and `create_histograms.py` to generate scatterplots and histograms comparing data from the uniform and nntm distributions.

