import mcmc


def error(actual, expected, word):
	raise Exception("expected {}, got {}, for word:\n{}".format(expected, actual, word))


for word, expected in (
	("111011101100000110110000", 8),
	("1110001100", 5),
	("1100111000", 5),
	("111000111000", 6)):

	actual = mcmc.ladder_distance(list(map(int, word)))
	if actual != expected:
		error(expected, actual, word)


for word, expected in (
	("1010", 2),
	("1100", 2),
	("11010010", 2.33),
	("10110100", 2.33),
	("11010100", 2),
	("1100110100", 2.7),
	("110100110100", 2.8)):
		actual = mcmc.average_ladder_distance(list(map(int, word)))
		if round(actual, 2) != expected:
			error(expected, actual, word)

print("Success!")
