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


def average_ladder_distance(word):
	ones = [i for i, char in enumerate(word) if char == 1]
	distances = []
	for i, start in enumerate(ones):
		for end in ones[i+1:]:
			curr_depth = 0
			min_depth = 0
			for char in word[start+1:end+1]:
				if char == 0:
					curr_depth -= 1
					min_depth = min(curr_depth, min_depth)
				else:
					curr_depth += 1
			# distance from start to root of path, plus distance from root to end
			distance = abs(min_depth) + abs(min_depth - curr_depth)
			if min_depth == 0:
				distance += 1
			distances.append(distance)
	return sum(distances) / (1.0 * len(distances))


for word, expected in (
	("1010", 2),
	("1100", 2),
	("11010010", 2.33),
	("10110100", 2.33),
	("11010100", 2),
	("1100110100", 2.7),
	("110100110100", 2.8)):
		actual = average_ladder_distance(list(map(int, word)))
		if round(actual, 2) != expected:
			error(expected, actual, word)

print("Success!")