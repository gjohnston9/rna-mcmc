import os


def write_word(list_word, vertex_number):
	filename = 'v{}.txt'.format(vertex_number)
	out_path = os.path.join(out_dir, filename)
	with open(out_path, 'w') as f:
		f.write('{}\n'.format(''.join(map(str, list_word))))

if __name__ == '__main__':
	out_dir = 'start_words'
	n = 1000

	w1 = [1]*n + [0]*n
	w2 = [1, 0]*n
	w3 = [1] + [1,0]*(n-1) + [0]
	w4 = [1,1,0]*(n/2) + [0]*(n/2)

	write_word(w1, 1)
	write_word(w2, 2)
	write_word(w3, 3)
	write_word(w4, 4)
