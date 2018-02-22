import argparse
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('n', type=int, help='value of n to use')
	parser.add_argument('num_words', type=int, help='number of random words to generate')

	args = parser.parse_args()

	out_dir = 'start_words'
	print('writing {} random words to {}/'.format(args.num_words, out_dir))
	for run in range(1, args.num_words+1):
		rand_word = DyckWords(args.n).random_element()
		out_name = 'random{}_n={}.txt'.format(run, args.n)
		out_path = os.path.join(out_dir, out_name)
		with open(out_path, 'w') as f:
			f.write('{}\n'.format(''.join(map(str, rand_word))))
