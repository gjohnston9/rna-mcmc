### starting points chosen uniformly at random
# for run in {1..10}; do
#   echo run: $run
#   pypy mcmc.py 1000 0 1 1000000 --uniform --prefix="run${run}_" --start_word_source="start_words/random${run}_n=1000.txt";
# done

### starting at each of the four vertices of bounding polyhedron
for run in {1..4}; do
  echo run: $run
  pypy mcmc.py 1000 0 1 10000000 --uniform --prefix="run${run}_" --start_word_source="start_words/v${run}.txt";
done

alert;
