for run in {1..4};
do
    echo "run: $run"
    pypy bin/main/mcmc.py 1000 0 1 40000000 --nntm --c1 -0.9 --c2 -1.8 --c3 -1.7 --c4 -8.8 --prefix="run${run}_" --start_word_source="start_words/v${run}.txt"
done
