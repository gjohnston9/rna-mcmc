for i in -0.9,-1.8,-1.7,-8.8 -0.9,-1.2,-1.7,-10.8 2.3,1.3,-0.4,-14.6 2.2,1.9,-0.4,-16.2 -2.8,-3.0,0.9,-10.3 -2.8,-2.2,0.9,-12.1;
do
    IFS=',' read c1 c2 c3 c4 <<< ${i};
    pypy bin/main/mcmc.py 1000 20000000 10000 10000 --nntm --c1 $c1 --c2 $c2 --c3 $c3 --c4 $c4
done
