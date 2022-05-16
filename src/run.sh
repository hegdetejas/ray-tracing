#!/bin/bash                                                                                                                                                                                                                                                                                                                                                                                                 

g++-11 ./ray_tracing.cpp -O3 -mtune=native -march=native -ffast-math
./a.out

rm ./a.out

python3 visualize.py

rm ../out/*.txt