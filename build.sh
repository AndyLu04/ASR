gcc -c 2d_array_allocation.c
gcc -c ASR.c
gcc -c datafiltering.c
gcc -c main.c
gcc -c sqrtm.c
gcc -o ASR 2d_array_allocation.o ASR.o datafiltering.o main.o sqrtm.o  -s -llapack -lm -lblas  /usr/local/lib/libblas.a /usr/local/lib/libfftw3.a /usr/local/lib/liblapack.a
