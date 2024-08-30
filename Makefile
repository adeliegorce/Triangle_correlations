# LDFLAGS= -lfftw3f_omp -lfftw3f -lm -L/usr/local/lib
# CPPFLAGS=-I/usr/local/include/
# CC=g++ -std=c++11 -fopenmp
CC = /opt/homebrew/opt/llvm/bin/clang++ -fopenmp -v
LDFLAGS= -lfftw3f_omp -lfftw3f -lm -L /opt/homebrew/lib #-L/opt/homebrew/opt/libomp/lib 
CPPFLAGS= -I /opt/homebrew/include # -I/opt/homebrew/opt/libomp/include 

2d: Spherical_correlations2d.cpp 
	$(CC) $(CPPFLAGS) -o ./SC_2d.o ./Spherical_correlations2d.cpp  $(LDFLAGS) 

3d: Spherical_correlations3d.cpp 
	$(CC) $(CPPFLAGS) -o ./SC_3d.o ./Spherical_correlations3d.cpp  $(LDFLAGS)

2d_cross: Spherical_correlations2d.cpp 
	$(CC) $(CPPFLAGS) -o ./CC_2d.o ./Cross_correlations2d.cpp  $(LDFLAGS) 

clean:
	rm -rf *.o
	rm -rf *~
