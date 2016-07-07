
CC = g++
LINK_FLAGS = -lm
CCFLAGS = -O3 -fopenmp

all:
	$(CC) $(CCFLAGS) main.c sparsehash.c MurmurHash3.cpp -o main $(LINK_FLAGS)

