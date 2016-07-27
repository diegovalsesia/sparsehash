
CC = g++
LINK_FLAGS = -lm
CCFLAGS = -O3 -fopenmp

all:
	$(CC) $(CCFLAGS) -DMULTITHREAD main.c sparsehash.c MurmurHash3.cpp utils.c -o main $(LINK_FLAGS)
