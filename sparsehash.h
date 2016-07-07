
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MurmurHash3.h"


// Compute m-bits sketch for data. data must be an array of num_elements integers of size element_size bytes (2 or 4) 
// or an array of num_elements strings of lengths str_len (element_size=1)
void sparsehash_sketch(void *data, uint32_t num_elements, uint16_t element_size, uint16_t *str_len, uint32_t seed, double gamma, uint16_t m, char *out);

// Fast version of sparsehash with window
void sparsehash_sketch_fast(void *data, uint32_t num_elements, uint16_t element_size, uint16_t *str_len, uint32_t seed, double gamma, uint16_t m, char *out);

// Compute Jaccard estimate from two sketches
double sparsehash_sim_J(const char *sketch_1, const char *sketch_2, uint16_t bit_len);

// Compute Hamming distance between two sketches
uint32_t sparsehash_dist_H(const char *sketch_1, const char *sketch_2, uint16_t bit_len);

// Compute gamma that maximizes the entropy of the sketch
double get_gamma(uint32_t sparsity);