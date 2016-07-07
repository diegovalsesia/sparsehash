
#include "sparsehash.h"


static inline uint64_t rand_64 (){
	return (((uint64_t) rand() <<  0) & 0x00000000FFFFFFFFull) | (((uint64_t) rand() << 32) & 0xFFFFFFFF00000000ull);
}


// Data is array of strings
static void sparsehash_compute_char(char **data, uint32_t num_elements, uint16_t *str_len, double gamma, uint16_t m, char *out){

	uint32_t h, i, ibyte;
	uint32_t *seeds;
	uint64_t hash[2];
	double tau;


	seeds = (uint32_t*) malloc(sizeof(uint32_t)*m);
	for (i = 0; i < m; i++)	{
		seeds[i] = rand();
	}
	
	tau = gamma*UINT64_MAX;

	#pragma omp parallel for private(i,ibyte,h,hash) shared(data, str_len)
	for (i = 0; i < m; i++) { 

		for ( h=0; h<num_elements; h++) {

			MurmurHash3_x64_128( data[h], str_len[h], seeds[i], hash );
			if (hash[0] < tau) {
				ibyte = i/8;
				out[ibyte] =  out[ibyte] | ( (0x80) >> (i%8) );
				break;
			}

		}

	}

	free(seeds);

}


static void sparsehash_compute_char_fast(char **data, uint32_t num_elements, uint16_t *str_len, double gamma, uint16_t m, char *out){

	uint32_t h, i, ibyte;
	uint64_t tau;
	uint64_t *bot, *top;
	uint64_t hash[2];
	uint64_t *hashes;


	tau = (uint64_t)(gamma*UINT64_MAX);

	bot = (uint64_t*) malloc(sizeof(uint64_t)*m);
	top = (uint64_t*) malloc(sizeof(uint64_t)*m);
	for (i = 0; i < m; i++)	{
		bot[i] = rand_64();
		if (bot[i]>UINT64_MAX-tau)
			top[i] = tau-(UINT64_MAX-bot[i]);
		else
			top[i] = bot[i] + tau;
	}
	
	hashes = (uint64_t*)malloc(sizeof(uint64_t)*num_elements);

	#pragma omp parallel for private(h,hash)
	for ( h=0; h<num_elements; ++h) {
		MurmurHash3_x64_128 ( data[h], str_len[h], bot[0], hash );
		hashes[h] = hash[0];
	}

	#pragma omp parallel for private(i,ibyte,h)
	for (i = 0; i < m; ++i) { 

		for ( h=0; h<num_elements; ++h) {

			if ((hashes[h] < top[i]) && (hashes[h] >= bot[i])) {
				ibyte = i/8;
				out[ibyte] =  out[ibyte] | ( (0x80) >> (i%8) );
				break;
			}

		}

	}

	free(bot);
	free(top);
	free(hashes);

}




// Data is array of 16-bit integers
static void sparsehash_compute_uint16(uint16_t *data, uint32_t num_elements, double gamma, uint16_t m, char *out){

	uint32_t h, i, ibyte;
	uint32_t *seeds;
	uint64_t hash[2];
	double tau;


	seeds = (uint32_t*) malloc(sizeof(uint32_t)*m);
	for (i = 0; i < m; i++)	{
		seeds[i] = rand();
	}
	
	tau = gamma*UINT64_MAX;

	#pragma omp parallel for private(i,ibyte,h,hash) shared(data)
	for (i = 0; i < m; i++) { 

		ibyte = i/8;

		for ( h=0; h<num_elements; h++) {

			MurmurHash3_x64_128( &(data[h]), 2, seeds[i], hash );
			if (hash[0] < tau) {
				ibyte = i/8;
				out[ibyte] =  out[ibyte] | ( (0x80) >> (i%8) );
				break;
			}
		}

	}

	
	free(seeds);

}

static void sparsehash_compute_uint16_fast(uint16_t *data, uint32_t num_elements, double gamma, uint16_t m, char *out){

	uint32_t h, i, ibyte;
	uint64_t tau;
	uint64_t *bot, *top;
	uint64_t hash[2];
	uint64_t *hashes;


	tau = (uint64_t)(gamma*UINT64_MAX);

	bot = (uint64_t*) malloc(sizeof(uint64_t)*m);
	top = (uint64_t*) malloc(sizeof(uint64_t)*m);
	for (i = 0; i < m; i++)	{
		bot[i] = rand_64();
		if (bot[i]>UINT64_MAX-tau)
			top[i] = tau-(UINT64_MAX-bot[i]);
		else
			top[i] = bot[i] + tau;
	}
	
	hashes = (uint64_t*)malloc(sizeof(uint64_t)*num_elements);

	#pragma omp parallel for private(h,hash)
	for ( h=0; h<num_elements; ++h) {
		MurmurHash3_x64_128 ( &(data[h]), 2, bot[0], hash );
		hashes[h] = hash[0];
	}

	#pragma omp parallel for private(i,ibyte,h)
	for (i = 0; i < m; i++) { 

		for ( h=0; h<num_elements; h++) {

			if (hashes[h] < top[i] && hashes[h] >= bot[i]) {
				ibyte = i/8;
				out[ibyte] =  out[ibyte] | ( (0x80) >> (i%8) );
				break;
			}

		}

	}

	free(bot);
	free(top);
	free(hashes);

}




// Data is array of 32-bit integers
static void sparsehash_compute_uint32(uint32_t *data, uint32_t num_elements, double gamma, uint16_t m, char *out){

	uint32_t h, i, ibyte;
	uint32_t *seeds;
	uint64_t hash[2];
	double tau;


	seeds = (uint32_t*) malloc(sizeof(uint32_t)*m);
	for (i = 0; i < m; i++)	{
		seeds[i] = rand();
	}
	
	tau = gamma*UINT64_MAX;

	#pragma omp parallel for private(i,ibyte,h,hash) shared(data)
	for (i = 0; i < m; i++) { 

		for ( h=0; h<num_elements; h++) {

			MurmurHash3_x64_128( &(data[h]), 4, seeds[i], hash );
			if (hash[0] < tau) {
				ibyte = i/8;
				out[ibyte] =  out[ibyte] | ( (0x80) >> (i%8) );
				break;
			}

		}

	}

	free(seeds);

}


static void sparsehash_compute_uint32_fast(uint32_t *data, uint32_t num_elements, double gamma, uint16_t m, char *out){

	uint32_t h, i, ibyte, hash_tmp;
	uint64_t tau;
	uint64_t *bot, *top;
	uint64_t hash[2];
	uint64_t *hashes;


	tau = (uint64_t)(gamma*UINT64_MAX);

	bot = (uint64_t*) malloc(sizeof(uint64_t)*m);
	top = (uint64_t*) malloc(sizeof(uint64_t)*m);
	for (i = 0; i < m; i++)	{
		bot[i] = rand_64();
		if (bot[i]>UINT64_MAX-tau)
			bot[i] = UINT64_MAX-tau;
		else
			top[i] = bot[i] + tau;
	}
	
	hashes = (uint64_t*)malloc(sizeof(uint64_t)*num_elements);

	#pragma omp parallel for private(h,hash)
	for ( h=0; h<num_elements; ++h) {
		MurmurHash3_x64_128 ( &(data[h]), 4, bot[0], hash );
		hashes[h] = hash[0];
	}


	#pragma omp parallel for private(i,ibyte,h)
	for (i = 0; i < m; i++) { 

		ibyte = i/8;

		for ( h=0; h<num_elements; h++) {

			if (hashes[h] < top[i] && hashes[h] >= bot[i]) {
				out[ibyte] =  out[ibyte] | ( (0x80) >> (i%8) );
				break;
			}

		}

	}

	free(bot);
	free(top);
	free(hashes);

}






void sparsehash_sketch(void *data, uint32_t num_elements, uint16_t element_size, uint16_t *str_len, uint32_t seed, double gamma, uint16_t m, char *out){

	uint16_t mbytes;


	srand(seed);

	mbytes = m/8;
	if (m%8!=0)
		mbytes++; 

	memset(out,0,mbytes);

	switch (element_size){

		case 1 : sparsehash_compute_char((char**)data, num_elements, str_len, gamma, m, out); break;
		case 2 : sparsehash_compute_uint16((uint16_t*)data, num_elements, gamma, m, out); break;
		case 4 : sparsehash_compute_uint32((uint32_t*)data, num_elements, gamma, m, out); break;

	}

}


void sparsehash_sketch_fast(void *data, uint32_t num_elements, uint16_t element_size, uint16_t *str_len, uint32_t seed, double gamma, uint16_t m, char *out){

	uint16_t mbytes;


	srand(seed);

	mbytes = m/8;
	if (m%8!=0)
		mbytes++; 

	memset(out,0,mbytes);

	switch (element_size){

		case 1 : sparsehash_compute_char_fast((char**)data, num_elements, str_len, gamma, m, out); break;
		case 2 : sparsehash_compute_uint16_fast((uint16_t*)data, num_elements, gamma, m, out); break;
		case 4 : sparsehash_compute_uint32_fast((uint32_t*)data, num_elements, gamma, m, out); break;

	}

}


double sparsehash_sim_J(const char *sketch_1, const char *sketch_2, uint16_t bit_len){

	uint32_t nzz=0, nz_1=0, nz_2=0;
	uint16_t i, byte_len, extra_bits;
	uint8_t not_temp_1, not_temp_2;
	double Jaccard;


	byte_len = bit_len/8;
	extra_bits = bit_len%8;

	if (extra_bits==0){
		not_temp_1 = ~(sketch_1[byte_len]);
		not_temp_2 = ~(sketch_2[byte_len]);
	}else{
		not_temp_1 = ~( sketch_1[byte_len] | (0xFF >> extra_bits) );
		not_temp_2 = ~( sketch_2[byte_len] | (0xFF >> extra_bits) );
	}


	nzz += __builtin_popcount( not_temp_1 & not_temp_2 );
	nz_1 += __builtin_popcount(not_temp_1);
	nz_2 += __builtin_popcount(not_temp_2);


	for (i=byte_len-1; i>=0; i--){
		
		not_temp_1 = ~(sketch_1[i]);
		not_temp_2 = ~(sketch_2[i]);

		//nzz += __builtin_popcount( ~(sketch_1[i] | sketch_2[i]) );
		nzz += __builtin_popcount( not_temp_1 & not_temp_2 );
		nz_1 += __builtin_popcount(not_temp_1);
		nz_2 += __builtin_popcount(not_temp_2);

	}

	Jaccard = log( ((double)(nz_1)*nz_2)/((double)(nzz)*bit_len) ) / log( (double)(nzz)/bit_len );

	return Jaccard;

}


uint32_t sparsehash_dist_H(const char *sketch_1, const char *sketch_2, uint16_t bit_len){

	uint32_t hamming=0;
	uint16_t i, byte_len, extra_bits;
	uint8_t temp;


	byte_len = bit_len/8;
	extra_bits = bit_len%8;

	if (extra_bits==0){
		temp = (sketch_1[byte_len]^sketch_2[byte_len]);
		hamming += __builtin_popcount(temp);
	}
	else{
		temp = (sketch_1[byte_len] | (0xFF >> extra_bits))^(sketch_2[byte_len] | (0xFF >> extra_bits));
		hamming += __builtin_popcount(temp);
	}


	for (i=byte_len-1;i>=0;i--){

		temp = (sketch_1[i]^sketch_2[i]);

		hamming += __builtin_popcount(temp);


	}

	return hamming;
	
}


double get_gamma(uint32_t sparsity){

	return ( 1 - (pow(2,-(1.0/(double)sparsity))) );

}