#include "stdint.h"
#include <stdlib.h>
#include "string.h"
#include "immintrin.h"
#include "emmintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>

struct bucket *bucket;


extern struct bucket *bucket;

/*starts the first bucket*/
void *init(ssize_t capacity );


/* fills the bucket */
void *fill(struct bucket *bucket, __m256i data);

/*emptys the bucket*/
void *empty(struct bucket *bucket);
