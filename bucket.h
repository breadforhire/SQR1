#include <stdlib.h>
#include "immintrin.h"
#include "emmintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>




/*starts the first bucket*/
extern struct table* init(ssize_t capacity );


/* fills the bucket */
extern void *fill(struct table *table, __m256i data, int key);


/*extract*/
extern __m256i extract(struct table *table, int key);

/*free memory*/
extern void *release(struct table* table);

/*print buckets*/
extern void b_p(struct table* table);
