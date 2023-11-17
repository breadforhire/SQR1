#include "stdint.h"
#include <stdlib.h>
#include "string.h"
#include "immintrin.h"
#include "emmintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>




struct table
{


 /*size of the table*/
 ssize_t size;

 /*the capacity for the entire */
 ssize_t capacity;

 __m256i bucket[8];

 /*allocation for memory size*/
 ssize_t alloc;

} ;

struct table *table;


extern struct table *table;



/*starts the first bucket*/
void *init(ssize_t capacity );


/* fills the bucket */
void *fill(struct table *table, __m256i data, int key);


/*checks if table is empty*/
void *empty(struct table *table);
