#include <stdlib.h>
#include "stdint.h"
#include "string.h"
#include "immintrin.h"
#include "emmintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>

/* we will fill 8 buckets  */





typedef struct table
{


 /*size of the table*/
 ssize_t size;

 ssize_t capacity;

 __m256i bucket[8];

 /*allocation for memory size*/
 ssize_t alloc;

} table;





void *init(ssize_t capacity)
{

     table* arthox = (table*)malloc(sizeof(table));
     arthox -> alloc = sizeof(table);
     /**/
     arthox -> size = 0;
     arthox -> size = capacity;


}

void fill(struct table* table, __m256i data, int key)
{

  /* the init will initialize values */


 //&table.bucket[key] =  _mm256_load_si256((__m256i*) data);
 //table -> size++;




}


void *empty()
{


}
