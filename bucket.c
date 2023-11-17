#include <stdlib.h>
#include "stdint.h"
#include "string.h"
#include "immintrin.h"
#include "emmintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>

/* we will fill 8 buckets  */
struct bucket
{

 /*An avx vector */
 __m256i arr;
 /*integer size of the heap*/
 ssize_t size;
 /*alloc memory size*/
 ssize_t alloc;
 /*capacity*/
 ssize_t capacity;


};


typedef struct bucket bucket;



void *init(ssize_t capacity )
{

     bucket* arthox = (bucket*)malloc(sizeof(bucket));
     arthox -> alloc = sizeof(bucket);
     /**/
     arthox -> size = 0;
     arthox -> capacity = capacity;


}

void *fill(struct bucket *bucket, __m256i* data)
{

  /* the init will initialize values */
 memcpy(&bucket -> arr, data, sizeof(&data));
 bucket -> size++;


}

void *empty(struct bucket *bucket)
{

  /*This is dangerous but I like dangerous*/
  free(bucket);

}
