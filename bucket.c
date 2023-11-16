#include <stdlib.h>
#include "stdint.h"
#include "string.h"



/* we will fill 8 buckets  */
struct bucket
{

 /*An avx vector */
 void * arr;
 /*integer size of the heap*/
 ssize_t size;
 /*alloc memory size*/
 ssize_t alloc;
 /*capacity*/
 ssize_t capacity;


};


typedef struct bucket bucket;



void *init(struct bucket, ssize_t capacity )
{

     bucket* arthox = (bucket*)malloc(sizeof(bucket));
     arthox -> alloc = sizeof(bucket);
     /**/
     arthox -> size = 0;
     arthox -> capacity = capacity;


}

void *fill(struct bucket *bucket, void* data)
{

  /* the init will initialize values */
 memcpy(&bucket -> arr, data, sizeof(data));
 bucket -> size++;


}

void *empty(struct bucket *bucket, void* data)
{


  for(;;)
  {
    if (bucket -> arr == NULL)
    {

      free(bucket);
      break;
    }
  }
}
