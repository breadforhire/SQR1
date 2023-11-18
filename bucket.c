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




struct table* init(ssize_t capacity)
{

     struct table* arthox = (table*)malloc(sizeof(table));
     arthox -> alloc = sizeof(table);
     arthox -> size = 0;
     arthox -> size = capacity;
     return arthox;


}

void* fill(struct table* table, __m256i data, int key)
{


  table -> bucket[key] =  data;
  table -> size++;

}


__m256i extract(struct table *table, int key)
{

  return table -> bucket[key];

}
