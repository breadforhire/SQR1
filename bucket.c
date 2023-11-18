#include <stdlib.h>
#include "stdint.h"
#include "string.h"
#include "immintrin.h"
#include "emmintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>
#include <stdio.h>
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
     size_t table_size = sizeof(struct table) - sizeof(__m256i) * 8;
     struct table* arthox = malloc(table_size + sizeof(__m256i) * 8);
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
  free(&table -> bucket[key]);

}

/*mostly for debugging this will be thrown away at the final stage*/
void b_p(struct table *table)
{

 for(int i  = 0; i < table -> size ; i++)
 {

 for(int k  = 0; k < table -> size / 2; k++ )
 {

  printf(" bucket %d --> %llu \n", i, table -> bucket[i][k]);

 }
 printf("%s \n", "______________");
}
}

void *release(struct table* table)
{

 free(table);

}
