#include "immintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>
#include <stdio.h>
#include "bucket.h"
#include <string.h>
/* This is going to a one month project in which we will implement the 64 bit granger moss primes in tables */

#define length 8
/*leaving this empty */
#define c_t 2


#define ROLL(start, end, buckets, vector, bucketindex, storeindex, result) \
    do { \
        __m256i sum = _mm256_setzero_si256(); \
        for (int i = start; i < end; i++) { \
            _mm256_storeu_si256((__m256i*)result, buckets[i]); \
            sum = _mm256_add_epi32(sum, _mm256_set1_epi32(result[bucketindex])); \
        } \
        vector[storeindex] = _mm256_extract_epi32(sum, 0); \
        memset(&result, 0x00, length); \
    } while (0)

#define UNROLL(vector, index, vector_assignment, index_assignment, mod, filler) \
    do { \
        filler = _mm256_and_si256(vector, _mm256_sub_epi32(mod, _mm256_set1_epi32(1))); \
        vector_assignment[index_assignment] = _mm256_extract_epi32(filler, 0); \
    } while (0)

#define PRINT(vector) \
    do { \
        for (int i = 0; i < length; i++) { \
            printf("%d\n", vector[i]); \
        } \
    } while (0)

#define STORE(vector, index , vector_assignment, index_assignment) \
 vector_assignment[index_assignment] = _mm256_extract_epi32(vector, index); \
 vector = _mm256_set1_epi32(0) \

/*
[
x0y0 + x1y1 ...

/*ex*
x0y0  + x1y8 + x2y7 + x3y6 + x4y5 + x5y4 + x6y3 + x7y2 + x8y1,
x0y1 + x1y0 + x2y8 + x3y7 + x4y6 + x5y5 + x6y4 + x7y3 + x8y2,
x0y2 + x1y1 + x2y0 + x3y8 + x4y7 + x5y6 + x6y5 + x7y4 + x8y3,
x0y3 + x1y2 + x2y1 + x3y0 + x4y8 + x5y7 + x6y6 + x7y5 + x8y4,
x0y4 + x1y3 + x2y2 + x3y1 + x4y0 + x5y8 + x6y7 + x7y6 + x8y5,
x0y5 + x1y4 + x2y3 + x3y2 + x4y1 + x5y0 + x6y8 + x7y7 + x8y6,
x0y6 + x1y5 + x2y4 + x3y3 + x4y2 + x5y1 + x6y0 + x7y8 + x8y7,
x0y7 + x1y6 + x2y5 + x3y4 + x4y3 + x5y2 + x6y1 + x7y0 + x8y8,
x0y8 + x1y7 + x2y6 + x3y5 + x4y4 + x5y3 + x6y2 + x7y1 + x8y0].  */




__attribute__((__aligned__(32)))
static const int32_t x_n[8] = {2, 4, 8, 16, 32, 64, 128, 256};

__attribute__((__aligned__(32)))
static int32_t y_n [8] = {2, 4, 8, 16, 32, 64, 128, 256};

__attribute__((__aligned__(32)))
static const int32_t twos [8] = {2, 2, 2, 2, 2, 2, 2, 2};

__attribute__((__aligned__(32)))
static const int32_t fours [8] = {4, 4, 4, 4, 4, 4, 4, 4};

/*this provides better efficiency */
__attribute__((__aligned__(32)))
static const int32_t sixteen_e_2 [8] = {256, 256, 256, 256, 256, 256, 256, 256};

__attribute__((__aligned__(32)))
static int32_t shift_mask [8] = {0, 1, 2, 3, 4, 5, 6, 7 };


int32_t temp[length];

void shiftArray(int arr[], int size, int shiftBy)
{

    shiftBy %= size;
    int buffer[shiftBy];
    memcpy(buffer, arr + (size - shiftBy), sizeof(int) * shiftBy);
    memmove(arr + shiftBy, arr, sizeof(int) * (size - shiftBy));
    memcpy(arr, buffer, sizeof(int) * shiftBy);
}


void A(struct table *table)
{




 /*load resiude*/
 int32_t result[length];
 __m256i e_x;
 __m256i e_y;
 __m256i z0;
 __m256i z1;
 __m256i buckets[length];

 __m256i filler = _mm256_set1_epi32(1);
 __m256i constant_two = _mm256_load_si256((__m256i*) twos);
 __m256i constant_sixteen_e2 = _mm256_load_si256((__m256i*) sixteen_e_2 );
 __m256i constant_fours =_mm256_load_si256((__m256i*) fours);



 /*vectors that will be applied but not used for storage*/
 __m256i t0;
 __m256i t1;


 __m256i t0_storage;
 __m256i t1_storage;
 __m256i t2_storage;
 __m256i z0_storage;






 /* bucket the entire x0y0 + x1y8 + x2y7 + x3y6 + x4y5 + x5y4 + x6y3 + x7y2 + x8y1  as*(x0y8 x1y7 x2y6 x3y5) +  (x4y4  x5y3  x6y2 x7y1  x8y0)  */
 /*shift then add to the bucket*/


 for(int i = 0; i < length; i++)
 {


                                    /*do this with every table*/
                                    /*x7y1   x8y0   x3y5*/
  /*fill every bucket with vector   [x0y8 + x1y7 + x2y6 + x3y5 + x4y4 + x5y3 + x6y2 + x7y1 + x8y0] */

   /*shift once permutations are fixed so this is fine*/

    e_x = _mm256_set1_epi32(x_n[i]);
    z0 = _mm256_mullo_epi32(_mm256_load_si256((__m256i*) y_n) , _mm256_set1_epi32(x_n[i]));
    shiftArray(y_n, length, i);
    _mm256_storeu_si256((__m256i*)temp, z0);
    shiftArray(temp, length, i);
     z1 = _mm256_loadu_si256((__m256i*)temp);

    _mm256_storeu_si256(&buckets[i], z1);






 }
 /*t1 = 2(x0x8 + x1x7 + x2x6 + x3x5) + x4 ^2*/
 /*start, end, buckets, vector, bucket_index, storeindex arr*/
 ROLL(0, 4, buckets, t1, length - 1, 0, result);
 t1 = _mm256_mullo_epi32(t1, constant_two);
 t1 = _mm256_add_epi32(t1, constant_sixteen_e2);
 /*vector, index , vector_assignment, index_assignment*/
 STORE(t1, 0, t1_storage, 0);


/*vector, index, vector_assignment, index_assignment, mod*/
/*takes the mod of t1[0] by constant two and stores it in t0_storage at index 0*/
/*t0 = t1 mod t*/
/*the storage is already built into unroll so we do not have to worry about it*/
 UNROLL(t1_storage, 0, t0_storage, 0, constant_two, filler);


 /*t2 = 4(x1x8 + x2x7 + x3x6 + x4x5) + x0^2 + 2(t1 >> 58)*/
 ROLL(1, 5, buckets, t1, 0, 1, result);
 t1 = _mm256_mullo_epi32(t1, constant_fours);
 t1 = _mm256_add_epi32(t1,  constant_fours);
 t1 = _mm256_add_epi32(t1, _mm256_mullo_epi32(_mm256_set1_epi32(_mm256_extract_epi32(t1_storage, 0) >> 58), constant_two)) ;

 STORE(t1, 2, t2_storage, 0);
 UNROLL(t2_storage, 0, z0, 0, constant_two, filler);

 PRINT(z0);

 ROLL(3, 6, buckets, t1, 1, 0, result);
 t1 = _mm256_mullo_epi32(t1, constant_fours);
 ROLL(0, 1, buckets, t1, 1, 1, result);
 t1 = _mm256_add_epi32(t1, _mm256_mullo_epi32(_mm256_set1_epi32(&t1[1]), constant_two));







}
int main()
{

 __m256i buckets[length];

 A(buckets);
 return 0;

}
