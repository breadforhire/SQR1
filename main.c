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

/*The table above will be z ≡ x y(mod t9 − 1)   and will contain our base values to actually run the program*/
/*2, 2, 4
4, 2, 8
8, 2, 16
16, 2, 32
32, 2, 64
64, 2, 128
128, 2, 256
256, 2, 512

256, 4, 1024
2, 4, 8
4, 4, 16
8, 4, 32
16, 4, 64
32, 4, 128
64, 4, 256
128, 4, 512

64, 8, 512
128, 8, 1024
256, 8, 2048
2, 8, 16
4, 8, 32
8, 8, 64
16, 8, 128
32, 8, 256

8, 16, 128
16, 16, 256
32, 16, 512
64, 16, 1024
128, 16, 2048
256, 16, 4096
2, 16, 32
4, 16, 64

128, 32, 4096
256, 32, 8192
2, 32, 64
4, 32, 128
8, 32, 256
16, 32, 512
32, 32, 1024
64, 32, 2048

4, 64, 256
8, 64, 512
16, 64, 1024
32, 64, 2048
64, 64, 4096
128, 64, 8192
256, 64, 16384
2, 64, 128

16, 128, 2048
32, 128, 4096
64, 128, 8192
128, 128, 16384
256, 128, 32768
2, 128, 256
4, 128, 512
8, 128, 1024

32, 256, 8192
64, 256, 16384
128, 256, 32768
256, 256, 65536
2, 256, 512
4, 256, 1024
8, 256, 2048
16, 256, 4096 */




/*we need align each int32_t array*/


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
 __m256i e_x;
 __m256i e_y;
 __m256i z0;
 __m256i z1;
 __m256i buckets[length];


 __m256i constant_two = _mm256_load_si256((__m256i*) twos);
 __m256i constant_sixteen_e2 = _mm256_load_si256((__m256i*) sixteen_e_2 );
 __m256i constant_fours =_mm256_load_si256((__m256i*) fours);



 __m256i t0;
 __m256i t1;
 __m256i t2;






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

    int32_t finalSum = temp[0] + temp[1] + temp[2] + temp[3];
    buckets[i] = _mm256_set1_epi32(finalSum);




 }


 //t1 = buckets[length];
 t1 = _mm256_mullo_epi16(buckets[length], constant_two);
 t1 = _mm256_add_epi16(buckets[length], constant_sixteen_e2);


 t2 = _mm256_mullo_epi16(buckets[length - 1], constant_fours);
 //t2 = _mm256_add_epi64(t2,  constant_fours);



}

void main()
{
 struct table* t = init(length);
 A(t);

}
