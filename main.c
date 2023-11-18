#include "immintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>
#include <stdio.h>
#include "bucket.h"
/* This is going to a one month project in which we will implement the 64 bit granger moss primes in tables */

#define length 8

/*
[
x0y0 + x1y1 ...


x0y0 + x1y8 + x2y7 + x3y6 + x4y5 + x5y4 + x6y3 + x7y2 + x8y1,
x0y1 + x1y0 + x2y8 + x3y7 + x4y6 + x5y5 + x6y4 + x7y3 + x8y2,
x0y2 + x1y1 + x2y0 + x3y8 + x4y7 + x5y6 + x6y5 + x7y4 + x8y3,
x0y3 + x1y2 + x2y1 + x3y0 + x4y8 + x5y7 + x6y6 + x7y5 + x8y4,
x0y4 + x1y3 + x2y2 + x3y1 + x4y0 + x5y8 + x6y7 + x7y6 + x8y5,
x0y5 + x1y4 + x2y3 + x3y2 + x4y1 + x5y0 + x6y8 + x7y7 + x8y6,
x0y6 + x1y5 + x2y4 + x3y3 + x4y2 + x5y1 + x6y0 + x7y8 + x8y7,
x0y7 + x1y6 + x2y5 + x3y4 + x4y3 + x5y2 + x6y1 + x7y0 + x8y8,
x0y8 + x1y7 + x2y6 + x3y5 + x4y4 + x5y3 + x6y2 + x7y1 + x8y0].  */

/*The table above will be z ≡ x y(mod t9 − 1)   and will contain our base values to actually run the program*/
 /*the buckets are filled in reverse for example the first index is x0y8 + x1y7 + x2y6 + x3y5 + x4y4 + x5y3 + x6y2 + x7y1 + x8y0 */





/*we need align each int32_t array*/


__attribute__((__aligned__(32)))
static const int32_t x_n[8] = {2, 4, 8, 16, 32, 64, 128, 256};

__attribute__((__aligned__(32)))
static const int32_t y_n [8] = {2, 4, 8, 16, 32, 64, 128, 256};

__attribute__((__aligned__(32)))
static const int32_t twos [8] = {2, 2, 2, 2, 2, 2, 2, 2};

void A(struct table *table)
{




 /*load resiude*/
 __m256i e_x = _mm256_load_si256((__m256i*) x_n);
 __m256i e_y = _mm256_load_si256((__m256i*) y_n);
 __m256i z0;
 __m256i z1;

 __m256i constant_two = _mm256_load_si256((__m256i*) twos);

 __m256i t;






 /* bucket the entire x0y0 + x1y8 + x2y7 + x3y6 + x4y5 + x5y4 + x6y3 + x7y2 + x8y1  as*(x0y8 x1y7 x2y6 x3y5) +  (x4y4  x5y3  x6y2 x7y1  x8y0)  */
 /*shift then add to the bucket*/

 for(int i = 0; i < length; ++i)
 {

   z0 = _mm256_mullo_epi16(e_x, e_y);
   z1 = _mm256_mulhi_epi16(e_x, e_y);
   
   /*shifts the entire vector*/
   e_y = _mm256_slli_epi64(e_y, 1) ;

  /*fill every bucket with vector  (x0y8 x1y7 x2y6 x3y5) +  (x4y4  x5y3  x6y2 x7y1  x8y0)  */
  fill(table, _mm256_add_epi64(z0, z1), i);
 }

 /*the buckets are filled in reverse for example the first index is x0y8 + x1y7 + x2y6 + x3y5 + x4y4 + x5y3 + x6y2 + x7y1 + x8y0 */
 /*t1 = 2(x0x8 + x1x7 + x2x6 + x3x5) + x4^2 */

 t = _mm256_mulhi_epi16(constant_two, extract(table, length % length));

}

void main()
{
 struct table* t = init(length);
 A(t);

}
