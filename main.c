#include "immintrin.h"
#include "emmintrin.h"
#include "avx2intrin.h"
#include "bucket.h"
/* This is going to a one month project in which we will implement the 64 bit granger moss primes in tables */

#define _LEN 8


#define ADD(A, B) _mm256_add_epi32(A, B)
#define SUB(A, B) _mm256_sub_epi32(A, B)
#define XOR(A, B) _mm256_xor_si256(A, B)


/*
[
x0y0 + x1y8 + x2y7 + x3y6 + x4y5 + x5y4 + x6y3 + x7y2 + x8y1,
x0y1 + x1y0 + x2y8 + x3y7 + x4y6 + x5y5 + x6y4 + x7y3 + x8y2,
x0y2 + x1y1 + x2y0 + x3y8 + x4y7 + x5y6 + x6y5 + x7y4 + x8y3,
x0y3 + x1y2 + x2y1 + x3y0 + x4y8 + x5y7 + x6y6 + x7y5 + x8y4,
x0y4 + x1y3 + x2y2 + x3y1 + x4y0 + x5y8 + x6y7 + x7y6 + x8y5,
x0y5 + x1y4 + x2y3 + x3y2 + x4y1 + x5y0 + x6y8 + x7y7 + x8y6,
x0y6 + x1y5 + x2y4 + x3y3 + x4y2 + x5y1 + x6y0 + x7y8 + x8y7,
x0y7 + x1y6 + x2y5 + x3y4 + x4y3 + x5y2 + x6y1 + x7y0 + x8y8,
x0y8 + x1y7 + x2y6 + x3y5 + x4y4 + x5y3 + x6y2 + x7y1 + x8y0].  */



/*1 /  8 */

#define M(z0, z1, z2, z3, z4, z5, z6, z7, z8) \
 MUL(z0, z1), MUL(z2, z3)              \
 MUL(z3, z4), MUL(z5, z6)              \
 MUL(z7, z8)                           \

#define A(z0, z1, z2, z3, z4, z5, z6, z7, z8) \
 ADD(z0, z1), ADD(z2, z3)              \
 ADD(z3, z4), ADD(z5, z6)              \
 Add(z7, z8)                           \

#define P(R) R

/* alignment */
__attribute__((__aligned__(32)))
static const int32_t r_y[8] = {2, 5, 17, 37, 101, 197, 257, 401};
/* alignment */
__attribute__((__aligned__(32)))
static const int32_t r_x [8] = {1297, 1601, 2917, 3137, 4357, 5477, 7057, 8101};


void E()
{



  __m256i e_x = _mm256_load_si256((__m256i*) r_x);
 __m256i e_y = _mm256_load_si256((__m256i*) r_y);
 __m256 arthox;

                                     // --> bucket
 /*x0y0 + x1y8 + x2y7 + x3y6 + x4y5 + x5y4 + x6y3 + x7y2 + x8y1 */

 for(;;)
 {
  /* fill up buckets and call it a day */
  __m256i z0 = _mm256_mulhi_epi16(e_x, e_y);
 _mm256_permutevar8x32_ps((__m256) e_x, _mm256_set_epi32(0,7,6,5,4,3,2,1)  );

 }

}
void main()
{

 E();
}
