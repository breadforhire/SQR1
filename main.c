#include "immintrin.h"
#include "avx2intrin.h"
#include <smmintrin.h>
#include <stdio.h>
#include <string.h>

#define length 8

/*This is for 64 bit, do not use for 32 bit*/
#define c_t 288230376151712000


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


#define PRINT(vector) \
    do { \
        for (int i = 0; i < 8; i++) { \
            printf("\n%d --> %d \n", vector[i], i); \
        } \
    } while (0)

#define PRINT_Z(vector) \
    do { \
        for (int i = 0; i < 4; i++) { \
            printf("\n%d --> %d \n", vector[i], i); \
        } \
    } while (0)






/*this residue or starting values are wrong please do not use as it was only used for tests and benchmarking*/
__attribute__((__aligned__(32)))
static const int32_t x_n[8] = {0, 1, 2, 3, 4, 5, 6, 7};

/*this residue or starting values are wrong please do not use as it was only used for tests and benchmarking*/
__attribute__((__aligned__(32)))
static int32_t y_n [8] = {0, 1, 2, 3, 4, 5, 6, 7};

__attribute__((__aligned__(32)))
static const int32_t twos [8] = {2, 2, 2, 2, 2, 2, 2, 2};

__attribute__((__aligned__(32)))
static const int32_t fours [8] = {4, 4, 4, 4, 4, 4, 4, 4};

/*this x_8 will act as our last value to the 7 elements*/
__attribute__((__aligned__(32)))
static const int32_t x_8_c[8] = {8, 8, 8, 8, 8, 8, 8, 8};

/*this provides better efficiency */
__attribute__((__aligned__(32)))
static const int32_t sixteen_e_2 [8] = {256, 256, 256, 256, 256, 256, 256, 256};

__attribute__((__aligned__(32)))
static int32_t shift_mask [8] = {0, 1, 2, 3, 4, 5, 6, 7 };



int32_t temp[length];

void shiftArray(int32_t arr[], int size, int shiftBy)
{

    shiftBy %= size;
    int buffer[shiftBy];
    memcpy(buffer, arr + (size - shiftBy), sizeof(int) * shiftBy);
    memmove(arr + shiftBy, arr, sizeof(int) * (size - shiftBy));
    memcpy(arr, buffer, sizeof(int) * shiftBy);
}


void A()
{




 /*load resiude*/
 int32_t result[length ];
 __m256i e_x;
 __m256i e_y;
 __m256i temp0;
 __m256i temp1;
 __m256i buckets[length + 1];


 __m256i constant_two = _mm256_load_si256((__m256i*) twos);
 __m256i constant_sixteen_e2 = _mm256_load_si256((__m256i*) sixteen_e_2 );
 __m256i constant_fours =_mm256_load_si256((__m256i*) fours);
 __m256i x_8 = _mm256_load_si256((__m256i*) x_8_c);



 /*vectors that will be applied but not used for storage*/
 __m256i t0;
 __m256i t1;
 __m256i t2;



 __m256i t0_storage;


 __m256i z0;
 __m256i z1;
 __m256i z2;
 __m256i z3;
 __m256i z4;
 __m256i z5;
 __m256i z6;
 __m256i z7;
 __m256i z8;

 __m256i filler;







 for(int i = 0; i < length; i++)
 {


                                    /*do this with every table*/
                                    /*x7y1   x8y0   x3y5*/
  /*fill every bucket with vector   [x0y8 + x1y7 + x2y6 + x3y5 + x4y4 + x5y3 + x6y2 + x7y1 + x8y0] */

   /*shift once permutations are fixed so this is fine*/

    e_x = _mm256_set1_epi32(x_n[i]);
    temp0 = _mm256_mullo_epi32(_mm256_load_si256((__m256i*) y_n) , _mm256_set1_epi32(x_n[i]));
    shiftArray(y_n, length, 1);
    _mm256_storeu_si256((__m256i*)temp, temp0);
    temp1 = _mm256_loadu_si256((__m256i*)temp);

    _mm256_storeu_si256(&buckets[i], temp1);


 }

 /*t1 = 2(x0x8 + x1x7 + x2x6 + x3x5) + x4 ^2 --> this works*/
 ROLL(1, 4, buckets, t1, 0, 0, result);
 t1 = _mm256_add_epi32(_mm256_mullo_epi32(x_8, _mm256_set1_epi32(x_n[0])), t1);
 t1 = _mm256_mullo_epi32(t1, constant_two);
 t1 = _mm256_add_epi32(t1, _mm256_set1_epi32(x_n[4] * x_n[4]));

  _mm256_storeu_si256((__m256i*)temp, t1);
 t0_storage = _mm256_set1_epi32(temp[0] % c_t);




 /*t2 = 4(x1x8 + x2x7 + x3x6 + x4x5) + x0^2 + 2(t1 >> 58)*/
 ROLL(2, 5, buckets, t2, 1, 0, result);
 t2 = _mm256_add_epi32(t2, _mm256_mullo_epi32(_mm256_set1_epi32(x_n[1]), x_8));
 t2 = _mm256_mullo_epi32(t2, constant_fours);
 t2 = _mm256_add_epi32(t2,  _mm256_set1_epi32(x_n[0] * x_n[0]));
 t2 = _mm256_add_epi32(t2, _mm256_mullo_epi32(_mm256_set1_epi32(temp[0] >> 58), constant_two)) ;

  _mm256_storeu_si256((__m256i*)temp, t2);
 z0 = _mm256_set1_epi32(temp[0] % c_t);



/*t1 = 4(x2x8 + x3x7 + x4x6) + 2(x0x1 + x5 ^ 2) + (t2 >> 58)*/
 ROLL(3, 5, buckets, t1, 2, 0, result);
 t1 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_set1_epi32(x_n[2]), x_8), t1);
 t1 = _mm256_mullo_epi32(t1, constant_fours);

 ROLL(0, 1, buckets, t0, 1, 0, result);
 t0 = _mm256_add_epi32(t0, _mm256_set1_epi32(x_n[5] * x_n[5]));
 t0 = _mm256_mullo_epi32(t0, constant_two);


 t1 = _mm256_add_epi32(t1, t0);
 t1 = _mm256_add_epi32(t1, _mm256_set1_epi32(temp[0]  >> 58));

 _mm256_storeu_si256((__m256i*)temp, t1);
 z1 = _mm256_set1_epi32(temp[0] % c_t);



 /*t2 = 4(x3x8 + x4x7 + x5x6) + 2(x0x2) + x1 ^ 2 + (t1 >> 58)*/
 ROLL(4, 6, buckets, t2, 3, 0 , result );
 t2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_set1_epi32(x_n[3]), x_8), t2);
 t2 = _mm256_mullo_epi32(t2, constant_fours);

 ROLL(0, 1, buckets, t0, 2, 0, result);
 t0 = _mm256_mullo_epi32(t0, constant_two);
 t2 = _mm256_add_epi32(t0, t2);
 t2 = _mm256_add_epi32(t2, _mm256_set1_epi32((x_n[1] * x_n[1])));
 t2 = _mm256_add_epi32(t2, _mm256_set1_epi32(temp[0]  >> 58));

 _mm256_storeu_si256((__m256i*)temp, t2);
 z2 = _mm256_set1_epi32(temp[0] % c_t);




 /*t1 = 4(x4x8 + x5x7) + 2(x0x3 + x1x2 + x6 ^ 2) + (t2 >> 58)*/
 ROLL(5, 6, buckets, t1, 4, 0, result);
 t1 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_set1_epi32(x_n[4]), x_8), t1);
 t1 = _mm256_mullo_epi32(t1, constant_fours);

 ROLL(0, 2, buckets, t2, 3, 0, result);
 t0 = _mm256_set1_epi32(x_n[6] * x_n[6]);
 t2 = _mm256_add_epi32(t2, t0);
 t2 = _mm256_mullo_epi32(t2, constant_two);

 t1 = _mm256_add_epi32(t2, t1);


 t1 = _mm256_add_epi32(_mm256_set1_epi32(temp[0]  >> 58), t1);

  _mm256_storeu_si256((__m256i*)temp, t1);
  z3 = _mm256_set1_epi32(temp[0] % c_t);



 /* t2 = 4(x5x8 + x6x7) + 2(x0x4 + x1x3) + x2^2 + (t1 >> 58)*/
 ROLL(6, 7, buckets, t2, 5, 0, result);
 t2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_set1_epi32(x_n[5]), x_8), t2);
 t2 = _mm256_mullo_epi32(t2, constant_fours);
 ROLL(0, 2, buckets, t1, 4, 0, result);
 t1 = _mm256_mullo_epi32(t1, constant_two);

 t2 = _mm256_add_epi32(t1, t2);
 t2 = _mm256_add_epi32(t2, _mm256_set1_epi32(x_n[2] * x_n[2]));
 t2 = _mm256_add_epi32(_mm256_set1_epi32(temp[0] >> 58), t2 );

  _mm256_storeu_si256((__m256i*)temp, t2);
  z4 = _mm256_set1_epi32(temp[0] % c_t);




 /*t1 = 4 (x6x8) + 2(x0x5 + x1x4 + x2x3 + x7 ^ 2) + (t2 >> 58)*/
 t1 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_set1_epi32(x_n[6]), x_8), _mm256_set1_epi32(0));
 t1 = _mm256_mullo_epi32(t1, constant_fours);


 ROLL(0, 3, buckets, t2, 5 , 0, result);
 t2 = _mm256_add_epi32(t2, _mm256_set1_epi32(x_n[7] * x_n[7]));
 t2 = _mm256_mullo_epi32(t2, constant_two);


 t1 = _mm256_add_epi32(t1, t2);
 t1 = _mm256_add_epi32(_mm256_set1_epi32(temp[0]  >> 58), t1);

 _mm256_storeu_si256((__m256i*)temp, t1);
 z5 = _mm256_set1_epi32(temp[0] % c_t);

 




 /*t2 = 4 (x7x8) + 2(x0x6 + x1x5 + x2x4) + x3  ^ 2 + (t1 >> 58)*/
 t0 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_set1_epi32(x_n[7]), x_8), _mm256_set1_epi32(0));
 t0 = _mm256_mullo_epi32(t0, constant_fours);

 ROLL(0, 3, buckets, t1, 6, 0, result);
 t1 = _mm256_mullo_epi32(t1, constant_two);

 ROLL(0, 1, buckets, t2, 6, 0, result);

 t2 = _mm256_add_epi32(t0, t1);
 t2 = _mm256_add_epi32(t2, _mm256_set1_epi32(x_n[3] * x_n[3]));

 t2 = _mm256_add_epi32(_mm256_set1_epi32( temp[0] >> 58), t2);

  _mm256_storeu_si256((__m256i*)temp, t2);
 z6 = _mm256_set1_epi32(temp[0] % c_t);



 /*t1 = 2(x0x7 + x1x6 + x2x5 + x3x4 + x8 ^ 2) + (t2 >> 58)*/
 ROLL(0, 4, buckets, t1, 7, 0, result);
 t1 = _mm256_add_epi32(t1, _mm256_mullo_epi32(x_8, x_8));
 t1 = _mm256_mullo_epi32(t1, constant_two);
 t1 = _mm256_add_epi32(_mm256_set1_epi32(temp[0] >> 58), t1);

  _mm256_storeu_si256((__m256i*)temp, t1);
  z7 = _mm256_set1_epi32(temp[0] % c_t);



 __m256i extracted_t0 = _mm256_set1_epi32(_mm256_extract_epi32(t0_storage, 0));
 __m256i extracted_t1 = _mm256_set1_epi32(_mm256_extract_epi32(t1, 0) >> 58);


 t2 = _mm256_add_epi32(t0_storage, _mm256_set1_epi32(temp[0]));
 _mm256_storeu_si256((__m256i*)temp, t2);
 z8 = _mm256_set1_epi32(temp[0] % c_t);



 z0 = _mm256_add_epi32(_mm256_set1_epi32(_mm256_extract_epi32(z0, 0)), _mm256_set1_epi32((temp[0] >> 58) * 2));


 PRINT_Z(z0);
 printf("\n %s \n", "________________");
 PRINT_Z(z1);
 printf("\n%s\n", "________________");
 PRINT_Z(z2);
 printf("\n%s\n", "________________");
 PRINT_Z(z3);
 printf("\n %s \n", "________________");
 PRINT_Z(z4);
 printf("\n%s\n", "________________");
 PRINT_Z(z5);
 printf("\n%s\n", "________________");
 PRINT_Z(z6);
 printf("\n %s \n", "________________");
 PRINT_Z(z7);
 printf("\n%s\n", "________________");
 PRINT_Z(z8);








}
int main()
{


 A();
 return 0;

}
