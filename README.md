# avx-primes
Faster ECC over F2^521 - 1 (SQR ALGORITHM )




INPUT: x = [x0, . . . , x8] ∈ [−2^59, 2^59 − 1] × [0, 2^58 − 1]^8
OUTPUT: z ∈ [−2^59, 259 − 1] × [0, 2^58 − 1]8
where z ≡ x^2 (mod t^9 − 2)




z ≡ x2(mod t9 − 2)



You can plugin in your own residue/values in the x_n, y_n, and eight make sure x_n = y_n and that they are aligned
























How to use
gcc -mavx -mavx2 main.c
./main

















https://eprint.iacr.org/2014/852.pdf
