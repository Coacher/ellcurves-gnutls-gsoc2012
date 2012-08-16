#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "ecc.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))

/* A local replacement for mpz_tstbit.
 * It is needed because for negative numbers original mpz_tstbit
 * returns an infinite number of `1`s after all bits of input number.
 * For positive numbers it returns zeros after all bits of input number.
 * This function mimics mpz_tstbit behavior for positive numbers in both cases.
 */
static int
mpz_unitstbit (mpz_srcptr u, mp_bitcnt_t bit_index) __GMP_NOTHROW
{
  mp_srcptr      u_ptr      = (u)->_mp_d;
  mp_size_t      size       = (u)->_mp_size;
  unsigned       abs_size   = ABS(size);
  mp_size_t      limb_index = bit_index / GMP_NUMB_BITS;
  mp_srcptr      p          = u_ptr + limb_index;
  mp_limb_t      limb;

  if (limb_index >= abs_size)
    return (size < 0);

  limb = *p;

  return (limb >> (bit_index % GMP_NUMB_BITS)) & 1;
}

/*
 * Return an array with wMNAF representation of given mpz_t number x
 * together with representation length.
 * The result is the array with elements from the set {0, +-1, +-3, +-5, ..., +-(2^w - 1)}
 * such that at most one of any (w + 1) consecutive digits is non-zero
 * with exception for the the most significant (w+1) bits.
 * With the last property it is modified version of wNAF.
 * Based on OpenSSL version of this function.
 * Overview of this algorithm can be found in
 * Bodo Moller, Improved Techniques for Fast Exponentiation.
 * Information Security and Cryptology – ICISC 2002, Springer-Verlag LNCS 2587, pp. 298–312
 */

signed char* ecc_wMNAF(mpz_t x, unsigned int w, size_t *ret_len) {
#if 0
    int window_val;
    signed char *ret = NULL;
    int sign = 1;
    int bit, next_bit, mask;
    size_t len = 0, j;

    if (!(sign = mpz_sgn(x))) {
        /* x == 0 */
        ret = malloc(1);
        if (!ret) {
            fprintf(stderr, "Unable to allocate memory for wMNAF repr.\n");
            return NULL;
        }

        ret[0] = 0;
        *ret_len = 1;
        return ret;
    }

#ifdef __WITH_EXTRA_SANITY_CHECKS
    if (w <= 0 || w > 7) {
        fprintf(stderr, "`signed char` can represent integers with absolute values less than 2^7.\n");
        return NULL;
    }
#endif

    /* 2^w, at most 128 */
    bit = 1 << w; 
    /* 2^(w + 1), at most 256 */
    next_bit = bit << 1;
    /* 2^(w + 1) - 1, at most 255 */
    mask = next_bit - 1;

    len = mpz_sizeinbase(x, 2);
    /* wMNAF may be one digit longer than binary representation
     * (*ret_len will be set to the actual length, i.e. at most
     * mpz_sizeinbase(x, 2) + 1) */
    ret = malloc(len + 1);

    if (!ret) {
        fprintf(stderr, "Unable to allocate memory for wMNAF repr.\n");
        return NULL;
    }

    window_val = (mpz_getlimbn(x, 0)) & mask;
    j = 0;

    /* if (j + w + 1) >= len, window_val will not increase */
    while ((window_val != 0) || (j + w + 1 < len)) {
        int digit = 0;

        /* 0 <= window_val <= 2^(w+1) */

        if (window_val & 1) {
            /* 0 < window_val < 2^(w+1) */

            if (window_val & bit) {
                digit = window_val - next_bit; /* -2^w < digit < 0 */

                if (j + w + 1 >= len) {
                    /* special case for generating modified wNAFs:
                     * no new bits will be added into window_val,
                     * so using a positive digit here will decrease
                     * the total length of the representation */

                    digit = window_val & (mask >> 1); /* 0 < digit < 2^w */
                }
            } else {
                digit = window_val; /* 0 < digit < 2^w */
            }

            window_val -= digit;
        }

        ret[j++] = sign * digit;

        window_val >>= 1;
        window_val += bit * mpz_unitstbit(x, j + w);
    }

    *ret_len = j;
    return ret;
#endif
    signed char *ret = NULL;
    int sign = 1;
    int mask = (1 << (w + 1)) - 1;
    size_t len = 0, j;
    mpz_t c;

    if (!(sign = mpz_sgn(x))) {
        /* x == 0 */
        ret = malloc(1);
        if (!ret) {
            fprintf(stderr, "Unable to allocate memory for wMOF repr.\n");
            return NULL;
        }

        ret[0] = 0;
        *ret_len = 1;
        return ret;
    }

    len = mpz_sizeinbase(x, 2);
    /* wMNAF may be one digit longer than binary representation
     * (*ret_len will be set to the actual length, i.e. at most
     * mpz_sizeinbase(x, 2) + 1) */
    ret = malloc(len + 1);

    if (!ret) {
        fprintf(stderr, "Unable to allocate memory for wMNAF repr.\n");
        return NULL;
    }

    mpz_init(c);
    mpz_abs(c, x);
    j = 0;

    while( j > w ) {

        if ( bit == next_bit ) {
            ret[j] = 0;
            --j;
        } else {
            j -= w;
        }

        next_bit = bit;
        bit = mpz_unitstbit(x, j);
}
