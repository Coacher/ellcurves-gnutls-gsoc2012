#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "ecc.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))

/*
 * A local replacement for mpz_tstbit.
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
 * Return an array with wMNAF representationtogether with its length.
 *
 * The result is the array with elements from the set {0, +-1, +-3, +-5, ..., +-(2^w - 1)}
 * such that at most one of any (w + 1) consecutive digits is non-zero
 * with exception for the the most significant (w + 1) bits.
 * With the last property it is modified version of wNAF.
 * Overview of this algorithm can be found, for exmaple, in
 * Bodo Moller, Improved Techniques for Fast Exponentiation.
 * Information Security and Cryptology – ICISC 2002, Springer-Verlag LNCS 2587, pp. 298–312
 */
/*
   @param x        The number to get wMNAF for
   @param w        Window size
   @param len      [out] Destination for the length of wMNAF
   @return         array with wMNAF representation
   @return         NULL in case of errors
 */
signed char* ecc_wMNAF(mpz_t x, unsigned int w, size_t* wmnaf_len) {
    int b, c;
    char sign = 1;
    size_t i, len;

    /* needed constants */
    int basew    = 1 << w;       /* 2^w */
    int baseww   = 1 << (w + 1); /* 2^(w+1) */
    int mask     = baseww - 1;

    signed char *ret = NULL;

    if (!(sign = mpz_sgn(x))) {
        /* x == 0 */
        ret = malloc(1);
        if (ret == NULL)
            goto done;

        ret[0] = 0;
        *wmnaf_len = 1;
        goto done;
    }

    /* total number of bits */
    len = mpz_sizeinbase(x, 2);

    /* wMNAF is at most (len + 1) bits long */
    ret = malloc(len + 1);
    if (ret == NULL)
        goto done;

    /* get (w + 1) Least Significant Bits */
    c = (mpz_getlimbn(x, 0)) & mask;

    /* how many bits we've already processed */
    i = 0;

    while ( (c != 0) || (i + w + 1 < len) ) {
        if (c & 1) {
            /* LSB == 1 */
            if (c >= basew) {
                b = c - baseww;
            } else {
                b = c;
            }

            c -= b;
        } else {
            b = 0;
        }

        ret[i++] = sign * b;

        /* fill c with next LSB */
        c >>= 1;
        c += basew * mpz_unitstbit(x, i + w);
    }

    *wmnaf_len = i--;

    /* modified wNAF */
    if ((ret[i] == 1) && (ret[i-w-1] < 0)) {
        ret[i-w-1] += basew;
        ret[i-1] = 1;
        ret[i] = 0;
        *wmnaf_len = i;
    }
done:
    return ret;
}
