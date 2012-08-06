#ifndef __WMNAF_REPR
#define __WMNAF_REPR
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

static signed char* wMNAF(mpz_t x, int w, size_t *ret_len) {
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
        window_val += bit * mpz_tstbit(x, j + w);
    }

    *ret_len = j;
    return ret;
}
#endif
