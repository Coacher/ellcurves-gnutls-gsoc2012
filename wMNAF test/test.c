#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "log.c"

/* return array containing wMNAF representation of given mpz_t number.
 * inspired by OpenSSL version of this function.
 * overview of this algorithm can be found in
 * Bodo Moller, Improved Techniques for Fast Exponentiation.
 * Information Security and Cryptology – ICISC 2002, Springer-Verlag LNCS 2587, pp. 298–312
 */
signed char* wMNAF(mpz_t x, int w, size_t *ret_len) {
    int window_val;
    signed char *ret = NULL;
    int sign = 1;
    int bit, next_bit, mask;
    size_t len = 0, j;
    
    if (!mpz_sgn(x)) { 
        ret = malloc(1);
        if (!ret) {
            dbg_msg("Unable to allocate memory in wMNAF.\n");
            exit(EXIT_FAILURE);
        }

        ret[0] = 0;
        *ret_len = 1;
        return ret;
    }
        
    if (w <= 0 || w > 7) {
        dbg_msg("signed char can represent integers with absolute values less than 2^7.\n");
        exit(EXIT_FAILURE);
    }

    bit = 1 << w;        /* at most 128 */
    next_bit = bit << 1; /* at most 256 */
    mask = next_bit - 1; /* at most 255 */

    if (mpz_sgn(x) < 0) {
        sign = -1;
    }

    len = mpz_sizeinbase(x, 2);
    /* modified wMNAF may be one digit longer than binary representation
     * (*ret_len will be set to the actual length, i.e. at most
     * mpz_sizeinbase(x, 2) + 1) */
    ret = malloc(len + 1);

    if (!ret) {
        dbg_msg("Unable to allocate memory in wMNAF.\n");
        exit(EXIT_FAILURE);
    }

    window_val = ((int) mpz_getlimbn(x, 0) ) & mask;
    j = 0;

    /* if j+w+1 >= len, window_val will not increase */
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
            
            if (digit <= -bit || digit >= bit || !(digit & 1)) {
                dbg_msg("Unexpected digit value: \"%i\".\n", digit);
                exit(EXIT_FAILURE);
            }

            window_val -= digit;

            /* now window_val is 0 or 2^(w+1) in standard wNAF generation;
             * for modified window NAFs, it may also be 2^w
             */
            if (window_val != 0 && window_val != next_bit && window_val != bit) {
                dbg_msg("Unexpected window_val value: \"%i\".\n", window_val);
                exit(EXIT_FAILURE);
            }
        }

        ret[j++] = sign * digit;

        window_val >>= 1;
        window_val += bit * mpz_tstbit(x, j + w);

        if (window_val > next_bit) {
            dbg_msg("Unexpected window_val and next_bit values: \"%i\", \"%i\".\n", window_val, next_bit);
            exit(EXIT_FAILURE);
        }
    }

    if (j > len + 1) {
        dbg_msg("Unexpected j and len + 1 values: \"%i\", \"%i\".\n", j, len+1);
        exit(EXIT_FAILURE);
    }
    len = j;
    
    *ret_len = len;
    return ret;
}

int main(void) {
    mpz_t x;
    signed char* wmnaf;
    size_t wmnaf_len;

    int i, j;

    mpz_init(x);

    for (i = 0; i <= 128; ++i) {
        mpz_set_si(x, i);
        wmnaf = wMNAF(x, 3, &wmnaf_len);

        printf("\"%3i\" has wMNAF repr. = [", i);
        for (j = wmnaf_len - 1; j > 0; --j) {
            printf(" %i,", wmnaf[j]);
        }
        printf(" %i ]\n", wmnaf[0]);

        free(wmnaf);
    }

    mpz_clear(x);

    return 0;
}
