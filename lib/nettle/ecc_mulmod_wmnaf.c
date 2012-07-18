/*
 * Copyright (C) 2011-2012 Free Software Foundation, Inc.
 *
 * This file is part of GNUTLS.
 *
 * The GNUTLS library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#include "ecc.h"

/* size of sliding window, don't change this! */
#ifndef WINSIZE
    #define WINSIZE 4
#endif

/* length of one array of precomputed values for ecc_mulmod_wmnaf 
 * we have two such arrays for positive and negative multipliers */
#define PRECOMPUTE_LENGTH (1 << (WINSIZE - 1))

/* returns an array with wMNAF representation of given mpz_t number x
 * together with length of the representation */
static signed char* wMNAF(mpz_t x, int w, size_t *ret_len);

/*
   Perform a point multiplication using wMNAF repr.
   @param k    The scalar to multiply by
   @param G    The base point
   @param R    [out] Destination for kG
   @param modulus  The modulus of the field the ECC curve is in
   @param map      Boolean whether to map back to affine or not (1 == map, 0 == leave in projective)
   @return CRYPT_OK on success
*/
int
ecc_mulmod_wmnaf (mpz_t k, ecc_point * G, ecc_point * R, mpz_t a, mpz_t modulus,
                int map)
{
    ecc_point *pos[PRECOMPUTE_LENGTH], *neg[PRECOMPUTE_LENGTH];
    int        i, j, err;

    signed char* wmnaf = NULL;
    size_t wmnaf_len;
    signed char digit;

    if (k == NULL || G == NULL || R == NULL || modulus == NULL)
        return -1;

    /* alloc ram for precomputed values */
    for (i = 0; i < PRECOMPUTE_LENGTH; ++i) {
        pos[i] = ecc_new_point();
        neg[i] = ecc_new_point();
        if (pos[i] == NULL || neg[i] == NULL) {
            for (j = 0; j < i; ++j) {
              ecc_del_point(pos[j]);
              ecc_del_point(neg[j]);
            }

            return -1;
        }
    }

    /* fill in pos and neg arrays with precomputed values
     * pos holds kG for k ==  1, 3, 5, ..., (2^w - 1)
     * neg holds kG for k == -1,-3,-5, ...,-(2^w - 1)
     */

    /* pos[0] == 2G for a while, later it will be set to the expected 1G */
    if ((err = ecc_projective_dbl_point(G, pos[0], a, modulus)) != 0)
        goto done;
   
    /* pos[1] == 3G */
    if ((err = ecc_projective_add_point(pos[0], G, pos[1], a, modulus)) != 0)
        goto done;

    /* fill in kG for k = 5, 7, ..., (2^w - 1) */
    for (j = 2; j < PRECOMPUTE_LENGTH; ++j) {
        if ((err = ecc_projective_add_point(pos[j-1], pos[0], pos[j], a, modulus)) != 0)
           goto done;
    }
   
    /* set pos[0] == 1G as expected
     * after this step we don't need G at all 
     * and can change it without worries even if R == G */
    mpz_set (pos[0]->x, G->x);
    mpz_set (pos[0]->y, G->y);
    mpz_set (pos[0]->z, G->z);

    /* neg[i] == -pos[i] */
    for (j = 0; j < PRECOMPUTE_LENGTH; ++j) {
        if ((err = ecc_projective_negate_point(pos[j], neg[j], modulus)) != 0)
            goto done;
    }

    /* calculate wMNAF */
    wmnaf = wMNAF(k, WINSIZE, &wmnaf_len);
    if (!wmnaf) {
        err = -2;
        goto done;
    }

    /* actual point computation */

    /* set R to neutral */
    mpz_set_ui(R->x, 1);
    mpz_set_ui(R->y, 1);
    mpz_set_ui(R->z, 0);
   
    /* perform ops */
    for (j = wmnaf_len - 1; j >= 0; --j) {
        if ((err = ecc_projective_dbl_point(R, R, a, modulus)) != 0)
            goto done;

        digit = wmnaf[j];

        if (digit) {
            if (digit > 0) {
                if ((err = ecc_projective_add_point(R, pos[( digit / 2)], R, a, modulus)) != 0)
                    goto done;
            } else {
                if ((err = ecc_projective_add_point(R, neg[(-digit / 2)], R, a, modulus)) != 0)
                    goto done;
            }
        }
    }


    /* map R back from projective space */
    if (map) {
        err = ecc_map(R, modulus);
    } else {
        err = 0;
    }
done:
    for (i = 0; i < PRECOMPUTE_LENGTH; ++i) {
        ecc_del_point(pos[i]);
        ecc_del_point(neg[i]);
    }
    if (wmnaf) free(wmnaf);
    return err;
}


/* 
 * Return an array with wMNAF representation of given mpz_t number x
 * together with representation length.
 * The result is the array with elements from set {0, +-1, +-3, +-5, ..., +-(2^w - 1)}
 * such that at most one of any (w + 1) consecutive digits is non-zero
 * with exception for the the most significant (w+1) bits.
 * With the last property it is modified version of wNAF.
 * Based on OpenSSL version of this function.
 * overview of this algorithm can be found in
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
