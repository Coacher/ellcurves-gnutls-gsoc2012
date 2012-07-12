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

/* length of (half)array of precomputed values for wMNAF */
#define PRECOMPUTE_LENGTH_SMALL (1 << (WINSIZE - 1))

/* returns wMNAF representation of given mpz_t number x */
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
    ecc_point *tG, *pos[PRECOMPUTE_LENGTH_SMALL], *neg[PRECOMPUTE_LENGTH_SMALL];
    int        i, j, err;

    signed char* wmnaf = NULL;
    size_t wmnaf_len;
    signed char digit;

    if (k == NULL || G == NULL || R == NULL || modulus == NULL)
        return -1;

    /* alloc ram for precomputed values */
    for (i = 0; i < PRECOMPUTE_LENGTH_SMALL; i++) {
        pos[i] = ecc_new_point();
        neg[i] = ecc_new_point();
        if (pos[i] == NULL || neg[i] == NULL) {
            for (j = 0; j < i; j++) {
              ecc_del_point(pos[j]);
              ecc_del_point(neg[j]);
            }

            return -1;
        }
    }

    /* make a copy of G in case R == G */
    tG = ecc_new_point();
    if (tG == NULL)
    { 
        err = -1;
        goto done; 
    }

    /* tG = G  and convert to montgomery */
    mpz_set (tG->x, G->x);
    mpz_set (tG->y, G->y);
    mpz_set (tG->z, G->z);

    /* 
     * calculate the pos and neg arrays
     * pos holds kG for k ==  1, 3, ..., (2^w - 1)
     * neg holds kG for k == -1,-3, ...,-(2^w - 1)
     */

    /* pos[0] == 2G for a while, later it will be set to expected 1G */
    if ((err = ecc_projective_dbl_point(tG, pos[0], a, modulus)) != 0)
        goto done;
   
    /* pos[1] == 3G */
    if ((err = ecc_projective_add_point(pos[1], tG, pos[0], a, modulus)) != 0)
        goto done;

    /* find kG for k = 5,7, ..., (2^w - 1) */
    for (j = 2; j < PRECOMPUTE_LENGTH_SMALL; ++j) {
        if ((err = ecc_projective_add_point(pos[j], pos[0], pos[j-1], a, modulus)) != 0)
           goto done;
    }
   
    /* set pos[0] == 1G as expected */
    mpz_set (pos[0]->x, tG->x);
    mpz_set (pos[0]->y, tG->y);
    mpz_set (pos[0]->z, tG->z);

    /* neg[i] == -pos[i] */
    for (j = 0; j < PRECOMPUTE_LENGTH_SMALL; ++j) {
        if ((err = ecc_projective_negate_point(pos[j], neg[j])) != 0)
            goto done;
    }

    /* calculate wMNAF */
    wmnaf = wMNAF(k, WINSIZE, &wmnaf_len);
    if (!wmnaf) {
        err = -2;
        goto done;
    }

    /* actual result computation */

    /* set R to 0 */
    mpz_set_ui(R->x, 0);
    mpz_set_ui(R->y, 0);
    mpz_set_ui(R->z, 1);
   
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
    ecc_del_point(tG);
    for (i = 0; i < PRECOMPUTE_LENGTH_SMALL; i++) {
        ecc_del_point(pos[i]);
        ecc_del_point(neg[i]);
    }
    if (wmnaf) free(wmnaf);
    return err;
}


/* 
 * Return array with wMNAF representation of given mpz_t number x.
 * based on OpenSSL version of this function.
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
    
    if (!mpz_sgn(x)) {
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
        
    if (w <= 0 || w > 7) {
        fprintf(stderr, "`signed char` can represent integers with absolute values less than 2^7.\n");
        return NULL;
    }

    /* 2^w, at most 128 */
    bit = 1 << w; 
    /* 2^(w + 1), at most 256 */
    next_bit = bit << 1;
    /* 2^(w + 1) - 1, at most 255 */
    mask = next_bit - 1;

    if (mpz_sgn(x) < 0) {
        sign = -1;
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
