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
#include "wmnaf.c"

/* size of sliding window, don't change this! */
#ifndef WINSIZE
    #define WINSIZE 4
#endif

/* length of one array of precomputed values for ecc_mulmod_wmnaf 
 * we have two such arrays for positive and negative multipliers */
#ifndef PRECOMPUTE_LENGTH
    #define PRECOMPUTE_LENGTH (1 << (WINSIZE - 1))
#endif

/* arrays with precomputed values of generator point */
static ecc_point *precomp_pos[PRECOMPUTE_LENGTH], *precomp_neg[PRECOMPUTE_LENGTH];

/* do precompute */
int ecc_precomp_init(ecc_point * G, mpz_t a, mpz_t modulus) {
    int i, j, err;

    if (G == NULL || modulus == NULL)
        return -1;

    /* alloc ram for precomputed values */
    for (i = 0; i < PRECOMPUTE_LENGTH; ++i) {
        precomp_pos[i] = ecc_new_point();
        precomp_neg[i] = ecc_new_point();
        if (precomp_pos[i] == NULL || precomp_neg[i] == NULL) {
            for (j = 0; j < i; ++j) {
              ecc_del_point(precomp_pos[j]);
              ecc_del_point(precomp_neg[j]);
            }

            return -1;
        }
    }

    /* fill in precomp_pos and precomp_neg arrays with precomputed values
     * precomp_pos holds kG for k ==  1, 3, 5, ..., (2^w - 1)
     * precomp_neg holds kG for k == -1,-3,-5, ...,-(2^w - 1)
     */

    /* precomp_pos[0] == 2G for a while, later it will be set to the expected 1G */
    if ((err = ecc_projective_dbl_point(G, precomp_pos[0], a, modulus)) != 0)
        goto done;
   
    /* precomp_pos[1] == 3G */
    if ((err = ecc_projective_add_point_ng(precomp_pos[0], G, precomp_pos[1], a, modulus)) != 0)
        goto done;

    /* fill in kG for k = 5, 7, ..., (2^w - 1) */
    for (j = 2; j < PRECOMPUTE_LENGTH; ++j) {
        if ((err = ecc_projective_add_point_ng(precomp_pos[j-1], precomp_pos[0], precomp_pos[j], a, modulus)) != 0)
           goto done;
    }
   
    /* set precomp_pos[0] == 1G as expected
     * after this step we don't need G at all 
     * and can change it without worries even if R == G */
    mpz_set (precomp_pos[0]->x, G->x);
    mpz_set (precomp_pos[0]->y, G->y);
    mpz_set (precomp_pos[0]->z, G->z);

    /* precomp_neg[i] == -precomp_pos[i] */
    for (j = 0; j < PRECOMPUTE_LENGTH; ++j) {
        if ((err = ecc_projective_negate_point(precomp_pos[j], precomp_neg[j], modulus)) != 0)
            goto done;
    }

    err = 0;

done:
    return err;
}

/* free memory allocated for precomputed values */
void ecc_precomp_free(void) {
    int i;

    for (i = 0; i < PRECOMPUTE_LENGTH; ++i) {
        ecc_del_point(precomp_pos[i]);
        ecc_del_point(precomp_neg[i]);
    }
}

/*
   Perform a point multiplication using wMNAF repr and precomputed values.
   @param k    The scalar to multiply by
   @param R    [out] Destination for kG
   @param modulus  The modulus of the field the ECC curve is in
   @param map      Boolean whether to map back to affine or not (1 == map, 0 == leave in projective)
   @return CRYPT_OK on success
*/
int
ecc_mulmod_precomp_wmnaf (mpz_t k, ecc_point * R, mpz_t a, mpz_t modulus, int map)
{
    int j, err;

    signed char* wmnaf = NULL;
    size_t wmnaf_len;
    signed char digit;

    if (k == NULL || R == NULL || modulus == NULL)
        return -1;

    /* calculate wMNAF */
    wmnaf = wMNAF(k, WINSIZE, &wmnaf_len);
    if (!wmnaf) {
        err = -2;
        goto done;
    }

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
                if ((err = ecc_projective_add_point_ng(R, precomp_pos[( digit / 2)], R, a, modulus)) != 0)
                    goto done;
            } else {
                if ((err = ecc_projective_add_point_ng(R, precomp_neg[(-digit / 2)], R, a, modulus)) != 0)
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
    if (wmnaf) free(wmnaf);
    return err;
}
