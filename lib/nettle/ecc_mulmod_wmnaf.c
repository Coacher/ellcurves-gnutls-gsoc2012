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
#ifndef PRECOMPUTE_LENGTH
    #define PRECOMPUTE_LENGTH (1 << (WINSIZE - 1))
#endif

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
    if ((err = ecc_projective_add_point_ng(pos[0], G, pos[1], a, modulus)) != 0)
        goto done;

    /* fill in kG for k = 5, 7, ..., (2^w - 1) */
    for (j = 2; j < PRECOMPUTE_LENGTH; ++j) {
        if ((err = ecc_projective_add_point_ng(pos[j-1], pos[0], pos[j], a, modulus)) != 0)
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
    wmnaf = ecc_wMNAF(k, WINSIZE, &wmnaf_len);
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
                if ((err = ecc_projective_add_point_ng(R, pos[( digit / 2)], R, a, modulus)) != 0)
                    goto done;
            } else {
                if ((err = ecc_projective_add_point_ng(R, neg[(-digit / 2)], R, a, modulus)) != 0)
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
