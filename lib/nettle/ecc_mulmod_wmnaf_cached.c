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

/* needed for gnutls_* types */
#include <gnutls_num.h>
#include <x509/x509_int.h>
#include <x509/common.h>

#include "ecc.h"

/* size of sliding window, don't change this! */
#ifndef WINSIZE
    #define WINSIZE 4
#endif

/* length of one array of precomputed values
 * we have two such arrays for positive and negative multipliers */
#ifndef PRECOMPUTE_LENGTH
    #define PRECOMPUTE_LENGTH (1 << (WINSIZE - 1))
#endif

/* per-curve cache structure */
struct gnutls_ecc_curve_cache_entry_t {
    gnutls_ecc_curve_t id;

    /** The prime that defines the field the curve is in */
    mpz_t modulus;

    /** The array of positive multipliers of G */
    ecc_point *pos[PRECOMPUTE_LENGTH];

    /** The array of positive multipliers of G */
    ecc_point *neg[PRECOMPUTE_LENGTH];
};
typedef struct gnutls_ecc_curve_cache_entry_t gnutls_ecc_curve_cache_entry_t;

/* global cache */
gnutls_ecc_curve_cache_entry_t* ecc_wmnaf_cache = NULL;

/* free single cache entry */
static void ecc_wmnaf_cache_entry_free(gnutls_ecc_curve_cache_entry_t *p) {
    int i;

    mpz_clear(p->modulus);

    for(i = 0; p->pos[i]; ++i) {
        ecc_del_point(p->pos[i]);
        ecc_del_point(p->neg[i]);
    }
}

/* free curves caches */
static void _ecc_wmnaf_cache_free(gnutls_ecc_curve_cache_entry_t *p) {
    if (p) {
        for (; p->id; ++p) {
            ecc_wmnaf_cache_entry_free(p);
        }
        free(p);
    }
}

/* free curves caches */
void ecc_wmnaf_cache_free(void) {
    _ecc_wmnaf_cache_free(ecc_wmnaf_cache);
}

/* initialize curves caches */
static int _ecc_wmnaf_cache_init(gnutls_ecc_curve_cache_entry_t **cache) {
    int i, j, k, err;
    ecc_point* G;
    mpz_t a;

    gnutls_ecc_curve_cache_entry_t* ret;

    const gnutls_ecc_curve_t *p;
    const gnutls_ecc_curve_entry_st *st;

    ret = (gnutls_ecc_curve_cache_entry_t*) malloc(MAX_ALGOS*sizeof(gnutls_ecc_curve_cache_entry_t));
    if (!ret) return 1;

    mpz_init(a);
    G = ecc_new_point();

    /* get supported curves' ids */
    p = gnutls_ecc_curve_list();

    mpz_set_si(a, -3);

    for (j = 0; *p; ++p, ++j) {
        st =  _gnutls_ecc_curve_get_params(*p);

        /* set id */
        ret[j].id = *p;

        /* set modulus*/
        mpz_init(ret[j].modulus);
        mpz_set_str(ret[j].modulus, st->prime, 16);

        /* get generator point*/
        mpz_set_str(G->x, st->Gx, 16);
        mpz_set_str(G->y, st->Gy, 16);
        mpz_set_ui (G->z, 1);        

        /* alloc ram for precomputed values */
        for (i = 0; i < PRECOMPUTE_LENGTH; ++i) {
            ret[j].pos[i] = ecc_new_point();
            ret[j].neg[i] = ecc_new_point();
            if (ret[j].pos[i] == NULL || ret[j].neg[i] == NULL) {
                for (k = 0; k < i; ++k) {
                  ecc_del_point(ret[j].pos[k]);
                  ecc_del_point(ret[j].neg[k]);
                }

                err = -11;
                goto done;
            }
        }

        /* fill in pos and neg arrays with precomputed values
         * pos holds kG for k ==  1, 3, 5, ..., (2^w - 1)
         * neg holds kG for k == -1,-3,-5, ...,-(2^w - 1)
         */

        /* pos[0] == 2G for a while, later it will be set to the expected 1G */
        if ((err = ecc_projective_dbl_point(G, ret[j].pos[0], a, ret[j].modulus)) != 0)
            goto done;

        /* pos[1] == 3G */
        if ((err = ecc_projective_add_point_ng(ret[j].pos[0], G, ret[j].pos[1], a, ret[j].modulus)) != 0)
            goto done;

        /* fill in kG for k = 5, 7, ..., (2^w - 1) */
        for (k = 2; k < PRECOMPUTE_LENGTH; ++k) {
            if ((err = ecc_projective_add_point_ng(ret[j].pos[k-1], ret[j].pos[0], ret[j].pos[k], a, ret[j].modulus)) != 0)
               goto done;
        }

        /* set pos[0] == 1G as expected
         * after this step we don't need G at all */
        mpz_set (ret[j].pos[0]->x, G->x);
        mpz_set (ret[j].pos[0]->y, G->y);
        mpz_set (ret[j].pos[0]->z, G->z);

        /* map to affine all elements in pos
         * this will allow to use ecc_projective_madd later
         * set neg[i] == -pos[i] */
        for (k = 0; k < PRECOMPUTE_LENGTH; ++k) {
            if ((err = ecc_map(ret[j].pos[k], ret[j].modulus)) != 0)
                goto done;

            if ((err = ecc_projective_negate_point(ret[j].pos[k], ret[j].neg[k], ret[j].modulus)) != 0)
                goto done;
        }
    }

    ret[++j].id = 0;

    err = 0;

    *cache = ret;
    goto done;
done:
    mpz_clear(a);
    ecc_del_point(G);
    if (err) {
        if (err == -11) {
            mpz_clear(ret[j].modulus);
        }

        for(k = 0; k < j; ++k) {
            ecc_wmnaf_cache_entry_free(&ret[k]);
        }

        free(ret);
        *cache = NULL;
    }
    return err;
}

/* initialize curves caches */
int ecc_wmnaf_cache_init(void) {
    return _ecc_wmnaf_cache_init(&ecc_wmnaf_cache);
}

/* perform cache lookup i.e. return curve's cache by its id */
static gnutls_ecc_curve_cache_entry_t * ecc_wmnaf_cache_lookup(gnutls_ecc_curve_t id) {
    unsigned int i, pid;
    static gnutls_ecc_curve_cache_entry_t * ret = NULL;

    if (!id) return NULL;

    if (ret->id == id) return ret;

    for(i = 0; (pid = ecc_wmnaf_cache[i].id); ++i) {
        if (pid == id) return &ecc_wmnaf_cache[i];
    }

    return NULL;
}

/*
   Perform a point multiplication utilizing cache
   @param k    The scalar to multiply by
   @param R    [out] Destination for kG
   @param id   Curve's id
   @return CRYPT_OK on success
*/
int
ecc_mulmod_wmnaf_cached (mpz_t k, ecc_point * R, gnutls_ecc_curve_t id, mpz_t a, int map)
{
    int j, err;
    
    gnutls_ecc_curve_cache_entry_t * cache = NULL;
    signed char* wmnaf = NULL;
    size_t wmnaf_len;
    signed char digit;

    if (k == NULL || R == NULL || id == 0)
        return -1;

    /* calculate wMNAF */
    wmnaf = ecc_wMNAF(k, WINSIZE, &wmnaf_len);
    if (!wmnaf) {
        err = -2;
        goto done;
    }

    /* set R to neutral */
    mpz_set_ui(R->x, 1);
    mpz_set_ui(R->y, 1);
    mpz_set_ui(R->z, 0);

    /* do cache lookup */
    cache = ecc_wmnaf_cache_lookup(id);
    if (!cache) {
        err = -1;
        goto done;
    }

    /* perform ops */
    for (j = wmnaf_len - 1; j >= 0; --j) {
        if ((err = ecc_projective_dbl_point(R, R, a, cache->modulus)) != 0)
            goto done;

        digit = wmnaf[j];

        if (digit) {
            if (digit > 0) {
                if ((err = ecc_projective_madd(R, cache->pos[( digit / 2)], R, a, cache->modulus)) != 0)
                    goto done;
            } else {
                if ((err = ecc_projective_madd(R, cache->neg[(-digit / 2)], R, a, cache->modulus)) != 0)
                    goto done;
            }
        }
    }


    /* map R back from projective space */
    if (map) {
        err = ecc_map(R, cache->modulus);
    } else {
        err = 0;
    }
done:
    if (wmnaf) free(wmnaf);
    return err;
}
