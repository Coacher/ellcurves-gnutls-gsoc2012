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

/* Based on public domain code of LibTomCrypt by Tom St Denis.
 * Adapted to gmp and nettle by Nikos Mavrogiannopoulos.
 */

#include "ecc.h"

/* size of sliding window, don't change this! */
#ifdef WITH_WMNAF
    #define WINSIZE 3
#else
    #define WINSIZE 4
#endif

/* length of (half)array of precomputed values for wMNAF */
#define PRECOMPUTE_LENGTH_SMALL (1 << (WINSIZE - 1))

/*
   Perform a point multiplication 
   @param k    The scalar to multiply by
   @param G    The base point
   @param R    [out] Destination for kG
   @param modulus  The modulus of the field the ECC curve is in
   @param map      Boolean whether to map back to affine or not (1==map, 0 == leave in projective)
   @return CRYPT_OK on success
*/
#ifndef WITH_WMNAF
int
ecc_mulmod (mpz_t k, ecc_point * G, ecc_point * R, mpz_t a, mpz_t modulus,
                int map)
#else
int
ecc_mulmod_classic (mpz_t k, ecc_point * G, ecc_point * R, mpz_t a, mpz_t modulus,
                int map)
#endif    
{
   ecc_point *tG, *M[8];
   int        i, j, err, bitidx;
   int        first, bitbuf, bitcpy, mode;

   if (k == NULL || G == NULL || R == NULL || modulus == NULL)
     return -1;

  /* alloc ram for window temps */
  for (i = 0; i < 8; i++) {
      M[i] = ecc_new_point();
      if (M[i] == NULL) {
         for (j = 0; j < i; j++) {
             ecc_del_point(M[j]);
         }

         return -1;
      }
  }

   /* make a copy of G incase R==G */
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

   /* calc the M tab, which holds kG for k==8..15 */
   /* M[0] == 8G */
   if ((err = ecc_projective_dbl_point (tG, M[0], a, modulus)) != 0)
     goto done;

   if ((err = ecc_projective_dbl_point (M[0], M[0], a, modulus)) != 0)
     goto done;

   if ((err = ecc_projective_dbl_point (M[0], M[0], a, modulus)) != 0)
     goto done;
 
   /* now find (8+k)G for k=1..7 */
   for (j = 9; j < 16; j++) {
     if (ecc_projective_add_point(M[j-9], tG, M[j-8], a, modulus) != 0)
       goto done;
   }

   /* setup sliding window */
   mode   = 0;
   bitidx = mpz_size (k) * GMP_NUMB_BITS - 1;
   bitcpy = bitbuf = 0;
   first  = 1;

   /* perform ops */
   for (;;) {
     /* grab next digit as required */
     if (bitidx == -1) {
       break;
     }

     /* grab the next msb from the ltiplicand */
     i = mpz_tstbit (k, bitidx--);

     /* skip leading zero bits */
     if (mode == 0 && i == 0) {
        continue;
     }

     /* if the bit is zero and mode == 1 then we double */
     if (mode == 1 && i == 0) {
        if ((err = ecc_projective_dbl_point(R, R, a, modulus)) != 0)
          goto done;
        continue;
     }

     /* else we add it to the window */
     bitbuf |= (i << (WINSIZE - ++bitcpy));
     mode = 2;

     if (bitcpy == WINSIZE) {
       /* if this is the first window we do a simple copy */
       if (first == 1) {
          /* R = kG [k = first window] */
          mpz_set(R->x, M[bitbuf-8]->x);
          mpz_set(R->y, M[bitbuf-8]->y);
          mpz_set(R->z, M[bitbuf-8]->z);
          first = 0;
       } else {
         /* normal window */
         /* ok window is filled so double as required and add  */
         /* double first */
         for (j = 0; j < WINSIZE; j++) {
           if ((err = ecc_projective_dbl_point(R, R, a, modulus)) != 0)
             goto done;
         }

         /* then add, bitbuf will be 8..15 [8..2^WINSIZE] guaranteed */
         if ((err = ecc_projective_add_point(R, M[bitbuf-8], R, a, modulus)) != 0)
           goto done;
       }
       /* empty window and reset */
       bitcpy = bitbuf = 0;
       mode = 1;
    }
  }

   /* if bits remain then double/add */
   if (mode == 2 && bitcpy > 0) {
     /* double then add */
     for (j = 0; j < bitcpy; j++) {
       /* only double if we have had at least one add first */
       if (first == 0) {
          if ((err = ecc_projective_dbl_point(R, R, a, modulus)) != 0)
            goto done;
       }

       bitbuf <<= 1;
       if ((bitbuf & (1 << WINSIZE)) != 0) {
         if (first == 1){
            /* first add, so copy */
            mpz_set(R->x, tG->x);
            mpz_set(R->y, tG->y);
            mpz_set(R->z, tG->z);
            first = 0;
         } else {
            /* then add */
            if ((err = ecc_projective_add_point(R, tG, R, a, modulus)) != 0)
              goto done;
         }
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
   for (i = 0; i < 8; i++) {
       ecc_del_point(M[i]);
   }
   return err;
}

/* return array containing wMNAF representation of given mpz_t number.
 * inspired by OpenSSL version of this function.
 * overview of this algorithm can be found in
 * Bodo Moller, Improved Techniques for Fast Exponentiation. Information Security and Cryptology – ICISC 2002, Springer-Verlag LNCS 2587, pp. 298–312
 */
static signed char* wMNAF(mpz_t x, int w, size_t *ret_len) {
    int window_val;
    signed char *ret = NULL;
    int sign = 1;
    int bit, next_bit, mask;
    size_t len = 0, j;
    
    if (!mpz_sgn(x)) { 
        ret = malloc(1);
        if (!ret) {
            fprintf(stderr, "Unable to allocate memory in wMNAF.\n");
            exit(EXIT_FAILURE);
        }

        ret[0] = 0;
        *ret_len = 1;
        return ret;
    }
        
    if (w <= 0 || w > 7) {
        fprintf(stderr, "`signed char` can represent integers with absolute values less than 2^7.\n");
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
        fprintf(stderr, "Unable to allocate memory in wMNAF.\n");
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
                fprintf(stderr, "Unexpected digit value: \"%i\".\n", digit);
                exit(EXIT_FAILURE);
            }

            window_val -= digit;

            /* now window_val is 0 or 2^(w+1) in standard wNAF generation;
             * for modified window NAFs, it may also be 2^w
             */
            if (window_val != 0 && window_val != next_bit && window_val != bit) {
                fprinf(stderr, "Unexpected window_val value: \"%i\".\n", window_val);
                exit(EXIT_FAILURE);
            }
        }

        ret[j++] = sign * digit;

        window_val >>= 1;
        window_val += bit * mpz_tstbit(x, j + w);

        if (window_val > next_bit) {
            fprintf(stderr, "Unexpected window_val and next_bit values: \"%i\", \"%i\".\n", window_val, next_bit);
            exit(EXIT_FAILURE);
        }
    }

    if (j > len + 1) {
        fprintf(stderr, "Unexpected j and len + 1 values: \"%i\", \"%i\".\n", j, len+1);
        exit(EXIT_FAILURE);
    }
    len = j;
    
    *ret_len = len;
    return ret;
}

#ifdef WITH_WMNAF
int
ecc_mulmod (mpz_t k, ecc_point * G, ecc_point * R, mpz_t a, mpz_t modulus,
                int map)
#else
int
ecc_mulmod_wMNAF (mpz_t k, ecc_point * G, ecc_point * R, mpz_t a, mpz_t modulus,
                int map)
#endif
{
   ecc_point *tG, *pos[PRECOMPUTE_LENGTH_SMALL], *neg[PRECOMPUTE_LENGTH_SMALL];
   int        i, j, err, bitidx;
   int        first, bitbuf, bitcpy, mode;

   signed char* wmnaf;
   size_t wmnaf_len;
   char digit;

   if (k == NULL || G == NULL || R == NULL || modulus == NULL)
     return -1;

  /* alloc ram for window temps */
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

   /* make a copy of G incase R==G */
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

   /* calc the pos and neg tabs.
    * pos holds kG for k ==  1, 3, ..., (2^w - 1)
    * neg holds kG for k == -1,-3, ...,-(2^w - 1) */
   /* pos[0] == 2G for a while, later it will be set to expected 1G */
   if ((err = ecc_projective_dbl_point (tG, pos[0], a, modulus)) != 0)
     goto done;
   
   /* pos[1] == 3G */
   if (ecc_projective_add_point(pos[1], tG, pos[0], a, modulus) != 0)
     goto done;

   /* now find kG for k = 5,7, ..., (2^w - 1) */
   for (j = 2; j < PRECOMPUTE_LENGTH_SMALL; ++j) {
     if (ecc_projective_add_point(pos[j], pos[0], pos[j-1], a, modulus) != 0)
       goto done;
   }
   
   /* pos[0] == 1G as expected */
   mpz_set (pos[0]->x, tG->x);
   mpz_set (pos[0]->y, tG->y);
   mpz_set (pos[0]->z, tG->z);

   /* neg[i] == -pos[i] */
   for (j = 0; j < PRECOMPUTE_LENGTH_SMALL; ++j) {
       mpz_set (neg[i]->x, pos[i]->x);
       mpz_neg (neg[i]->y, pos[i]->y);
       mpz_set (neg[i]->z, pos[i]->z);
   }

   /* calculate wMNAF */
   wmnaf = wMNAF(k, WINDOWSIZE, &wmnaf_len);

   /* set R to 0 */
   mpz_set_ui(R->x, 0);
   mpz_set_ui(R->y, 0);
   mpz_set_ui(R->z, 1);
   
   /* perform ops */
   for (j = 0; j < wmnaf_len; ++j) {
     if ((err = ecc_projective_dbl_point(R, R, a, modulus)) != 0)
       goto done;

     digit = wmnaf[j];
     if (digit) {
         if (digit > 0) {
            if ((err = ecc_projective_add_point(R, pos[(digit / 2) + (digit % 2) - 1], R, a, modulus)) != 0)
                goto done;
         } else {
             digit = -digit;
            if ((err = ecc_projective_add_point(R, neg[(digit / 2) + (digit % 2) - 1], R, a, modulus)) != 0)
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
   return err;
}

