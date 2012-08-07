#ifndef _MY_ECC_H
#define _MY_ECC_H

/** A point on a ECC curve, stored in Jacbobian format such that (x,y,z) => (x/z^2, y/z^3, 1) when interpretted as affine */
typedef struct {
    /** The x co-ordinate */
    mpz_t x;

    /** The y co-ordinate */
    mpz_t y;

    /** The z co-ordinate */
    mpz_t z;
} ecc_point;

int mp_init_multi(mpz_t *a, ...);
void mp_clear_multi(mpz_t *a, ...);

ecc_point *ecc_new_point (void);
void ecc_del_point (ecc_point * p);
#endif
