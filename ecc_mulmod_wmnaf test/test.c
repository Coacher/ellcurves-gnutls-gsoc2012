#include <stdlib.h> 
#include <stdio.h> 
#include <gmp.h>

#include "my_ecc.h"

int main(void) {
    mpz_t k, a, modulus;
    ecc_point *G, *Rclas, *Rwmnaf;
    int map = 0;

    mpz_inits(k, a, modulus, NULL);

    G = ecc_new_point();
    Rclas = ecc_new_point();
    Rwmnaf = ecc_new_point();

    /* SECP256R1 */
    mpz_set_str(G->x, "0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 0);
    mpz_set_str(G->y, "0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 0);
    mpz_set_str(a,    "0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 0);
    mpz_set_str(modulus, "0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 0);

    mpz_set_ui(k, 100);

    ecc_mulmod_wmnaf(k, G, Rclas, a, modulus, map);

    return 0;
}
