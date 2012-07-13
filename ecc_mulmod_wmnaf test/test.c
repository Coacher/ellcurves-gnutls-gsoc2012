#include <stdlib.h> 
#include <stdio.h> 
#include <time.h>
#include <gmp.h>

#include "ecc.c"
#include "my_ecc.h"

#define NUM_OF_TRIES 5

int main(void) {
    mpz_t k, a, modulus;
    ecc_point *G, *Rclas, *Rwmnaf;
    int map = 1;
    
    int rand, i;

    clock_t start;
    double wmnaf_time, classic_time;

    char check;

    mpz_inits(k, a, modulus, NULL);

    G = ecc_new_point();
    Rclas = ecc_new_point();
    Rwmnaf = ecc_new_point();
    
    GNUTLS_ECC_CURVE_LOOP (
        printf("Running test sequence for curve %s\n", p->name);

        mpz_set_str(G->x,    p->Gx,     16);
        mpz_set_str(G->y,    p->Gy,     16);
        mpz_set_ui (G->z,    1);

        mpz_set_str(modulus, p->prime,  16);
        mpz_set_si(a, -3);

        for (i = 0; i < NUM_OF_TRIES; ++i) {
            rand = random();
            printf("Testing for k = %i\n", rand);

            mpz_set_ui(k, rand);

            start = clock();
            ecc_mulmod(k, G, Rclas, a, modulus, map);
            classic_time = ((double) (clock() - start)) / CLOCKS_PER_SEC;

            start = clock();
            ecc_mulmod_wmnaf(k, G, Rwmnaf, a, modulus, map);
            wmnaf_time = ((double) (clock() - start)) / CLOCKS_PER_SEC;

            check = (!mpz_cmp(Rwmnaf->x, Rclas->x)) && (!mpz_cmp(Rwmnaf->y, Rclas->y));

            printf("Check: %i; Classic time: %.15f; wMNAF time: %.15f\n", check, classic_time, wmnaf_time);
        }

        printf("\n");

    );

    return 0;
}
