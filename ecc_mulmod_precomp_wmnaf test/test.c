#include <stdlib.h> 
#include <stdio.h> 
#include <time.h>
#include <gmp.h>

#include "ecc.c"
#include "my_ecc.h"

#define NUM_OF_TRIES 10

int main(void) {
    mpz_t k, a, modulus;
    ecc_point *G, *Rclas, *Rcached;
    int map = 1;

    int rand, i;

    clock_t start;
    double cached_time, classic_time;

    char check;

    mpz_inits(k, a, modulus, NULL);

    G = ecc_new_point();
    Rclas = ecc_new_point();
    Rcached = ecc_new_point();

    ecc_wmnaf_cache_init();

    GNUTLS_ECC_CURVE_LOOP (
        printf("Running test sequence for curve %s with id:%i\n", p->name, p->id);

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
            ecc_mulmod_wmnaf_cached(k, Rcached, p->id, a, map);
            cached_time = ((double) (clock() - start)) / CLOCKS_PER_SEC;

            check = (!mpz_cmp(Rcached->x, Rclas->x)) && (!mpz_cmp(Rcached->y, Rclas->y));

            printf("Classic time: %.15f; Cached time: %.15f; Check: %i\n", classic_time, cached_time, check);
        }

        printf("\n");

    );

    ecc_wmnaf_cache_free();
    ecc_del_point(Rcached);
    ecc_del_point(Rclas);
    ecc_del_point(G);

    mpz_clears(k, a, modulus, NULL);

    return 0;
}
