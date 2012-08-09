#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "ecc.h"
#include "ecc_point.h"

#define NUM_OF_TRIES_PER_CURVE 10

#define MAP 1

int main(void) {

    mpz_t k, a, modulus;
    ecc_point *G, *Rclassic;
    ecc_point *Rcached, *Rlookup;

    int rand, i;
    char check;

    mpz_inits(k, a, modulus, NULL);

    G       = ecc_new_point();
    Rclassic= ecc_new_point();
    Rcached = ecc_new_point();
    Rlookup = ecc_new_point();

    ecc_wmnaf_cache_init();

    printf("\nRunning tests for ecc_mulmod_wmnaf_cached() and ecc_mulmod_wmnaf_cached_lookup()\n\n");

    GNUTLS_ECC_CURVE_LOOP (
        printf("Running test sequence for curve %s with id:%i\n", p->name, p->id);

        mpz_set_str(G->x,    p->Gx,     16);
        mpz_set_str(G->y,    p->Gy,     16);
        mpz_set_ui (G->z,    1);

        mpz_set_str(modulus, p->prime,  16);
        mpz_set_str(a,       p->A,      16);

        for (i = 0; i < NUM_OF_TRIES_PER_CURVE; ++i) {

            rand = random();
            mpz_set_ui(k, rand);
            printf("Testing k = %-i\t", rand);

            ecc_mulmod                    (k, G,     Rclassic, a, modulus, MAP);
            ecc_mulmod_wmnaf_cached       (k, p->id, Rcached,  a, modulus, MAP);
            ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

            check = ( (!mpz_cmp(Rcached->x, Rclassic->x)) && (!mpz_cmp(Rlookup->x, Rclassic->x)) \
                   && (!mpz_cmp(Rcached->y, Rclassic->y)) && (!mpz_cmp(Rlookup->y, Rclassic->y)) );

            printf("check: %i\n", check);
        }
        printf("\n");
    );
    printf("\nDoing sanity checks for ecc_mulmod_wmnaf_cached() and ecc_mulmod_wmnaf_cached_lookup()\n");
    printf("Nothing will be printed out, but the test utility should't crash.\n");

    /* all other id values like -100 or 1000 will crash the program
     * it is not an issue it is due to gnutls_ecc_curve_t definition
     * as long as you pass id of this type you are safe (and we do) */
    ecc_mulmod_wmnaf_cached       (k, 0,     Rcached,  a, modulus, MAP);

    mpz_set_ui(G->x,    1);
    mpz_set_ui(G->y,    1);
    mpz_set_ui(G->z,    1);
    ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

    mpz_set_ui(G->x,    1);
    mpz_set_ui(G->y,    1);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

    mpz_set_ui(G->x,    1);
    mpz_set_ui(G->y,    0);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

    mpz_set_ui(G->x,    0);
    mpz_set_ui(G->y,    1);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

    mpz_set_ui(G->x,    0);
    mpz_set_ui(G->y,    0);
    mpz_set_ui(G->z,    1);
    ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

    mpz_set_ui(G->x,    0);
    mpz_set_ui(G->y,    0);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

    mpz_set_ui(G->x,    100);
    mpz_set_ui(G->y,    200);
    mpz_set_ui(G->z,    300);
    ecc_mulmod_wmnaf_cached_lookup(k, G,     Rlookup,  a, modulus, MAP);

    printf("Sanity checks successfully passed.\n");

    ecc_wmnaf_cache_free();
    ecc_del_point(Rlookup);
    ecc_del_point(Rcached);
    ecc_del_point(Rclassic);
    ecc_del_point(G);

    mpz_clears(k, a, modulus, NULL);

    return 0;
}
