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
    ecc_point *Rwmnaf;

    int rand, i;
    char check;

    mpz_inits(k, a, modulus, NULL);

    G       = ecc_new_point();
    Rclassic= ecc_new_point();
    Rwmnaf  = ecc_new_point();

    printf("\nRunning tests for ecc_mulmod_wmnaf()\n\n");
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

            ecc_mulmod      (k, G, Rclassic, a, modulus, MAP);
            ecc_mulmod_wmnaf(k, G, Rwmnaf,   a, modulus, MAP);

            check = ( (!mpz_cmp(Rwmnaf->x, Rclassic->x)) && (!mpz_cmp(Rwmnaf->y, Rclassic->y)) );

            printf("check: %i\n", check);
        }
        printf("\n");
    );
    printf("\nDoing sanity checks for ecc_mulmod_wmnaf()\n");
    printf("Nothing will be printed out, but the test utility should't crash.\n");

    mpz_set_ui(G->x,    1);
    mpz_set_ui(G->y,    1);
    mpz_set_ui(G->z,    1);
    ecc_mulmod_wmnaf(k, G,     Rwmnaf,  a, modulus, MAP);

    mpz_set_ui(G->x,    1);
    mpz_set_ui(G->y,    1);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf(k, G,     Rwmnaf,  a, modulus, MAP);

    mpz_set_ui(G->x,    1);
    mpz_set_ui(G->y,    0);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf(k, G,     Rwmnaf,  a, modulus, MAP);

    mpz_set_ui(G->x,    0);
    mpz_set_ui(G->y,    1);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf(k, G,     Rwmnaf,  a, modulus, MAP);

    mpz_set_ui(G->x,    0);
    mpz_set_ui(G->y,    0);
    mpz_set_ui(G->z,    1);
    ecc_mulmod_wmnaf(k, G,     Rwmnaf,  a, modulus, MAP);

    mpz_set_ui(G->x,    0);
    mpz_set_ui(G->y,    0);
    mpz_set_ui(G->z,    0);
    ecc_mulmod_wmnaf(k, G,     Rwmnaf,  a, modulus, MAP);

    mpz_set_ui(G->x,    100);
    mpz_set_ui(G->y,    200);
    mpz_set_ui(G->z,    300);
    ecc_mulmod_wmnaf(k, G,     Rwmnaf,  a, modulus, MAP);

    printf("Sanity checks successfully passed.\n");

    ecc_del_point(Rwmnaf);
    ecc_del_point(Rclassic);
    ecc_del_point(G);

    mpz_clears(k, a, modulus, NULL);

    return 0;
}
