#include <stdlib.h> 
#include <stdio.h> 
#include <gmp.h>

#include "ecc.h"
#include "ecc_point.h"

#define NUM_OF_TRIES_PER_CURVE 10

#define MAP 1

#define TEST_ROUTINE(RESULT, FUNC)                                \
do {                                                              \
    ecc_projective_add_point(R1, R2, Rclassic, a, modulus);       \
    ecc_map(Rclassic, modulus);                                   \
    FUNC(R1, R2, RESULT, a, modulus);                             \
    ecc_map(RESULT, modulus);                                     \
    check = (!mpz_cmp(Rclassic->x, RESULT->x)) && (!mpz_cmp(Rclassic->y, RESULT->y)); \
    printf("check: %i\n", check);                                 \
} while(0)

int main(void) {

    mpz_t k, a, modulus;
    ecc_point *G, *R1, *R2, *Rclassic;
    ecc_point *Raddng, *Rmadd;

    int rand1, rand2, i;
    char check;

    mpz_inits(k, a, modulus, NULL);

    G          = ecc_new_point();
    R1         = ecc_new_point();
    R2         = ecc_new_point();
    Rclassic   = ecc_new_point();
    Raddng     = ecc_new_point();
    Rmadd      = ecc_new_point();

    printf("\nRunning tests for ecc_projective_add_point_ng() and ecc_projective_madd()\n\n");

    GNUTLS_ECC_CURVE_LOOP (
        printf("Running test sequence for curve %s with id:%i\n", p->name, p->id);

        mpz_set_str(G->x,    p->Gx,     16);
        mpz_set_str(G->y,    p->Gy,     16);
        mpz_set_ui (G->z,    1);

        mpz_set_str(modulus, p->prime,  16);
        mpz_set_str(a,       p->A,      16);

        for (i = 0; i < NUM_OF_TRIES_PER_CURVE; ++i) {
            rand1 = random();
            rand2 = random();
            printf("Testing r1 = %-i, r2 = %-i:\n", rand1, rand2);

            mpz_set_ui(k, rand1);
            ecc_mulmod_wmnaf(k, G, R1, a, modulus, MAP);

            mpz_set_ui(k, rand2);
            ecc_mulmod_wmnaf(k, G, R2, a, modulus, MAP);

            printf("\tTesting madd\t\t");
            /* make Z1 != 1 */
            mpz_mul_ui(R1->z, R1->z, rand1);
            mpz_mul_ui(R1->x, R1->x, rand1*rand1);
            mpz_mul_ui(R1->y, R1->y, rand1*rand1*rand1);
            TEST_ROUTINE(Rmadd, ecc_projective_madd);

            printf("\tTesting add_point_ng\t");
            /* make Z2 != 1 */
            mpz_mul_ui(R2->z, R2->z, rand2);
            mpz_mul_ui(R2->x, R2->x, rand2*rand2);
            mpz_mul_ui(R2->y, R2->y, rand2*rand2*rand2);
            TEST_ROUTINE(Raddng, ecc_projective_add_point_ng);
        }
        printf("\n");
    );
    printf("\nDoing sanity checks for ecc_projective_add_point_ng() and ecc_projective_madd()\n");
    printf("Nothing will be printed out, but the test utility should't crash.\n");

    mpz_set_ui(R1->x,    1);
    mpz_set_ui(R1->y,    1);
    mpz_set_ui(R1->z,    1);
    ecc_projective_add_point_ng (G, R1, R2,  a, modulus);
    ecc_projective_madd         (G, R1, R2,  a, modulus);

    mpz_set_ui(R1->x,    1);
    mpz_set_ui(R1->y,    1);
    mpz_set_ui(R1->z,    0);
    ecc_projective_add_point_ng (G, R1, R2,  a, modulus);
    ecc_projective_madd         (G, R1, R2,  a, modulus);

    mpz_set_ui(R1->x,    1);
    mpz_set_ui(R1->y,    0);
    mpz_set_ui(R1->z,    0);
    ecc_projective_add_point_ng (G, R1, R2,  a, modulus);
    ecc_projective_madd         (G, R1, R2,  a, modulus);

    mpz_set_ui(R1->x,    0);
    mpz_set_ui(R1->y,    1);
    mpz_set_ui(R1->z,    0);
    ecc_projective_add_point_ng (G, R1, R2,  a, modulus);
    ecc_projective_madd         (G, R1, R2,  a, modulus);

    mpz_set_ui(R1->x,    0);
    mpz_set_ui(R1->y,    0);
    mpz_set_ui(R1->z,    1);
    ecc_projective_add_point_ng (G, R1, R2,  a, modulus);
    ecc_projective_madd         (G, R1, R2,  a, modulus);

    mpz_set_ui(R1->x,    0);
    mpz_set_ui(R1->y,    0);
    mpz_set_ui(R1->z,    0);
    ecc_projective_add_point_ng (G, R1, R2,  a, modulus);
    ecc_projective_madd         (G, R1, R2,  a, modulus);

    mpz_set_ui(R1->x,    100);
    mpz_set_ui(R1->y,    200);
    mpz_set_ui(R1->z,    300);
    ecc_projective_add_point_ng (G, R1, R2,  a, modulus);
    ecc_projective_madd         (G, R1, R2,  a, modulus);

    printf("Sanity checks successfully passed.\n");

    ecc_del_point(Rmadd);
    ecc_del_point(Raddng);
    ecc_del_point(Rclassic);
    ecc_del_point(R2);
    ecc_del_point(R1);
    ecc_del_point(G);

    mpz_clears(k, a, modulus, NULL);

    return 0;
}
