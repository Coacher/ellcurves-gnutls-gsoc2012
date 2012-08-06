#include <stdlib.h> 
#include <stdio.h> 
#include <time.h>
#include <gmp.h>

#include "ecc.c"
#include "my_ecc.h"

#define NUM_OF_TRIES 10

#define TEST_ROUTINE(RESULT, FUNC)                                \
do {                                                              \
    start = clock();                                              \
    ecc_projective_add_point(R1, R2, Raddclas, a, modulus);       \
    classic_time = ((double) (clock() - start)) / CLOCKS_PER_SEC; \
    ecc_map(Raddclas, modulus);                                   \
    start = clock();                                              \
    FUNC(R1, R2, RESULT, a, modulus);                             \
    ng_time = ((double) (clock() - start)) / CLOCKS_PER_SEC;      \
    ecc_map(RESULT, modulus);                                     \
    check = (!mpz_cmp(Raddclas->x, RESULT->x)) && (!mpz_cmp(Raddclas->y, RESULT->y));             \
    printf("    Classic time: %.15f; NG time: %.15f; Check: %i\n", classic_time, ng_time, check); \
} while(0)

int main(void) {
    mpz_t k, a, modulus;
    ecc_point *G, *R1, *R2, *Raddclas;
    ecc_point *Raddng_gen, *Raddng_Zs_equal, *Raddng_Zs_equal_one, *Raddng_Z2_is_one;
    int map = 1;

    int rand1, rand2, i;

    clock_t start;
    double ng_time, classic_time;

    char check;

    mpz_inits(k, a, modulus, NULL);

    G = ecc_new_point();
    R1 = ecc_new_point();
    R2 = ecc_new_point();
    Raddclas = ecc_new_point();
    Raddng_gen = ecc_new_point();
    Raddng_Zs_equal = ecc_new_point();
    Raddng_Zs_equal_one = ecc_new_point();
    Raddng_Z2_is_one = ecc_new_point();
    
    GNUTLS_ECC_CURVE_LOOP (
        printf("Running test sequence for curve %s\n", p->name);

        mpz_set_str(G->x,    p->Gx,     16);
        mpz_set_str(G->y,    p->Gy,     16);
        mpz_set_ui (G->z,    1);

        mpz_set_str(modulus, p->prime,  16);
        mpz_set_si(a, -3);

        for (i = 0; i < NUM_OF_TRIES; ++i) {
            rand1 = random();
            rand2 = random();
            printf("Testing for m = %i, n = %i\n", rand1, rand2);

            mpz_set_ui(k, rand1);
            ecc_mulmod_wmnaf(k, G, R1, a, modulus, map);
            
            mpz_set_ui(k, rand2);
            ecc_mulmod_wmnaf(k, G, R2, a, modulus, map);
/*            
            printf("  Testing special case Z1 == Z2 == 1\n");
            ecc_map(R1, modulus);
            ecc_map(R2, modulus);
            TEST_ROUTINE(Raddng_Zs_equal_one);
*/
            printf("  Testing special case Z1 != 1, Z2 == 1\n");
            mpz_mul_ui(R1->z, R1->z, rand1);
            mpz_mul_ui(R1->x, R1->x, rand1*rand1);
            mpz_mul_ui(R1->y, R1->y, rand1*rand1*rand1);
            TEST_ROUTINE(Raddng_Z2_is_one, ecc_projective_madd);
/*
            printf("  Testing special case Z1 == Z2 != 1\n");
            mpz_mul_ui(R2->z, R2->z, rand1);
            mpz_mul_ui(R2->x, R2->x, rand1*rand1);
            mpz_mul_ui(R2->y, R2->y, rand1*rand1*rand1);
            TEST_ROUTINE(Raddng_Zs_equal);
*/
            printf("  Testing general case\n");
            mpz_mul_ui(R1->z, R1->z, rand2);
            mpz_mul_ui(R1->x, R1->x, rand2*rand2);
            mpz_mul_ui(R1->y, R1->y, rand2*rand2*rand2);
            TEST_ROUTINE(Raddng_gen, ecc_projective_add_point_ng);
        }

        printf("\n");
    );

    ecc_del_point(Raddng_Z2_is_one);
    ecc_del_point(Raddng_Zs_equal_one);
    ecc_del_point(Raddng_Zs_equal);
    ecc_del_point(Raddng_gen);
    ecc_del_point(Raddclas);
    ecc_del_point(R2);
    ecc_del_point(R1);
    ecc_del_point(G);
    
    mpz_clears(k, a, modulus, NULL);

    return 0;
}
