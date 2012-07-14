#include <stdlib.h> 
#include <stdio.h> 
#include <time.h>
#include <gmp.h>

#include "ecc.c"
#include "my_ecc.h"

#define NUM_OF_TRIES 5

int main(void) {
    mpz_t k, a, modulus;
    ecc_point *G, *R1, *R2, *Raddclas, *Raddng;
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
    Raddng = ecc_new_point();
    
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
            
            start = clock();
            ecc_projective_add_point(R2, R1, Raddclas, a, modulus);
            classic_time = ((double) (clock() - start)) / CLOCKS_PER_SEC;

            ecc_map(Raddclas, modulus);

            start = clock();
            ecc_projective_add_point_ng(R2, R1, Raddng, a, modulus);
            ng_time = ((double) (clock() - start)) / CLOCKS_PER_SEC;

            ecc_map(Raddng, modulus);

            check = (!mpz_cmp(Raddclas->x, Raddng->x)) && (!mpz_cmp(Raddclas->y, Raddng->y));

            printf("Check: %i; Classic time: %.15f; NG time: %.15f\n", check, classic_time, ng_time);
        }

        printf("\n");

    );

    ecc_del_point(Raddng);
    ecc_del_point(Raddclas);
    ecc_del_point(R2);
    ecc_del_point(R1);
    ecc_del_point(G);
    
    mpz_clears(k, a, modulus, NULL);

    return 0;
}
