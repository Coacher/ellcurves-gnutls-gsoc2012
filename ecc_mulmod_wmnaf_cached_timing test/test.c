#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>

#include "ecc.h"
#include "ecc_point.h"

/* 540 is near 512 and divisible by 3,4,5
 * and +1 bit for storing trailing zero */
#define BIT_LENGTH 541

/* strings to store multipliers
 * their names represent bits density */
char zero           [BIT_LENGTH];
char one_fifth      [BIT_LENGTH];
char one_fourth     [BIT_LENGTH];
char one_third      [BIT_LENGTH];
char two_fifths     [BIT_LENGTH];
char half           [BIT_LENGTH];
char three_fifths   [BIT_LENGTH];
char two_thirds     [BIT_LENGTH];
char three_fourths  [BIT_LENGTH];
char four_fifths    [BIT_LENGTH];
char one            [BIT_LENGTH];

#define NUM_OF_MULTIPLIERS 11

struct multiplier {
    char* bitstring;
    char* desc;
} multipliers[NUM_OF_MULTIPLIERS] = {
    zero,           "zero\t\t",
    one_fifth,      "one fifth\t",
    one_fourth,     "one fourth\t",
    one_third,      "one third\t",
    two_fifths,     "two fifths\t",
    half,           "half\t\t",
    three_fifths,   "three fifths\t",
    two_thirds,     "two thirds\t",
    three_fourths,  "three fourths\t",
    four_fifths,    "four fifths\t",
    one,            "one\t\t",
};

/* set bitstrings */
static void init_multipliers(void);

#define NUM_OF_TRIES_PER_CURVE 1

#define NUM_OF_TRIES_PER_MULTIPLIER 1000

#define MAP 1

int main(void) {

    mpz_t k, a, modulus;
    ecc_point *G, *Rclassic;
    ecc_point *Rtiming;

    clock_t classic_start, classic_stop;
    clock_t wmnaf_start, wmnaf_stop;
    double wmnaf_time, classic_time;

    int i, j, l;
    char check;

    mpz_inits(k, a, modulus, NULL);

    G       = ecc_new_point();
    Rtiming = ecc_new_point();
    Rclassic= ecc_new_point();

    init_multipliers();

    ecc_wmnaf_cache_init();

    printf("\nRunning timing tests for ecc_mulmod_timing() and ecc_mulmod_wmnaf_cached_timing()\n\n");

    GNUTLS_ECC_CURVE_LOOP (
        printf("Running test sequence for curve %s with id:%i\n", p->name, p->id);

        mpz_set_str(G->x,    p->Gx,     16);
        mpz_set_str(G->y,    p->Gy,     16);
        mpz_set_ui (G->z,    1);

        mpz_set_str(modulus, p->prime,  16);
        mpz_set_str(a,       p->A,      16);

        for (i = 0; i < NUM_OF_TRIES_PER_CURVE; ++i) {
                for (l = 0; l < NUM_OF_MULTIPLIERS; ++l) {
                    printf("Testing %s bit density\t", multipliers[l].desc);

                    mpz_set_str(k, multipliers[l].bitstring, 2);

                    /* run original code */
                    classic_start = clock();
                    for (j = 0; j < NUM_OF_TRIES_PER_MULTIPLIER; ++j) {
                        ecc_mulmod_timing ( k,  G,  Rclassic,   a,  modulus,    MAP);
                    }
                    classic_stop = clock();

                    classic_time = (double) (classic_stop - classic_start) / CLOCKS_PER_SEC;

                    /* run wmnaf code */
                    wmnaf_start = clock();
                    for (j = 0; j < NUM_OF_TRIES_PER_MULTIPLIER; ++j) {
                        ecc_mulmod_wmnaf_cached_timing ( k, p->id,  Rtiming,  a,  modulus,    MAP);
                    }
                    wmnaf_stop = clock();

                    wmnaf_time = (double) (wmnaf_stop - wmnaf_start) / CLOCKS_PER_SEC;


                    check = ( (!mpz_cmp(Rtiming->x, Rclassic->x)) && (!mpz_cmp(Rtiming->y, Rclassic->y)) );

                    printf("wMNAF time: %.3f, classic time: %.3f, check: %i\n", wmnaf_time, classic_time, check);
                }
        }
        printf("\n");
    );

    ecc_wmnaf_cache_free();

    ecc_del_point(Rtiming);
    ecc_del_point(Rclassic);
    ecc_del_point(G);

    mpz_clears(k, a, modulus, NULL);

    return 0;
}

static void init_multipliers(void) {
    int i;

    for (i = 0; i < (BIT_LENGTH - 1); ++i) {

        zero[i] = '0';

        one[i]  = '1';

        /* fifths */
        if ((i % 5) == 1) {
            one_fifth[i] = '1';
        } else {
            one_fifth[i] = '0';
        }

        if ((i % 5) == 2) {
            two_fifths[i] = '1';
        } else {
            two_fifths[i] = '0';
        }

        if ((i % 5) == 3) {
            three_fifths[i] = '1';
        } else {
            three_fifths[i] = '0';
        }

        if ((i % 5) == 4) {
            four_fifths[i] = '1';
        } else {
            four_fifths[i] = '0';
        }

        /* fourths */
        if ((i % 4) == 1) {
            one_fourth[i] = '1';
        } else {
            one_fourth[i] = '0';
        }

        if ((i % 4) == 3) {
            three_fourths[i] = '1';
        } else {
            three_fourths[i] = '0';
        }

        /* thirds */
        if ((i % 3) == 1) {
            one_third[i] = '1';
        } else {
            one_third[i] = '0';
        }

        if ((i % 3) == 2) {
            two_thirds[i] = '1';
        } else {
            two_thirds[i] = '0';
        }

        /* half */
        if ((i % 2) == 0) {
            half[i] = '1';
        } else {
            half[i] = '0';
        }

    }

    zero[0] = '1';

    zero            [BIT_LENGTH - 1] = '\0';
    one_fifth       [BIT_LENGTH - 1] = '\0';
    one_fourth      [BIT_LENGTH - 1] = '\0';
    one_third       [BIT_LENGTH - 1] = '\0';
    two_fifths      [BIT_LENGTH - 1] = '\0';
    half            [BIT_LENGTH - 1] = '\0';
    three_fifths    [BIT_LENGTH - 1] = '\0';
    two_thirds      [BIT_LENGTH - 1] = '\0';
    three_fourths   [BIT_LENGTH - 1] = '\0';
    four_fifths     [BIT_LENGTH - 1] = '\0';
    one             [BIT_LENGTH - 1] = '\0';
}
