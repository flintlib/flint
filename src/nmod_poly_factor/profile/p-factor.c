/*
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>

#include "flint.h"
#include "nmod_poly.h"

#define NP 100  /* number of moduli */
#define ND 8   /* number of degrees */

/*
    Benchmarking code for factorisation in nmod_poly.

    Test how the relation between n (degree of polynomial) and p
    affects working time for Cantor-Zassenhaus, Berlekamp and
    Kaltofen-Shoup algorithms. p and n are chosen independently.
*/

int main(void)
{
    FLINT_TEST_INIT(state);
    nmod_poly_t f, g;
    nmod_poly_factor_t res;
    mp_limb_t modulus;
    int i, j, k, n, num;
    double t, T1, T2, T3, T4;

    const slong degs[] = {8, 16, 32, 64, 128, 256, 512, 1024};
    const int iter_count[] = {10000, 5000, 1000, 500, 300, 100, 50, 20};

    

    flint_printf("Random polynomials\n");
    for (i = 0; i < NP; i++)
    {
        modulus = n_randtest_prime(state, 0);
        flint_printf("========== p: %wu\n", modulus);
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = degs[j];
            flint_printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                nmod_poly_init(f, modulus);
                nmod_poly_randtest_not_zero(f, state, n);

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_cantor_zassenhaus(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_berlekamp(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_kaltofen_shoup(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                nmod_poly_clear(f);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;

        }
    }

    /* This code checks whether nmod_poly_factor
       made a correct choice between CZ, B and KS */

    flint_printf("Check choice correctness\n");
    for (i = 0; i < NP; i++)
    {
        modulus = n_randtest_prime(state, 0);
        flint_printf("========== p: %wu\n", modulus);
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = degs[j];
            flint_printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            T4 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                nmod_poly_init(f, modulus);
                nmod_poly_randtest_not_zero(f, state, n);

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_cantor_zassenhaus(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_berlekamp(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_kaltofen_shoup(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T4 += t;

                nmod_poly_clear(f);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf F: %.2lf\n", T1, T2, T3, T4);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;

        }
    }

    flint_printf("Irreducible polynomials\n");
    for (i = 0; i < NP; i++)
    {
        modulus = n_randtest_prime(state, 0);
        flint_printf("========== p: %wu\n", modulus);
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = degs[j];
            flint_printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                nmod_poly_init(f, modulus);
                nmod_poly_randtest_irreducible(f, state, n);

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_cantor_zassenhaus(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_berlekamp(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_kaltofen_shoup(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                nmod_poly_clear(f);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }
    }

    flint_printf("Product of two irreducible polynomials\n");
    for (i = 0; i < NP; i++)
    {
        modulus = n_randtest_prime(state, 0);
        flint_printf("========== p: %wu\n", modulus);
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = (degs[j] >> 1);
            flint_printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                nmod_poly_init(f, modulus);
                nmod_poly_init(g, modulus);
                nmod_poly_randtest_irreducible(f, state, n);
                nmod_poly_randtest_irreducible(g, state, n);
                nmod_poly_mul(f, f, g);

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_cantor_zassenhaus(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_berlekamp(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_kaltofen_shoup(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                nmod_poly_clear(f);
                nmod_poly_clear(g);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }
    }

    flint_printf("Product of 8 small irreducible polynomials\n");
    for (i = 0; i < NP; i++)
    {
        modulus = n_randtest_prime(state, 0);
        flint_printf("========== p: %wu\n", modulus);
        fflush(stdout);

        for (j = 1; j < ND; j++)
        {
            n = (degs[j] >> 3);
            flint_printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                nmod_poly_init(f, modulus);
                nmod_poly_init(g, modulus);
                nmod_poly_randtest_irreducible(f, state, n);
                for (num = 1; num < 8; num++)
                {
                    nmod_poly_randtest_irreducible(g, state, n);
                    nmod_poly_mul(f, f, g);
                }

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_cantor_zassenhaus(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_with_berlekamp(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                nmod_poly_factor_init(res);
                nmod_poly_factor_kaltofen_shoup(res, f);
                nmod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                nmod_poly_clear(f);
                nmod_poly_clear(g);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }
    }

    flint_randclear(state);
    return EXIT_SUCCESS;
}

