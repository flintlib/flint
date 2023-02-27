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
#include "fmpz_mod_poly.h"

#define NP 20  /* number of moduli */
#define ND 8   /* number of degrees */

/*
    Benchmarking code for factorisation in fmpz_mod_poly.

    Test how the relation between n (degree of polynomial) and p
    affects working time for Cantor-Zassenhaus, Berlekamp and
    Kaltofen-Shoup algorithms. p and n are chosen independently.
*/

int main(void)
{
    FLINT_TEST_INIT(state);
    fmpz_mod_poly_t f, g;
    fmpz_mod_poly_factor_t res;
    fmpz_t p;
    fmpz_mod_ctx_t ctx;
    mpz_t pz, curr;
    int i, j, k, n, num;
    double t, T1, T2, T3;
    const slong degs[] = {8, 16, 32, 64, 128, 256, 512, 1024};
    const int iter_count[] = {10000, 5000, 1000, 500, 300, 100, 50, 20};

    mpz_init(pz);
    mpz_init(curr);
    fmpz_init(p);
    fmpz_mod_ctx_init_ui(ctx, 2);

    flint_printf("Random polynomials\n");
    flint_mpz_set_ui(pz, 2);
    flint_mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        fmpz_mod_ctx_set_modulus(ctx, p);
        flint_printf("========== p: "); fmpz_print(p); flint_printf(" ==========\n");
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
                fmpz_mod_poly_init(f, ctx);
                fmpz_mod_poly_randtest_not_zero(f, state, n, ctx);

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_berlekamp(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f, ctx);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        flint_mpz_mul_ui(curr, curr, 10);
    }

    /* This code checks whether fmpz_mod_poly_factor
       made a correct choice between CZ and KS */

    flint_printf("Check choice correctness\n");
    flint_mpz_set_ui(pz, 2);
    flint_mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        fmpz_mod_ctx_set_modulus(ctx, p);
        flint_printf("========== p: "); fmpz_print(p); flint_printf(" ==========\n");
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
                fmpz_mod_poly_init(f, ctx);
                fmpz_mod_poly_randtest_not_zero(f, state, n, ctx);

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f, ctx);
            }

            flint_printf("CZ: %.2lf F: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        flint_mpz_mul_ui(curr, curr, 10);
    }

    flint_printf("Irreducible polynomials\n");
    flint_mpz_set_ui(pz, 2);
    flint_mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        fmpz_mod_ctx_set_modulus(ctx, p);
        flint_printf("========== p: "); fmpz_print(p); flint_printf(" ==========\n");
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
                fmpz_mod_poly_init(f, ctx);
                fmpz_mod_poly_randtest_irreducible(f, state, n, ctx);

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_berlekamp(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f, ctx);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        flint_mpz_mul_ui(curr, curr, 10);
    }

    flint_printf("Product of two irreducible polynomials\n");
    flint_mpz_set_ui(pz, 2);
    flint_mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        fmpz_mod_ctx_set_modulus(ctx, p);
        flint_printf("========== p: "); fmpz_print(p); flint_printf(" ==========\n");
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
                fmpz_mod_poly_init(f, ctx);
                fmpz_mod_poly_init(g, ctx);
                fmpz_mod_poly_randtest_irreducible(f, state, n, ctx);
                fmpz_mod_poly_randtest_irreducible(g, state, n, ctx);
                fmpz_mod_poly_mul(f, f, g, ctx);

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_berlekamp(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f, ctx);
                fmpz_mod_poly_clear(g, ctx);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        flint_mpz_mul_ui(curr, curr, 10);
    }

    flint_printf("Product of 8 small irreducible polynomials\n");
    flint_mpz_set_ui(pz, 2);
    flint_mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        fmpz_mod_ctx_set_modulus(ctx, p);
        flint_printf("========== p: "); fmpz_print(p); flint_printf(" ==========\n");
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
                fmpz_mod_poly_init(f, ctx);
                fmpz_mod_poly_init(g, ctx);
                fmpz_mod_poly_randtest_irreducible(f, state, n, ctx);
                for (num = 1; num < 8; num++)
                {
                    fmpz_mod_poly_randtest_irreducible(g, state, n, ctx);
                    fmpz_mod_poly_mul(f, f, g, ctx);
                }

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_berlekamp(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res, ctx);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f, ctx);
                fmpz_mod_poly_factor_clear(res, ctx);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f, ctx);
                fmpz_mod_poly_clear(g, ctx);
            }

            flint_printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        flint_mpz_mul_ui(curr, curr, 10);
    }

    mpz_clear(pz);
    mpz_clear(curr);
    fmpz_clear(p);
    fmpz_mod_ctx_clear(ctx);
    FLINT_TEST_CLEANUP(state);
    
    return EXIT_SUCCESS;
}

