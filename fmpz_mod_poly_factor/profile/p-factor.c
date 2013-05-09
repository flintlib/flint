/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz_mod_poly_factor.h"

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
    flint_rand_t state;
    fmpz_mod_poly_t f, g;
    fmpz_mod_poly_factor_t res;
    fmpz_t p;
    mpz_t pz, curr;
    int i, j, k, n, num;
    double t, T1, T2, T3;

    const len_t degs[] = {8, 16, 32, 64, 128, 256, 512, 1024};
    const int iter_count[] = {10000, 5000, 1000, 500, 300, 100, 50, 20};

    flint_randinit(state);

    mpz_init(pz);
    mpz_init(curr);
    fmpz_init(p);

    printf("Random polynomials\n");
    mpz_set_ui(pz, 2);
    mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        printf("========== p: "); fmpz_print(p); printf(" ==========\n");
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = degs[j];
            printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                fmpz_mod_poly_init(f, p);
                fmpz_mod_poly_randtest_not_zero(f, state, n);

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_berlekamp(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f);
            }

            printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        mpz_mul_ui(curr, curr, 10);
    }

    /* This code checks whether fmpz_mod_poly_factor
       made a correct choice between CZ and KS */

    printf("Check choice correctness\n");
    mpz_set_ui(pz, 2);
    mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        printf("========== p: "); fmpz_print(p); printf(" ==========\n");
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = degs[j];
            printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                fmpz_mod_poly_init(f, p);
                fmpz_mod_poly_randtest_not_zero(f, state, n);

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f);
            }

            printf("CZ: %.2lf F: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        mpz_mul_ui(curr, curr, 10);
    }

    printf("Irreducible polynomials\n");
    mpz_set_ui(pz, 2);
    mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        printf("========== p: "); fmpz_print(p); printf(" ==========\n");
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = degs[j];
            printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                fmpz_mod_poly_init(f, p);
                fmpz_mod_poly_randtest_irreducible(f, state, n);

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_berlekamp(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f);
            }

            printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        mpz_mul_ui(curr, curr, 10);
    }

    printf("Product of two irreducible polynomials\n");
    mpz_set_ui(pz, 2);
    mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        printf("========== p: "); fmpz_print(p); printf(" ==========\n");
        fflush(stdout);

        for (j = 0; j < ND; j++)
        {
            n = (degs[j] >> 1);
            printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                fmpz_mod_poly_init(f, p);
                fmpz_mod_poly_init(g, p);
                fmpz_mod_poly_randtest_irreducible(f, state, n);
                fmpz_mod_poly_randtest_irreducible(g, state, n);
                fmpz_mod_poly_mul(f, f, g);

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_berlekamp(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f);
                fmpz_mod_poly_clear(g);
            }

            printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        mpz_mul_ui(curr, curr, 10);
    }

    printf("Product of 8 small irreducible polynomials\n");
    mpz_set_ui(pz, 2);
    mpz_set_ui(curr, 10);
    for (i = 0; i < NP; i++)
    {
        fmpz_set_mpz(p, pz);
        printf("========== p: "); fmpz_print(p); printf(" ==========\n");
        fflush(stdout);

        for (j = 1; j < ND; j++)
        {
            n = (degs[j] >> 3);
            printf(">>>>>n: %d\n", n);
            fflush(stdout);

            T1 = 0;
            T2 = 0;
            T3 = 0;
            for (k = 0; k < iter_count[j]; k++)
            {
                fmpz_mod_poly_init(f, p);
                fmpz_mod_poly_init(g, p);
                fmpz_mod_poly_randtest_irreducible(f, state, n);
                for (num = 1; num < 8; num++)
                {
                    fmpz_mod_poly_randtest_irreducible(g, state, n);
                    fmpz_mod_poly_mul(f, f, g);
                }

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_cantor_zassenhaus(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T1 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_berlekamp(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T2 += t;

                t = clock();
                fmpz_mod_poly_factor_init(res);
                fmpz_mod_poly_factor_kaltofen_shoup(res, f);
                fmpz_mod_poly_factor_clear(res);
                t = (clock() - t) / CLOCKS_PER_SEC;
                T3 += t;

                fmpz_mod_poly_clear(f);
                fmpz_mod_poly_clear(g);
            }

            printf("CZ: %.2lf B: %.2lf KS: %.2lf\n", T1, T2, T3);
            fflush(stdout);

            if (T1 > T3 + 1)
                break;
        }

        mpz_nextprime(pz, curr);
        mpz_mul_ui(curr, curr, 10);
    }

    mpz_clear(pz);
    mpz_clear(curr);
    fmpz_clear(p);
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

