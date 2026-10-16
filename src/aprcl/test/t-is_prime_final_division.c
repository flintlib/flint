/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "fmpz.h"
#include "aprcl.h"

/* naive version: walk n^i mod s and test the residues 1 < a <= sqrt(n) */
static int
final_division_naive(const fmpz_t n, const fmpz_t s, ulong r)
{
    int result = 1;
    ulong i;
    fmpz_t npow, nmul, rem, bound;

    fmpz_init(rem);
    fmpz_init(bound);
    fmpz_init(nmul);
    fmpz_init(npow);
    fmpz_sqrt(bound, n);
    fmpz_mod(nmul, n, s);
    fmpz_set(npow, nmul);

    for (i = 1; i < r; i++)
    {
        if (fmpz_is_one(npow))
            break;

        if (fmpz_cmp(npow, bound) <= 0)
        {
            fmpz_mod(rem, n, npow);
            if (fmpz_is_zero(rem))
            {
                result = 0;
                break;
            }
        }

        fmpz_mul(npow, npow, nmul);
        fmpz_mod(npow, npow, s);
    }

    fmpz_clear(npow);
    fmpz_clear(nmul);
    fmpz_clear(rem);
    fmpz_clear(bound);

    return result;
}

TEST_FUNCTION_START(aprcl_is_prime_final_division, state)
{
    slong iter;
    int old_threads = flint_get_num_threads();

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_t n, a, b, g;
        aprcl_config conf;
        int r1, r2, type;
        ulong bits = 20 + n_randint(state, 300);

        fmpz_init(n);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(g);

        type = n_randint(state, 4);
        if (type == 0)
        {
            /* prime */
            fmpz_randprime(n, state, bits, 0);
        }
        else if (type == 1)
        {
            /* product of two primes of about the same size */
            fmpz_randprime(a, state, bits / 2 + 1, 0);
            fmpz_randprime(b, state, bits / 2 + 1 + n_randint(state, 3), 0);
            fmpz_mul(n, a, b);
        }
        else if (type == 2)
        {
            /* square of a prime */
            fmpz_randprime(a, state, bits / 2 + 1, 0);
            fmpz_mul(n, a, a);
        }
        else
        {
            /* random odd composite */
            fmpz_randtest_unsigned(n, state, bits);
            fmpz_setbit(n, 0);
            fmpz_add_ui(n, n, 4);
        }

        /* s, R with s^2 > n as chosen by the Jacobi configuration */
        aprcl_config_jacobi_init(conf, n);

        /* the final division requires gcd(n, s) = 1 */
        fmpz_gcd(g, n, conf->s);
        if (fmpz_is_one(g))
        {
            flint_set_num_threads(1 + n_randint(state, 4));

            r1 = aprcl_is_prime_final_division(n, conf->s, conf->R);
            r2 = final_division_naive(n, conf->s, conf->R);

            if (r1 != r2)
            {
                flint_printf("FAIL:\n");
                flint_printf("n = "); fmpz_print(n);
                flint_printf("\ns = "); fmpz_print(conf->s);
                flint_printf("\nR = %wu\n", conf->R);
                flint_printf("result = %d, naive = %d\n", r1, r2);
                fflush(stdout);
                flint_abort();
            }

            /* a prime always passes */
            if (type == 0 && r1 != 1)
            {
                flint_printf("FAIL: prime rejected\n");
                fflush(stdout);
                flint_abort();
            }
        }

        aprcl_config_jacobi_clear(conf);
        fmpz_clear(n);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(g);
    }

    /*
        Constructed cases with a divisor at a known position: for a prime s,
        any a < s and an exponent i coprime to s - 1, choose b with
        b^i = a^{1 - i} mod s; then n = a b satisfies n^i = a mod s, so the
        walk must find the divisor a at step i (and must not if r <= i).
        Sizes cover the fmpz, plain, Shoup and matrix multiplication paths;
        an s of exactly 128 bits (2 limbs, no normalisation shift) reaches
        the plain mulmod_preinvn path.
    */
    {
        fmpz_t n, s, a, b, sm1, i, e, t;
        slong j;
        const slong sizes[] = {40, 64, 128, 130, 200, 400, 700, 1000};

        fmpz_init(n); fmpz_init(s); fmpz_init(a); fmpz_init(b);
        fmpz_init(sm1); fmpz_init(i); fmpz_init(e); fmpz_init(t);

        for (j = 0; j < 4 * (slong) (sizeof(sizes) / sizeof(sizes[0])); j++)
        {
            slong sbits = sizes[j % (sizeof(sizes) / sizeof(sizes[0]))];
            ulong r = 1000 + n_randint(state, 100000), ii;
            int r1, r2, r3;

            fmpz_randprime(s, state, sbits, 0);
            fmpz_sub_ui(sm1, s, 1);

            /* a: random odd number below s, at least 3 */
            do
            {
                fmpz_randm(a, state, s);
                fmpz_setbit(a, 0);
            } while (fmpz_cmp_ui(a, 3) < 0 || fmpz_cmp(a, s) >= 0);

            /* exponent i coprime to s - 1, with 1 <= i < r */
            do
            {
                ii = 1 + n_randint(state, r - 1);
                fmpz_set_ui(i, ii);
                fmpz_gcd(t, i, sm1);
            } while (!fmpz_is_one(t));

            /* b = b0 + k s with b0 = a^{(1 - i) e} mod s, e = 1/i mod (s - 1) */
            fmpz_invmod(e, i, sm1);
            fmpz_sub(t, sm1, i);
            fmpz_add_ui(t, t, 1);           /* t = 1 - i mod (s - 1) */
            fmpz_mul(t, t, e);
            fmpz_mod(t, t, sm1);
            fmpz_powm(b, a, t, s);
            fmpz_addmul_ui(b, s, 1 + n_randint(state, 4));
            fmpz_mul(n, a, b);              /* a < s < b, so a <= sqrt(n) */

            flint_set_num_threads(1 + n_randint(state, 4));

            r1 = aprcl_is_prime_final_division(n, s, r);
            r2 = final_division_naive(n, s, r);
            r3 = aprcl_is_prime_final_division(n, s, ii);   /* divisor out of range */

            if (r1 != 0 || r2 != 0 || r3 != 1)
            {
                flint_printf("FAIL (constructed divisor):\n");
                flint_printf("s bits = %wd, i = %wu, r = %wu, results %d %d %d\n",
                             sbits, ii, r, r1, r2, r3);
                flint_printf("n = "); fmpz_print(n);
                flint_printf("\ns = "); fmpz_print(s); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(n); fmpz_clear(s); fmpz_clear(a); fmpz_clear(b);
        fmpz_clear(sm1); fmpz_clear(i); fmpz_clear(e); fmpz_clear(t);
    }

    flint_set_num_threads(old_threads);

    TEST_FUNCTION_END(state);
}
