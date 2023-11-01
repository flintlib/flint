/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"

TEST_FUNCTION_START(fmpz_mod_discrete_log_pohlig_hellman, state)
{
    slong i, j, k;

    /* check random word-sized primes */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_t p, pm1;
        fmpz_mod_ctx_t fpctx;
        fmpz_mod_discrete_log_pohlig_hellman_t L;

        fmpz_init_set_ui(p, 2);
        fmpz_init_set_ui(pm1, 1);
        fmpz_mod_ctx_init(fpctx, p);
        fmpz_mod_discrete_log_pohlig_hellman_init(L);

        for (j = 0; j < 10; j++)
        {
            double score;

            fmpz_set_ui(p, n_randtest_prime(state, 1));
            fmpz_sub_ui(pm1, p, 1);
            fmpz_mod_ctx_set_modulus(fpctx, p);
            score = fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(L, p);
            if (score > 10000)
            {
                continue;
            }

            for (k = 0; k < 10; k++)
            {
                fmpz_t x, y, xr;
                const fmpz * alpha = fmpz_mod_discrete_log_pohlig_hellman_primitive_root(L);
                fmpz_init(x);
                fmpz_init(y);
                fmpz_init(xr);

                fmpz_randm(x, state, pm1);
                fmpz_mod_pow_fmpz(y, alpha, x, fpctx);
                fmpz_mod_discrete_log_pohlig_hellman_run(xr, L, y);
                if (!fmpz_equal(x, xr))
                {
                    printf("FAIL\ncheck random word-sized primes\n");
                    flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_clear(x);
                fmpz_clear(y);
                fmpz_clear(xr);
            }
        }
        fmpz_mod_discrete_log_pohlig_hellman_clear(L);
        fmpz_mod_ctx_clear(fpctx);
        fmpz_clear(p);
        fmpz_clear(pm1);
    }

    /* check smooth primes */
    for (i = 0; i < flint_test_multiplier(); i++)
    {
        fmpz_t p;
        int success;
        double score;
        fmpz_mod_ctx_t fpctx;
        fmpz_mod_discrete_log_pohlig_hellman_t L;

        fmpz_init_set_ui(p, 2);
        fmpz_mod_ctx_init(fpctx, p);
        fmpz_mod_discrete_log_pohlig_hellman_init(L);

        for (j = 0; j < 1000; j++)
        {
            success = fmpz_next_smooth_prime(p, p);
            if (!success)
            {
                break;
            }
            if (1 != fmpz_is_prime(p))
            {
                printf("FAIL\nprimality test failed p = ");
                fmpz_print(p);
                fflush(stdout);
                flint_abort();
            }

            /*
                The primes returned by fmpz_next_smooth_prime should not
                require too many operations for discrete logarithms.
            */
            fmpz_mod_ctx_set_modulus(fpctx, p);
            score = fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(L, p);
            if (score > 10000)
            {
                printf("FAIL\nsmooth prime score %f too high\np = ", score);
                fmpz_print(p);
                fflush(stdout);
                flint_abort();
            }

            for (k = 0; k < 10; k++)
            {
                fmpz_t x, y, xr, pm1;
                const fmpz * alpha = fmpz_mod_discrete_log_pohlig_hellman_primitive_root(L);
                fmpz_init(x);
                fmpz_init(y);
                fmpz_init(xr);
                fmpz_init(pm1);

                fmpz_sub_ui(pm1, p, 1);

                fmpz_randm(x, state, pm1);
                fmpz_mod_pow_fmpz(y, alpha, x, fpctx);
                fmpz_mod_discrete_log_pohlig_hellman_run(xr, L, y);
                if (!fmpz_equal(x, xr))
                {
                    printf("FAIL\n");
                    flint_printf("i = %wd, j = %kw, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_clear(x);
                fmpz_clear(y);
                fmpz_clear(xr);
                fmpz_clear(pm1);
            }
        }

        fmpz_clear(p);
        fmpz_mod_ctx_clear(fpctx);
        fmpz_mod_discrete_log_pohlig_hellman_clear(L);
    }

    TEST_FUNCTION_END(state);
}
