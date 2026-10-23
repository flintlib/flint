/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "aprcl.h"

TEST_FUNCTION_START(aprcl_config_jacobi, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t n, s2;
        aprcl_config conf;

        fmpz_init(n);
        fmpz_init(s2);
        fmpz_randtest_unsigned(n, state, 1000);

        aprcl_config_jacobi_init(conf, n);

        fmpz_mul(s2, conf->s, conf->s);
        if (fmpz_cmp(s2, n) <= 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("s^2 = ");
            fmpz_print(s2);
            flint_printf(" <= ");
            fmpz_print(n);
            fflush(stdout);
            flint_abort();
        }

        aprcl_config_jacobi_clear(conf);
        fmpz_clear(s2);
        fmpz_clear(n);
    }

#if FLINT_BITS == 64
    /* beyond the tabulated values of R (about 20004 bits), R is extended */
    for (i = 0; i < 3; i++)
    {
        fmpz_t n, s2;
        aprcl_config conf;
        slong j;

        fmpz_init(n);
        fmpz_init(s2);
        fmpz_one(n);
        fmpz_mul_2exp(n, n, (i == 0) ? 20003 : (i == 1) ? 20004 : 30000);
        fmpz_add_ui(n, n, 1);

        aprcl_config_jacobi_init(conf, n);

        fmpz_mul(s2, conf->s, conf->s);
        if (fmpz_cmp(s2, n) <= 0)
        {
            flint_printf("FAIL (large n):\n");
            flint_printf("bits = %wu, R = %wu\n", fmpz_bits(n), conf->R);
            fflush(stdout);
            flint_abort();
        }

        /* s is a product of powers of primes q with q - 1 | R */
        fmpz_one(s2);
        for (j = 0; j < conf->qs->num; j++)
        {
            ulong q = fmpz_get_ui(conf->qs->p + j);

            if (!n_is_prime(q) || conf->R % (q - 1) != 0)
            {
                flint_printf("FAIL (large n): q = %wu, R = %wu\n", q, conf->R);
                fflush(stdout);
                flint_abort();
            }

            if (conf->qs_used[j])
            {
                fmpz_t t;
                fmpz_init(t);
                fmpz_pow_ui(t, conf->qs->p + j, conf->qs->exp[j]);
                fmpz_mul(s2, s2, t);
                fmpz_clear(t);
            }
        }

        if (!fmpz_equal(s2, conf->s))
        {
            flint_printf("FAIL (large n): s is not the product of the used q\n");
            fflush(stdout);
            flint_abort();
        }

        aprcl_config_jacobi_clear(conf);
        fmpz_clear(s2);
        fmpz_clear(n);
    }
#endif

    TEST_FUNCTION_END(state);
}
