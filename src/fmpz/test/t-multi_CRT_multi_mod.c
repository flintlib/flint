/*
    Copyright (C) 2019-2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_multi_CRT_multi_mod, state)
{
    slong i, j, k;

    /* test internal interface */
    {
        fmpz_multi_CRT_t P;
        fmpz * moduli, * inputs, * outputs;
        slong moduli_count = 1000;

        moduli = _fmpz_vec_init(moduli_count);
        inputs = _fmpz_vec_init(moduli_count);
        outputs = _fmpz_vec_init(moduli_count);

        for (k = 0; k < moduli_count; k++)
        {
            fmpz_set_ui(moduli + k, n_nth_prime(k + 2));
            fmpz_sub_ui(inputs + k, moduli + k, 1);
        }

        fmpz_multi_CRT_init(P);
        if (!fmpz_multi_CRT_precompute(P, moduli, moduli_count))
        {
            flint_printf("FAIL: Check internal interface precompute\n");
            fflush(stdout);
            flint_abort();
        }

        FLINT_ASSERT(P->localsize <= moduli_count);

        for (k = 0; k < 1; k++)
        {
            _fmpz_multi_CRT_precomp(outputs, P, inputs, 1);
        }

        _fmpz_vec_clear(moduli, moduli_count);
        _fmpz_vec_clear(inputs, moduli_count);
        _fmpz_vec_clear(outputs, moduli_count);

        fmpz_multi_CRT_clear(P);
    }

    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_multi_CRT_t P;
        fmpz_multi_mod_t M;
        fmpz_t t, p;
        slong moduli_count;
        flint_bitcnt_t moduli_bits;
        fmpz * moduli, * inputs, * values;
        fmpz_t output;

        fmpz_init(t);
        fmpz_init(p);
        fmpz_init(output);

        fmpz_multi_CRT_init(P);
        fmpz_multi_mod_init(M);

        for (j = 0; j < 5; j++)
        {
            moduli_count = n_randint(state, 40) + 1;

            moduli = _fmpz_vec_init(moduli_count);
            inputs = _fmpz_vec_init(moduli_count);
            values = _fmpz_vec_init(moduli_count);
            for (k = 0; k < moduli_count; k++)
            {
                moduli_bits = n_randint(state, 150) + 1;
                fmpz_init(moduli + k);
                fmpz_init(inputs + k);
                fmpz_randtest_unsigned(moduli + k, state, moduli_bits);
                fmpz_randtest(inputs + k, state, moduli_bits);
                if (n_randint(state, 2))
                    fmpz_nextprime(moduli + k, moduli + k, 1);
            }

            if (fmpz_multi_CRT_precompute(P, moduli, moduli_count))
            {
                int sign = 1;

                if (!fmpz_multi_mod_precompute(M, moduli, moduli_count))
                {
                    flint_printf("FAIL: Check multi mod precompute\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_multi_CRT_precomp(output, P, inputs, sign);

                fmpz_multi_mod_precomp(values, M, output, sign);

                fmpz_one(p);
                for (k = 0; k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);

                    fmpz_sub(t, output, inputs + k);
                    if (!fmpz_divisible(t, moduli + k))
                    {
                        flint_printf("FAIL: Check output remainder\n");
                        flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }

                    fmpz_sub(t, inputs + k, values + k);
                    if (!fmpz_divisible(t, moduli + k))
                    {
                        flint_printf("FAIL: Check multi mod\n");
                        flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }

                    if (sign ? (fmpz_cmp2abs(moduli + k, values + k) < 0) :
                               (fmpz_cmp(values + k, moduli + k) >= 0 ||
                                                     fmpz_sgn(values + k) < 0))
                    {
                        flint_printf("FAIL: Check mod output size\n");
                        flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }

                if (sign ? (fmpz_cmp2abs(p, output) < 0) :
                           (fmpz_cmp(output, p) >= 0 || fmpz_sgn(output) < 0))
                {
                    flint_printf("FAIL: Check CRT output size\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_abs(p, p);
                fmpz_add(p, p, p);
                fmpz_add_ui(p, p, 3);
                fmpz_randtest_mod_signed(output, state, p);

                fmpz_multi_mod_precomp(values, M, output, sign);

                fmpz_one(p);
                for (k = 0; k < moduli_count; k++)
                {
                    fmpz_sub(t, output, values + k);
                    if (!fmpz_divisible(t, moduli + k))
                    {
                        flint_printf("FAIL: Check multi mod 2\n");
                        flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                fmpz_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);
                    ok |= (fmpz_is_zero(moduli + k) || fmpz_is_pm1(moduli + k));
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_divexact(t, p, moduli + k);
                    fmpz_gcd(t, t, moduli + k);
                    ok |= (fmpz_cmp_ui(t, 1) > 0);
                }

                if (!ok)
                {
                    flint_printf("FAIL: Check flat CRT failure\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            _fmpz_vec_clear(moduli, moduli_count);
            _fmpz_vec_clear(inputs, moduli_count);
            _fmpz_vec_clear(values, moduli_count);
        }

        fmpz_clear(t);
        fmpz_clear(p);
        fmpz_clear(output);

        fmpz_multi_CRT_clear(P);
        fmpz_multi_mod_clear(M);
    }

    /* test lazy interface */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_t t, p;
        slong moduli_count;
        flint_bitcnt_t moduli_bits;
        fmpz * moduli, * inputs;
        fmpz_t output;

        fmpz_init(t);
        fmpz_init(p);
        fmpz_init(output);

        for (j = 0; j < 4; j++)
        {
            moduli_bits = n_randint(state, 100) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = (fmpz *) flint_malloc(moduli_count*sizeof(fmpz));
            inputs = (fmpz *) flint_malloc(moduli_count*sizeof(fmpz));
            for (k = 0; k < moduli_count; k++)
            {
                fmpz_init(moduli + k);
                fmpz_init(inputs + k);
                fmpz_randtest_unsigned(moduli + k, state, moduli_bits);
                fmpz_randtest(inputs + k, state, moduli_bits);
            }

            if (fmpz_multi_CRT(output, moduli, inputs, moduli_count, 1))
            {
                fmpz_one(p);
                for (k = 0; k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);

                    fmpz_sub(t, output, inputs + k);
                    if (!fmpz_divisible(t, moduli + k))
                    {
                        printf("FAIL\n");
                        flint_printf("Check remainder lazy "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                fmpz_add(output, output, output);
                if (fmpz_cmpabs(output, p) > 0)
                {
                    printf("FAIL\n");
                    flint_printf("Check output size lazy "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                fmpz_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);
                    if (fmpz_is_zero(moduli + k) || fmpz_is_pm1(moduli + k))
                    {
                        ok = 1;
                    }
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_divexact(t, p, moduli + k);
                    fmpz_gcd(t, t, moduli + k);
                    if (fmpz_cmp_ui(t, 1) > 0)
                    {
                        ok = 1;
                    }
                }

                if (!ok)
                {
                    printf("FAIL\n");
                    flint_printf("Check lazy CRT failure i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            for (k = 0; k < moduli_count; k++)
            {
                fmpz_clear(moduli + k);
                fmpz_clear(inputs + k);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        fmpz_clear(t);
        fmpz_clear(p);
        fmpz_clear(output);
    }

    TEST_FUNCTION_END(state);
}
