/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_reconstruct_fmpz_2, state)
{
    int i;

    /* check successful reconstructions */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int success;
        fmpq_t x, y;
        fmpz_t mod, res, N, D, t;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(mod);
        fmpz_init(res);
        fmpz_init(N);
        fmpz_init(D);
        fmpz_init(t);

        if (i % 2)
        {
            fmpq_randtest(x, state, 1 + n_randint(state, 5000));
        }
        else
        {
            slong k, bound;
            fmpz * c1;

            bound = 1 + n_randint(state, 100);
            c1 = _fmpz_vec_init(bound);

            fmpz_randtest(c1 + 0, state, 2*FLINT_BITS);
            for (k = 1; k < bound; k++)
            {
                fmpz_randtest_unsigned(c1 + k, state, 3*FLINT_BITS);
                fmpz_add_ui(c1 + k, c1 + k, 1);
            }

            fmpq_set_cfrac(x, c1, bound);

            _fmpz_vec_clear(c1, bound);
        }

        fmpz_abs(N, fmpq_numref(x));
        fmpz_set(D, fmpq_denref(x));

        /* Randomly generate larger bounds */
        if (n_randint(state, 2))
        {
            fmpz_randtest_not_zero(t, state, 2000);
            fmpz_abs(t, t);
            fmpz_mul(N, N, t);
        }
        if (n_randint(state, 2))
        {
            fmpz_randtest_not_zero(t, state, 2000);
            fmpz_abs(t, t);
            fmpz_mul(D, D, t);
        }

        if (n_randint(state, 10))
        {
            fmpz_mul(mod, N, D);
            fmpz_mul_ui(mod, mod, 2);
        }
        else
        {
            fmpz_randbits(mod, state, FLINT_BITS * FMPQ_RECONSTRUCT_HGCD_CUTOFF + fmpz_bits(N) + n_randint(state, 10));
            fmpz_abs(mod, mod);
        }

        do fmpz_add_ui(mod, mod, 1);
        while (!fmpq_mod_fmpz(res, x, mod));

        success = fmpq_reconstruct_fmpz_2(y, res, mod, N, D);

        if (!success || !fmpq_equal(x, y))
            flint_throw(FLINT_TEST_FAIL,
                    "success = %d\n"
                    "x = %{fmpq}\n"
                    "N = %{fmpz}\n"
                    "D = %{fmpz}\n"
                    "mod = %{fmpz}\n"
                    "res = %{fmpz}\n"
                    "y = %{fmpq}\n",
                    success, x, N, D, mod, res, y);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(mod);
        fmpz_clear(res);
        fmpz_clear(N);
        fmpz_clear(D);
        fmpz_clear(t);
    }

    /* check random reconstructions */
    for (i = 0; i < 1000*flint_test_multiplier(); i++)
    {
        int success1, success2;
        fmpq_t x, y;
        fmpz_t mod, res, N, D, t;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(mod);
        fmpz_init(res);
        fmpz_init(N);
        fmpz_init(D);
        fmpz_init(t);

        if (i % 2)
        {
            fmpz_randtest_not_zero(mod, state, 5000);
            fmpz_randtest_not_zero(res, state, 5000);
        }
        else
        {
            slong k, bound;
            fmpz * c1;

            bound = 2 + n_randint(state, 100);
            c1 = _fmpz_vec_init(bound);

            fmpz_zero(c1 + 0);
            for (k = 1; k < bound; k++)
            {
                fmpz_randtest_unsigned(c1 + k, state, 2*FLINT_BITS);
                fmpz_add_ui(c1 + k, c1 + k, 1);
            }

            fmpq_set_cfrac(x, c1, bound);

            _fmpz_vec_clear(c1, bound);

            fmpz_swap(res, fmpq_numref(x));
            fmpz_swap(mod, fmpq_denref(x));
        }

        if (fmpz_cmp_ui(mod, 3) < 0)
            fmpz_set_ui(mod, 3);

        fmpz_mod(res, res, mod);

        fmpz_sub_ui(t, mod, 1);
        fmpz_fdiv_q_ui(t, t, 2);

        fmpz_randtest_mod(N, state, t);
        fmpz_add_ui(N, N, 1);

        fmpz_fdiv_q(D, t, N);
        fmpz_randtest_mod(D, state, D);
        fmpz_add_ui(D, D, 1);

        success1 = _fmpq_reconstruct_fmpz_2(fmpq_numref(x),
                                               fmpq_denref(x), res, mod, N, D);

        success2 = _fmpq_reconstruct_fmpz_2_naive(fmpq_numref(y),
                                               fmpq_denref(y), res, mod, N, D);

        if (success1 != success2 || (success1 && !fmpq_equal(x, y)))
            flint_throw(FLINT_TEST_FAIL, "check match with naive: i = %wd\n", i);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(mod);
        fmpz_clear(res);
        fmpz_clear(N);
        fmpz_clear(D);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
