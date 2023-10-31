/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "mpoly.h"

TEST_FUNCTION_START(mpoly_monomial_halves, state)
{
    slong i;

    for (i = 0; i < 10000*flint_test_multiplier(); i++)
    {
        mpoly_ctx_t ctx;
        slong N, v, nvars = 1 + n_randint(state, 20);
        fmpz* e;
        ulong* e1, * e2;
        flint_bitcnt_t bits;
        int all_are_even, res;

        mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        e = _fmpz_vec_init(nvars);
        _fmpz_vec_randtest_unsigned(e, state, nvars, n_randint(state, 200));

        bits = mpoly_exp_bits_required_ffmpz(e, ctx);
        bits = mpoly_fix_bits(bits, ctx);

        N = mpoly_words_per_exp(bits, ctx);
        e1 = FLINT_ARRAY_ALLOC(N, ulong);
        e2 = FLINT_ARRAY_ALLOC(N, ulong);

        mpoly_set_monomial_ffmpz(e1, e, bits, ctx);

        if (bits <= FLINT_BITS)
            res = mpoly_monomial_overflows(e1, N, mpoly_overflow_mask_sp(bits));
        else
            res = mpoly_monomial_overflows_mp(e1, N, bits);

        if (res)
        {
            flint_printf("FAIL: input should not overflow\n");
            fflush(stdout);
            flint_abort();
        }

        all_are_even = 1;
        for (v = 0; v < nvars; v++)
            if (fmpz_is_odd(e + v))
                all_are_even = 0;

        if (bits <= FLINT_BITS)
            res = mpoly_monomial_halves(e2, e1, N, mpoly_overflow_mask_sp(bits));
        else
            res = mpoly_monomial_halves_mp(e2, e1, N, bits);

        if (all_are_even)
        {
            if (!res)
            {
                flint_printf("FAIL: should have been able to halve exponent\n");
                flint_printf("nvars = %wd, bits = %wu\n", nvars, bits);
                fflush(stdout);
                flint_abort();
            }

            if (bits <= FLINT_BITS)
                res = mpoly_monomial_overflows(e2, N, mpoly_overflow_mask_sp(bits));
            else
                res = mpoly_monomial_overflows_mp(e2, N, bits);

            if (res)
            {
                flint_printf("FAIL: half overflows\n");
                flint_printf("nvars = %wd, bits = %wu\n", nvars, bits);
                fflush(stdout);
                flint_abort();
            }

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(e2, e2, e2, N);
            else
                mpoly_monomial_add_mp(e2, e2, e2, N);

            res = mpoly_monomial_equal(e1, e2, N);

            if (!res)
            {
                flint_printf("FAIL: double of half doesn't match\n");
                flint_printf("nvars = %wd, bits = %wu\n", nvars, bits);
                fflush(stdout);
                flint_abort();
            }

        }
        else
        {
            if (res)
            {
                flint_printf("FAIL: shouldn't have been able to halve exponent\n");
                flint_printf("nvars = %wd, bits = %wu\n", nvars, bits);
                fflush(stdout);
                flint_abort();
            }
        }

        mpoly_ctx_clear(ctx);
        _fmpz_vec_clear(e, nvars);
        flint_free(e1);
        flint_free(e2);
    }

    TEST_FUNCTION_END(state);
}
