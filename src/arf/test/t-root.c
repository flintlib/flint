/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "test_helpers.h"
#include "arf.h"

int
arf_root_naive(arf_t z, const arf_t x, int k, slong prec, arf_rnd_t rnd)
{
    fmpz_t t;
    mpfr_t u, v;
    int res;

    if ((k == 0) || (k > 1 && arf_sgn(x) < 0))
    {
        arf_nan(z);
        return 0;
    }

    fmpz_init(t);

    if (arf_is_special(x))
        fmpz_zero(t);
    else
        fmpz_fdiv_q_ui(t, ARF_EXPREF(x), k);
    fmpz_mul_ui(t, t, k);
    fmpz_neg(t, t);
    arf_mul_2exp_fmpz(z, x, t);

    mpfr_init2(u, FLINT_MAX(2, arf_bits(z)));
    mpfr_init2(v, prec);

    arf_get_mpfr(u, z, MPFR_RNDD);

    res = (mpfr_rootn_ui(v, u, k, arf_rnd_to_mpfr(rnd)) != 0);
    arf_set_mpfr(z, v);

    fmpz_neg(t, t);
    fmpz_tdiv_q_ui(t, t, k);
    arf_mul_2exp_fmpz(z, z, t);

    mpfr_clear(u);
    mpfr_clear(v);

    fmpz_clear(t);
    return res;
}

TEST_FUNCTION_START(arf_root, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, z, v;
        slong prec, r1, r2;
        ulong k;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(z);
        arf_init(v);

        for (iter2 = 0; iter2 < 10; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            prec = 2 + n_randint(state, 2000);
            k = n_randint(state, 50);

            if (n_randint(state, 20) == 0)
                arf_mul(x, x, x, prec, ARF_RND_DOWN);
            else if (n_randint(state, 20) == 0)
                arf_mul(x, x, x, prec, ARF_RND_UP);

            switch (n_randint(state, 4))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                default: rnd = ARF_RND_CEIL; break;
            }

            switch (n_randint(state, 2))
            {
            case 0:
                r1 = arf_root(z, x, k, prec, rnd);
                r2 = arf_root_naive(v, x, k, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("k = %wu, prec = %wd, rnd = %d\n\n", k, prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r2 = arf_root_naive(v, x, k, prec, rnd);
                r1 = arf_root(x, x, k, prec, rnd);
                if (!arf_equal(v, x) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing)!\n");
                    flint_printf("k = %wu, prec = %wd, rnd = %d\n\n", k, prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        arf_clear(x);
        arf_clear(z);
        arf_clear(v);
    }

    TEST_FUNCTION_END(state);
}
