/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

int
arf_add_naive(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    if (!arf_is_special(x) && !arf_is_special(y))
    {
        fmpz_t t;
        int inexact;

        fmpz_init(t);
        fmpz_sub(t, ARF_EXPREF(x), ARF_EXPREF(y));

        if (fmpz_cmp_si(t, 5000) > 0)
        {
            arf_t eps;
            arf_init(eps);
            arf_mul_2exp_si(eps, x, -5000);
            arf_abs(eps, eps);
            if (arf_sgn(y) < 0)
                arf_neg(eps, eps);
            arf_add(z, x, eps, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_set_round(z, z, prec, rnd);
            arf_clear(eps);
        }
        else if (fmpz_cmp_si(t, -5000) < 0)
        {
            arf_t eps;
            arf_init(eps);
            arf_mul_2exp_si(eps, y, -5000);
            arf_abs(eps, eps);
            if (arf_sgn(x) < 0)
                arf_neg(eps, eps);
            arf_add(z, eps, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_set_round(z, z, prec, rnd);
            arf_clear(eps);
        }
        else
        {
            arf_add(z, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_set_round(z, z, prec, rnd);
        }

        fmpz_clear(t);
        return inexact;
    }
    else
    {
        arf_add(z, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
        return arf_set_round(z, z, prec, rnd);
    }
}

TEST_FUNCTION_START(arf_add, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, v;
        slong prec, r1, r2;
        arf_rnd_t rnd;
        fmpz_t t;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(v);
        fmpz_init(t);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            arf_randtest_special(y, state, 2000, 100);
            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                case 3:  rnd = ARF_RND_CEIL; break;
                default: rnd = ARF_RND_NEAR; break;
            }

            rnd = ARF_RND_DOWN;

            if (arf_is_normal(x) && arf_is_normal(y))
            {
                fmpz_sub(t, ARF_EXPREF(x), ARF_EXPREF(y));

                /* if not too far apart, sometimes test exact addition */
                if (fmpz_bits(t) < 10)
                {
                    if (n_randint(state, 10) == 0)
                        prec = ARF_PREC_EXACT;
                }
            }

            switch (n_randint(state, 5))
            {
            case 0:
                r1 = arf_add(z, x, y, prec, rnd);
                r2 = arf_add_naive(v, x, y, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 1:
                r1 = arf_add(z, x, x, prec, rnd);
                r2 = arf_add_naive(v, x, x, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 2:
                r2 = arf_add_naive(v, x, x, prec, rnd);
                r1 = arf_add(x, x, x, prec, rnd);
                if (!arf_equal(v, x) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 3:
                r2 = arf_add_naive(v, x, y, prec, rnd);
                r1 = arf_add(x, x, y, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r2 = arf_add_naive(v, x, y, prec, rnd);
                r1 = arf_add(x, y, x, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(v);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
