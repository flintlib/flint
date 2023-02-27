/*
    Copyright (C) 2012, 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

static slong
fmpr_mul_fmpz_naive(fmpr_t z, const fmpr_t x, const fmpz_t y,
    slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_fmpz(t, y);
    r = fmpr_mul(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

int main()
{
    slong iter, iter2;
    flint_rand_t state;

    flint_printf("mul_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 300 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, z, v;
        fmpz_t y;
        slong prec, r1, r2;
        fmpr_rnd_t rnd;

        fmpr_init(x);
        fmpr_init(z);
        fmpr_init(v);
        fmpz_init(y);

        for (iter2 = 0; iter2 < 30; iter2++)
        {
            fmpr_randtest_special(x, state, 2000, 200);
            fmpz_randtest(y, state, 2000);
            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 4))
            {
                case 0:  rnd = FMPR_RND_DOWN; break;
                case 1:  rnd = FMPR_RND_UP; break;
                case 2:  rnd = FMPR_RND_FLOOR; break;
                default: rnd = FMPR_RND_CEIL; break;
            }

            switch (n_randint(state, 2))
            {
            case 0:
                r1 = fmpr_mul_fmpz(z, x, y, prec, rnd);
                r2 = fmpr_mul_fmpz_naive(v, x, y, prec, rnd);
                if (!fmpr_equal(z, v) || r1 != r2 || !fmpr_check_ulp(z, r1, prec))
                {
                    flint_printf("FAIL!\n");
                    flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
                    flint_printf("v = "); fmpr_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                fmpr_set(v, x);
                fmpr_set(z, x);
                r1 = fmpr_mul_fmpz(z, z, y, prec, rnd);
                r2 = fmpr_mul_fmpz_naive(v, v, y, prec, rnd);
                if (!fmpr_equal(z, v) || r1 != r2 || !fmpr_check_ulp(z, r1, prec))
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
                    flint_printf("v = "); fmpr_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        fmpr_clear(x);
        fmpr_clear(z);
        fmpr_clear(v);
        fmpz_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
