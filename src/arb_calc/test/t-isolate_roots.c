/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"
#include "arb_calc.h"

/* sin((pi/2)x) */
static int
sin_pi2_x(arb_ptr out, const arb_t inp, void * params, slong order, slong prec)
{
    arb_ptr x;

    x = _arb_vec_init(2);

    arb_set(x, inp);
    arb_one(x + 1);

    arb_const_pi(out, prec);
    arb_mul_2exp_si(out, out, -1);
    _arb_vec_scalar_mul(x, x, 2, out, prec);
    _arb_poly_sin_series(out, x, order, order, prec);

    _arb_vec_clear(x, 2);

    return 0;
}

TEST_FUNCTION_START(arb_calc_isolate_roots, state)
{
    slong iter;

    for (iter = 0; iter < 40 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, r, a, b, maxdepth, maxeval, maxfound, prec, i, j, num;
        arf_interval_ptr blocks;
        int * info;
        arf_interval_t interval;
        arb_t t;
        fmpz_t nn;

        prec = 2 + n_randint(state, 50);

        m = n_randint(state, 80);
        r = 1 + n_randint(state, 80);
        a = m - r;
        b = m + r;

        maxdepth = 1 + n_randint(state, 60);
        maxeval = 1 + n_randint(state, 5000);
        maxfound = 1 + n_randint(state, 100);

        arf_interval_init(interval);
        arb_init(t);
        fmpz_init(nn);

        arf_set_si(&interval->a, a);
        arf_set_si(&interval->b, b);

        num = arb_calc_isolate_roots(&blocks, &info, sin_pi2_x, NULL,
            interval, maxdepth, maxeval, maxfound, prec);

        /* check that all roots are accounted for */
        for (i = a; i <= b; i++)
        {
            if (i % 2 == 0)
            {
                int found = 0;

                for (j = 0; j < num; j++)
                {
                    arf_interval_get_arb(t, blocks + j, ARF_PREC_EXACT);

                    if (arb_contains_si(t, i))
                    {
                        found = 1;
                        break;
                    }
                }

                if (!found)
                {
                    flint_printf("FAIL: missing root %wd\n", i);
                    flint_printf("a = %wd, b = %wd, maxdepth = %wd, maxeval = %wd, maxfound = %wd, prec = %wd\n",
                        a, b, maxdepth, maxeval, maxfound, prec);

                    for (j = 0; j < num; j++)
                    {
                        arf_interval_printd(blocks + j, 15);
                        flint_printf("   %d \n", info[i]);
                    }

                    flint_abort();
                }
            }
        }

        /* check that all reported single roots are good */
        for (i = 0; i < num; i++)
        {
            if (info[i] == 1)
            {
                /* b contains unique 2n -> b/2 contains unique n */
                arf_interval_get_arb(t, blocks + i, ARF_PREC_EXACT);
                arb_mul_2exp_si(t, t, -1);

                if (!arb_get_unique_fmpz(nn, t))
                {
                    flint_printf("FAIL: bad root %wd\n", i);
                    flint_printf("a = %wd, b = %wd, maxdepth = %wd, maxeval = %wd, maxfound = %wd, prec = %wd\n",
                        a, b, maxdepth, maxeval, maxfound, prec);

                    for (j = 0; j < num; j++)
                    {
                        arf_interval_printd(blocks + j, 15);
                        flint_printf("   %d \n", info[i]);
                    }

                    flint_abort();
                }
            }
        }

        _arf_interval_vec_clear(blocks, num);
        flint_free(info);

        arf_interval_clear(interval);
        arb_clear(t);
        fmpz_clear(nn);
    }

    TEST_FUNCTION_END(state);
}
