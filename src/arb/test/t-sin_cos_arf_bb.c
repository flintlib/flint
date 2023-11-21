/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

void arb_sin_cos_fmpz_div_2exp_bsplit(arb_t wsin, arb_t wcos,
    const fmpz_t x, flint_bitcnt_t r, slong prec);

TEST_FUNCTION_START(arb_sin_cos_arf_bb, state)
{
    slong iter;

    /* test the series evaluation code directly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_t x;
        slong r;
        arb_t t, s1, s2, c1, c2;
        slong prec1, prec2;

        fmpz_init(x);
        arb_init(t);
        arb_init(s1);
        arb_init(s2);
        arb_init(c1);
        arb_init(c2);

        if (n_randint(state, 100) == 0)
            flint_set_num_threads(1 + n_randint(state, 4));
        else
            flint_set_num_threads(1);

        prec1 = 10 + n_randint(state, 2000);
        prec2 = 10 + n_randint(state, 2000);
        r = n_randint(state, prec2);
        r = FLINT_MAX(r, 4);
        fmpz_randtest(x, state, r - 3);

        arb_set_fmpz(t, x);
        arb_mul_2exp_si(t, t, -r);

        arb_sin_cos(s1, c1, t, prec1);
        arb_sin_cos_fmpz_div_2exp_bsplit(s2, c2, x, r, prec2);

        if (!arb_overlaps(s1, s2) || !arb_overlaps(c1, c2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
            flint_printf("t = "); arb_printn(t, 500, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 500, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 500, 0); flint_printf("\n\n");
            flint_printf("c1 = "); arb_printn(c1, 500, 0); flint_printf("\n\n");
            flint_printf("c2 = "); arb_printn(c2, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(x);
        arb_clear(t);
        arb_clear(s1);
        arb_clear(s2);
        arb_clear(c1);
        arb_clear(c2);
    }

    /* test the bb algorithm */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, s1, s2, c1, c2;
        slong prec = 2 + n_randint(state, 4000);

        arb_init(x);
        arb_init(s1);
        arb_init(s2);
        arb_init(c1);
        arb_init(c2);

        if (n_randint(state, 100) == 0)
            flint_set_num_threads(1 + n_randint(state, 4));
        else
            flint_set_num_threads(1);

        arb_randtest(x, state, 1 + n_randint(state, 4000), 0);
        mag_zero(arb_radref(x));

        if (n_randint(state, 2))
            arb_mul_2exp_si(x, x, -n_randint(state, 2 * prec));

        while (arf_cmpabs_d(arb_midref(x), 3.141) > 0)
            arb_mul_2exp_si(x, x, -1);

        arb_sin_cos(s1, c1, x, prec);

        switch (n_randint(state, 6))
        {
            case 0:
                arb_sin_cos_arf_bb(s2, c2, arb_midref(x), prec);
                break;
            case 1:
                arb_sin_cos_arf_bb(s2, NULL, arb_midref(x), prec);
                arb_sin_cos_arf_bb(NULL, c2, arb_midref(x), prec);
                break;
            case 2:
                arb_set(s2, x);
                arb_sin_cos_arf_bb(NULL, c2, arb_midref(s2), prec);
                arb_sin_cos_arf_bb(s2, NULL, arb_midref(s2), prec);
                break;
            case 3:
                arb_set(c2, x);
                arb_sin_cos_arf_bb(s2, NULL, arb_midref(c2), prec);
                arb_sin_cos_arf_bb(NULL, c2, arb_midref(c2), prec);
                break;
            case 4:
                arb_set(c2, x);
                arb_sin_cos_arf_bb(s2, c2, arb_midref(c2), prec);
                break;
            default:
                arb_set(s2, x);
                arb_sin_cos_arf_bb(s2, c2, arb_midref(s2), prec);
                break;
        }

        if (!arb_overlaps(s1, s2) || !arb_overlaps(c1, c2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("x = "); arb_printn(x, 500, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 500, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 500, 0); flint_printf("\n\n");
            flint_printf("c1 = "); arb_printn(c1, 500, 0); flint_printf("\n\n");
            flint_printf("c2 = "); arb_printn(c2, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(s2) < prec - 2 || arb_rel_accuracy_bits(c2) < prec - 2)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec = %wd,  acc1 = %wd,  acc2 = %wd\n\n",
                prec, arb_rel_accuracy_bits(s2), arb_rel_accuracy_bits(c2));
            flint_printf("x = "); arb_printn(x, 500, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 500, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 500, 0); flint_printf("\n\n");
            flint_printf("c1 = "); arb_printn(c1, 500, 0); flint_printf("\n\n");
            flint_printf("c2 = "); arb_printn(c2, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(s1);
        arb_clear(s2);
        arb_clear(c1);
        arb_clear(c2);
    }

    TEST_FUNCTION_END(state);
}
