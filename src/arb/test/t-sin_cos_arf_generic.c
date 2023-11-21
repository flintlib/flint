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

/* these functions are not exposed to the public for now,
   but it still makes sense to test them explicitly */
void arb_sin_cos_arf_rs_generic(arb_t res_sin, arb_t res_cos, const arf_t x, slong prec);
void arb_sin_cos_taylor_sum_rs(arb_t s, const arb_t x, slong N, int cosine, slong prec);

TEST_FUNCTION_START(arb_sin_cos_arf_generic, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, s1, s2, c1, c2;
        slong prec = 2 + n_randint(state, 2000);

        arb_init(x);
        arb_init(s1);
        arb_init(s2);
        arb_init(c1);
        arb_init(c2);

        if (n_randint(state, 2))
            arb_randtest(x, state, 1 + n_randint(state, 2000), 15);
        else
            arb_randtest(x, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 10000));

        mag_zero(arb_radref(x));

        arb_sin_cos(s1, c1, x, prec);

        switch (n_randint(state, 6))
        {
            case 0:
                arb_sin_cos_arf_generic(s2, c2, arb_midref(x), prec);
                break;
            case 1:
                arb_sin_cos_arf_generic(s2, NULL, arb_midref(x), prec);
                arb_sin_cos_arf_generic(NULL, c2, arb_midref(x), prec);
                break;
            case 2:
                arb_set(s2, x);
                arb_sin_cos_arf_generic(NULL, c2, arb_midref(s2), prec);
                arb_sin_cos_arf_generic(s2, NULL, arb_midref(s2), prec);
                break;
            case 3:
                arb_set(c2, x);
                arb_sin_cos_arf_generic(s2, NULL, arb_midref(c2), prec);
                arb_sin_cos_arf_generic(NULL, c2, arb_midref(c2), prec);
                break;
            case 4:
                arb_set(c2, x);
                arb_sin_cos_arf_generic(s2, c2, arb_midref(c2), prec);
                break;
            default:
                arb_set(s2, x);
                arb_sin_cos_arf_generic(s2, c2, arb_midref(s2), prec);
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

        if (arf_cmpabs_ui(arb_midref(x), 1) <= 0 &&
            (arb_rel_accuracy_bits(s2) < prec - 2 || arb_rel_accuracy_bits(c2) < prec - 2))
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

    /* test the rs algorithm explicitly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, s1, s2, c1, c2;
        slong prec = 2 + n_randint(state, 2000);

        arb_init(x);
        arb_init(s1);
        arb_init(s2);
        arb_init(c1);
        arb_init(c2);

        arb_randtest(x, state, 1 + n_randint(state, 2000), 0);
        mag_zero(arb_radref(x));

        if (n_randint(state, 2))
            arb_mul_2exp_si(x, x, -n_randint(state, 2 * prec));

        while (arf_cmpabs_d(arb_midref(x), 3.141) > 0)
            arb_mul_2exp_si(x, x, -1);

        arb_sin_cos(s1, c1, x, prec);

        switch (n_randint(state, 6))
        {
            case 0:
                arb_sin_cos_arf_rs_generic(s2, c2, arb_midref(x), prec);
                break;
            case 1:
                arb_sin_cos_arf_rs_generic(s2, NULL, arb_midref(x), prec);
                arb_sin_cos_arf_rs_generic(NULL, c2, arb_midref(x), prec);
                break;
            case 2:
                arb_set(s2, x);
                arb_sin_cos_arf_rs_generic(NULL, c2, arb_midref(s2), prec);
                arb_sin_cos_arf_rs_generic(s2, NULL, arb_midref(s2), prec);
                break;
            case 3:
                arb_set(c2, x);
                arb_sin_cos_arf_rs_generic(s2, NULL, arb_midref(c2), prec);
                arb_sin_cos_arf_rs_generic(NULL, c2, arb_midref(c2), prec);
                break;
            case 4:
                arb_set(c2, x);
                arb_sin_cos_arf_rs_generic(s2, c2, arb_midref(c2), prec);
                break;
            default:
                arb_set(s2, x);
                arb_sin_cos_arf_rs_generic(s2, c2, arb_midref(s2), prec);
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

    /* test the series evaluation code directly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        slong prec;
        slong N;
        int cosine;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        prec = 2 + n_randint(state, 2000);
        N = n_randint(state, 100);
        cosine = n_randint(state, 2);
        arb_randtest(x, state, 1 + n_randint(state, 2000), 1);
        mag_zero(arb_radref(x));

        if (n_randint(state, 2))
            arb_mul_2exp_si(x, x, -n_randint(state, 2 * prec));

        if (cosine)
        {
            arb_cos(y, x, prec);
            arb_sin_cos_taylor_sum_rs(z, x, N, 1, prec);
        }
        else
        {
            arb_sin(y, x, prec);
            arb_sin_cos_taylor_sum_rs(z, x, N, 0, prec);
        }

        if (!arb_overlaps(z, y))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec = %wd, N = %wd, cosine = %d\n\n", prec, N, cosine);
            flint_printf("x = "); arb_printn(x, 500, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 500, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
