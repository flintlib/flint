/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_powsum_series_tree, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, a, q;
        acb_ptr z1, z2;
        slong i, n, len, prec;
        int nice;

        acb_init(s);
        acb_init(a);
        acb_init(q);

        /* "nice" inputs: moderate exact parameters with Re(a) > 0, for
           which we also check that the accuracy is not much worse than
           that of the naive algorithm */
        nice = n_randint(state, 2);

        prec = 2 + n_randint(state, 300);

        if (nice)
        {
            switch (n_randint(state, 4))
            {
                case 0:
                    acb_set_si(s, (slong) n_randint(state, 8) - 2);
                    break;
                case 1:
                    acb_set_d_d(s, 0.5, (double) n_randint(state, 40));
                    break;
                default:
                    arb_set_d(acb_realref(s), (double) n_randint(state, 1000) / 100 - 3);
                    arb_set_d(acb_imagref(s), n_randint(state, 2) ? 0.0 : (double) n_randint(state, 1000) / 100 - 5);
            }

            if (n_randint(state, 2))
                acb_one(a);
            else if (n_randint(state, 2))
                acb_set_d(a, (double) (1 + n_randint(state, 1000)) / 64);
            else
                acb_set_d_d(a, (double) (1 + n_randint(state, 1000)) / 64, (double) n_randint(state, 1000) / 64 - 8);

            if (n_randint(state, 2))
                acb_one(q);
            else if (n_randint(state, 2))
                acb_set_si(q, -1);
            else
                acb_set_d_d(q, (double) (1 + n_randint(state, 100)) / 64, (double) n_randint(state, 100) / 64 - 0.75);

            /* sometimes make them inexact */
            if (n_randint(state, 4) == 0)
            {
                mag_set_ui_2exp_si(arb_radref(acb_realref(s)), 1, -(slong) n_randint(state, prec + 20));
                if (n_randint(state, 2))
                    mag_set_ui_2exp_si(arb_radref(acb_imagref(s)), 1, -(slong) n_randint(state, prec + 20));
                if (n_randint(state, 2))
                    mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -(slong) n_randint(state, prec + 20) - 8);
                if (n_randint(state, 2))
                    mag_set_ui_2exp_si(arb_radref(acb_imagref(q)), 1, -(slong) n_randint(state, prec + 20) - 8);
            }
        }
        else
        {
            acb_randtest(s, state, 1 + n_randint(state, 200), 3);

            if (n_randint(state, 2))
                acb_one(a);
            else
                acb_randtest(a, state, 1 + n_randint(state, 200), 3);

            if (n_randint(state, 2))
                acb_one(q);
            else
                acb_randtest(q, state, 1 + n_randint(state, 200), 3);
        }

        n = n_randint(state, 100);
        len = 1 + n_randint(state, 30);

        if (n_randint(state, 10) == 0)
        {
            n = n_randint(state, 400);
            len = 1 + n_randint(state, 200);
        }

        z1 = _acb_vec_init(len);
        z2 = _acb_vec_init(len);

        _acb_poly_powsum_series_naive(z1, s, a, q, n, len, prec);
        flint_set_num_threads(1 + n_randint(state, 3));
        _acb_poly_powsum_series_tree(z2, s, a, q, n, len, prec);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(z1 + i, z2 + i) ||
                (nice && acb_rel_accuracy_bits(z2 + i) < acb_rel_accuracy_bits(z1 + i) - 10
                      && acb_rel_accuracy_bits(z1 + i) > 10))
            {
                flint_printf("FAIL: overlap or accuracy\n\n");
                flint_printf("iter = %wd\n", iter);
                flint_printf("n = %wd, prec = %wd, len = %wd, i = %wd\n\n", n, prec, len, i);
                flint_printf("s = "); acb_printd(s, prec / 3.33); flint_printf("\n\n");
                flint_printf("a = "); acb_printd(a, prec / 3.33); flint_printf("\n\n");
                flint_printf("q = "); acb_printd(q, prec / 3.33); flint_printf("\n\n");
                flint_printf("z1 = "); acb_printd(z1 + i, prec / 3.33); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2 + i, prec / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(s);
        acb_clear(q);
        _acb_vec_clear(z1, len);
        _acb_vec_clear(z2, len);
    }

    /* zeta series with enough derivatives that the Euler-Maclaurin
       sum may use the tree algorithm */
    for (iter = 0; iter < 10 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, a;
        acb_ptr z1, z2;
        slong i, len, prec1, prec2;

        acb_init(s);
        acb_init(a);

        len = 100 + n_randint(state, 100);
        prec1 = 400 + n_randint(state, 400);
        prec2 = prec1 + 100;

        arb_set_d(acb_realref(s), (double) n_randint(state, 400) / 100 - 1.0);
        if (n_randint(state, 2))
            arb_set_d(acb_imagref(s), (double) n_randint(state, 3000) / 100);
        if (acb_is_one(s))
            acb_set_si(s, 2);

        if (n_randint(state, 2))
            acb_one(a);
        else
        {
            acb_set_ui(a, 1 + n_randint(state, 10));
            acb_div_ui(a, a, 1 + n_randint(state, 10), prec2);
        }

        z1 = _acb_vec_init(len);
        z2 = _acb_vec_init(len);

        flint_set_num_threads(1 + n_randint(state, 3));
        _acb_poly_zeta_cpx_series(z1, s, a, 0, len, prec1);
        _acb_poly_zeta_cpx_series(z2, s, a, 0, len, prec2);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(z1 + i, z2 + i) || !acb_is_finite(z1 + i))
            {
                flint_printf("FAIL: zeta overlap\n\n");
                flint_printf("iter = %wd, len = %wd, prec1 = %wd, i = %wd\n\n", iter, len, prec1, i);
                flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("z1 = "); acb_printd(z1 + i, 30); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2 + i, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(s);
        acb_clear(a);
        _acb_vec_clear(z1, len);
        _acb_vec_clear(z2, len);
    }

    /* polylogarithm series with small |z|, which may use the tree algorithm */
    for (iter = 0; iter < 10 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, z;
        acb_ptr w1, w2;
        slong i, len, prec1, prec2;

        acb_init(s);
        acb_init(z);

        len = 100 + n_randint(state, 100);
        prec1 = 400 + n_randint(state, 400);
        prec2 = prec1 + 100;

        arb_set_d(acb_realref(s), (double) n_randint(state, 400) / 100 - 1.0);
        if (n_randint(state, 2))
            arb_set_d(acb_imagref(s), (double) n_randint(state, 3000) / 100);

        acb_set_si(z, 1 + n_randint(state, 40));
        acb_div_ui(z, z, 100, prec2);
        if (n_randint(state, 2))
        {
            arb_set_si(acb_imagref(z), (slong) n_randint(state, 40) - 20);
            arb_div_ui(acb_imagref(z), acb_imagref(z), 100, prec2);
        }

        w1 = _acb_vec_init(len);
        w2 = _acb_vec_init(len);

        flint_set_num_threads(1 + n_randint(state, 3));
        _acb_poly_polylog_cpx(w1, s, z, len, prec1);
        _acb_poly_polylog_cpx(w2, s, z, len, prec2);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(w1 + i, w2 + i) || !acb_is_finite(w1 + i))
            {
                flint_printf("FAIL: polylog overlap\n\n");
                flint_printf("iter = %wd, len = %wd, prec1 = %wd, i = %wd\n\n", iter, len, prec1, i);
                flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("w1 = "); acb_printd(w1 + i, 30); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2 + i, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(s);
        acb_clear(z);
        _acb_vec_clear(w1, len);
        _acb_vec_clear(w2, len);
    }

    TEST_FUNCTION_END(state);
}
