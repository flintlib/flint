/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_beta_lower_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, a1, b, b1, t, u;
        acb_poly_t z, w, wa1, wb1, pt, pu;
        slong prec, len, lena1, lenb1;
        int regularized;

        acb_init(a);
        acb_init(a1);
        acb_init(b);
        acb_init(b1);
        acb_init(t);
        acb_init(u);

        acb_poly_init(z);
        acb_poly_init(w);
        acb_poly_init(wa1);
        acb_poly_init(wb1);
        acb_poly_init(pt);
        acb_poly_init(pu);

        regularized = n_randint(state, 2);

        prec = 2 + n_randint(state, 100);

        acb_randtest_param(a, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_randtest_param(b, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_poly_randtest(z, state, 10, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_poly_randtest(w, state, 10, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_poly_randtest(wa1, state, 10, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_poly_randtest(wb1, state, 10, 1 + n_randint(state, 100), 1 + n_randint(state, 10));

        len = n_randint(state, 10);
        lena1 = n_randint(state, 10);
        lenb1 = n_randint(state, 10);

        acb_add_ui(a1, a, 1, prec);
        acb_add_ui(b1, b, 1, prec);

        acb_hypgeom_beta_lower_series(w, a, b, z, regularized, len, prec);
        acb_hypgeom_beta_lower_series(wa1, a1, b, z, regularized, lena1, prec);
        acb_hypgeom_beta_lower_series(wb1, a, b1, z, regularized, lenb1, prec);

        if (regularized)
        {
            acb_poly_scalar_mul(pt, wa1, a, prec);
            acb_poly_scalar_mul(pu, wb1, b, prec);
            acb_poly_add(pt, pt, pu, prec);
            acb_add(u, a, b, prec);
            acb_poly_scalar_div(pt, pt, u, prec);
        }
        else
        {
            acb_poly_add(pt, wa1, wb1, prec);
        }

        acb_poly_set(pu, w);
        acb_poly_truncate(pu, FLINT_MIN(FLINT_MIN(len, lena1), lenb1));
        acb_poly_truncate(pt, FLINT_MIN(FLINT_MIN(len, lena1), lenb1));

        if (!acb_poly_overlaps(pu, pt))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("len = %wd, lena1 = %wd, lenb1 = %wd\n\n", len, lena1, lenb1);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_poly_printd(w, 30); flint_printf("\n\n");
            flint_printf("wa1 = "); acb_poly_printd(wa1, 30); flint_printf("\n\n");
            flint_printf("wb1 = "); acb_poly_printd(wb1, 30); flint_printf("\n\n");
            flint_printf("pt = "); acb_poly_printd(pt, 30); flint_printf("\n\n");
            flint_abort();
        }

        /* test I(a,b;z) = 1-I(b,a,1-z) */

        if (!regularized)
        {
            acb_add(t, a, b, prec);
            acb_gamma(t, t, prec);
            acb_poly_scalar_mul(w, w, t, prec);
            acb_rgamma(t, a, prec);
            acb_poly_scalar_mul(w, w, t, prec);
            acb_rgamma(t, b, prec);
            acb_poly_scalar_mul(w, w, t, prec);
        }

        acb_poly_add_si(pt, z, -1, prec);
        acb_poly_neg(pt, pt);
        acb_hypgeom_beta_lower_series(pt, b, a, pt, 1, lena1, prec);
        acb_poly_add_si(pt, pt, -1, prec);
        acb_poly_neg(pt, pt);

        acb_poly_set(pu, w);
        acb_poly_truncate(pu, FLINT_MIN(len, lena1));
        acb_poly_truncate(pt, FLINT_MIN(len, lena1));

        if (!acb_poly_overlaps(pu, pt))
        {
            flint_printf("FAIL: symmetry\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("len = %wd, lena1 = %wd\n\n", len, lena1);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_poly_printd(w, 30); flint_printf("\n\n");
            flint_printf("pt = "); acb_poly_printd(pt, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(a1);
        acb_clear(b);
        acb_clear(b1);
        acb_clear(t);
        acb_clear(u);

        acb_poly_clear(z);
        acb_poly_clear(w);
        acb_poly_clear(wa1);
        acb_poly_clear(wb1);
        acb_poly_clear(pt);
        acb_poly_clear(pu);
    }

    TEST_FUNCTION_END(state);
}
