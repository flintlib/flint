/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"
#include "arith.h"

TEST_FUNCTION_START(arb_hypgeom_legendre_p_ui_root, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        ulong n, k;
        slong prec;
        arb_ptr roots, weights;
        arb_poly_t pol;
        arb_t s;
        fmpq_poly_t pol2;

        n = 1 + n_randint(state, 100);
        prec = 20 + n_randint(state, 500);

        roots = _arb_vec_init(n);
        weights = _arb_vec_init(n);
        arb_poly_init(pol);
        fmpq_poly_init(pol2);
        arb_init(s);

        for (k = 0; k < n; k++)
        {
            if (k > n / 2 && n_randint(state, 2))
            {
                arb_neg(roots + k, roots + n - k - 1);
                arb_set(weights + k, weights + n - k - 1);
            }
            else
            {
                arb_hypgeom_legendre_p_ui_root(roots + k, weights + k, n, k, prec);
            }
        }

        arb_poly_product_roots(pol, roots, n, prec);
        /* fmpq_poly_legendre_p(pol2, n); */
        arith_legendre_polynomial(pol2, n);
        arb_set_fmpz(s, pol2->coeffs + n);
        arb_div_fmpz(s, s, pol2->den, prec);
        arb_poly_scalar_mul(pol, pol, s, prec);

        if (!arb_poly_contains_fmpq_poly(pol, pol2))
        {
            flint_printf("FAIL: polynomial containment\n\n");
            flint_printf("n = %wu, prec = %wd\n\n", n, prec);
            flint_printf("pol = "); arb_poly_printd(pol, 30); flint_printf("\n\n");
            flint_printf("pol2 = "); fmpq_poly_print(pol2); flint_printf("\n\n");
            flint_abort();
        }

        arb_zero(s);
        for (k = 0; k < n; k++)
        {
            arb_add(s, s, weights + k, prec);
        }

        if (!arb_contains_si(s, 2))
        {
            flint_printf("FAIL: sum of weights\n\n");
            flint_printf("n = %wu, prec = %wd\n\n", n, prec);
            flint_printf("s = "); arb_printn(s, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        _arb_vec_clear(roots, n);
        _arb_vec_clear(weights, n);
        arb_poly_clear(pol);
        fmpq_poly_clear(pol2);
        arb_clear(s);
    }

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x1, x2, w1, w2;
        ulong n, k;
        slong prec1, prec2;

        arb_init(x1);
        arb_init(x2);
        arb_init(w1);
        arb_init(w2);

        n = 1 + n_randtest(state) % 100000;
        if (n_randint(state, 2) || n == 1)
            k = n_randtest(state) % n;
        else
            k = n / 2 - (n_randtest(state) % (n / 2));

        prec1 = 2 + n_randtest(state) % 2000;
        prec2 = 2 + n_randtest(state) % 2000;

        arb_hypgeom_legendre_p_ui_root(x1, w1, n, k, prec1);
        if (n_randint(state, 10) == 0)
            arb_hypgeom_legendre_p_ui_root(x1, NULL, n, k, prec1);

        arb_hypgeom_legendre_p_ui_root(x2, w2, n, k, prec2);

        if (!arb_overlaps(x1, x2) || !arb_overlaps(w1, w2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu, k = %wu, prec1 = %wd, prec2 = %wd\n\n", n, k, prec1, prec2);
            flint_printf("x1 = "); arb_printn(x1, 100, 0); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printn(x2, 100, 0); flint_printf("\n\n");
            flint_printf("w1 = "); arb_printn(w1, 100, 0); flint_printf("\n\n");
            flint_printf("w2 = "); arb_printn(w2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(x1) < prec1 - 3 || arb_rel_accuracy_bits(w1) < prec1 - 3)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("n = %wu, k = %wu, prec1 = %wd\n\n", n, k, prec1);
            flint_printf("acc(x1) = %wd, acc(w1) = %wd\n\n", arb_rel_accuracy_bits(x1), arb_rel_accuracy_bits(w1));
            flint_printf("x1 = "); arb_printn(x1, prec1, ARB_STR_CONDENSE * 30); flint_printf("\n\n");
            flint_printf("w1 = "); arb_printn(w1, prec1, ARB_STR_CONDENSE * 30); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(x2) < prec2 - 3 || arb_rel_accuracy_bits(w2) < prec2 - 3)
        {
            flint_printf("FAIL: accuracy 2\n\n");
            flint_printf("n = %wu, k = %wu, prec2 = %wd\n\n", n, k, prec2);
            flint_printf("acc(x2) = %wd, acc(w2) = %wd\n\n", arb_rel_accuracy_bits(x2), arb_rel_accuracy_bits(w2));
            flint_printf("x2 = "); arb_printn(x2, prec2, ARB_STR_CONDENSE * 30); flint_printf("\n\n");
            flint_printf("w2 = "); arb_printn(w2, prec2, ARB_STR_CONDENSE * 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x1);
        arb_clear(x2);
        arb_clear(w1);
        arb_clear(w2);
    }

    TEST_FUNCTION_END(state);
}
