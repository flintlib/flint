/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_interpolate_fast, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i, n, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t P;
        arb_poly_t R, S;
        fmpq_t t, u;
        arb_ptr xs, ys;

        fmpq_poly_init(P);
        arb_poly_init(R);
        arb_poly_init(S);
        fmpq_init(t);
        fmpq_init(u);

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 5);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        fmpq_poly_randtest(P, state, 1 + n_randint(state, 20), qbits1);
        n = P->length;

        xs = _arb_vec_init(n);
        ys = _arb_vec_init(n);

        arb_poly_set_fmpq_poly(R, P, rbits1);

        if (n > 0)
        {
            fmpq_randtest(t, state, qbits2);
            arb_set_fmpq(xs, t, rbits2);

            for (i = 1; i < n; i++)
            {
                fmpq_randtest_not_zero(u, state, qbits2);
                fmpq_abs(u, u);
                fmpq_add(t, t, u);
                arb_set_fmpq(xs + i, t, rbits2);
            }
        }

        for (i = 0; i < n; i++)
            arb_poly_evaluate(ys + i, R, xs + i, rbits2);

        arb_poly_interpolate_fast(S, xs, ys, n, rbits3);

        if (!arb_poly_contains_fmpq_poly(S, P))
        {
            flint_printf("FAIL:\n");
            flint_printf("P = "); fmpq_poly_print(P); flint_printf("\n\n");
            flint_printf("R = "); arb_poly_printd(R, 15); flint_printf("\n\n");
            flint_printf("S = "); arb_poly_printd(S, 15); flint_printf("\n\n");
            flint_abort();
        }

        fmpq_poly_clear(P);
        arb_poly_clear(R);
        arb_poly_clear(S);
        fmpq_clear(t);
        fmpq_clear(u);
        _arb_vec_clear(xs, n);
        _arb_vec_clear(ys, n);
    }

    /* Degenerate nodes (repeated exact nodes, or overlapping inexact
       nodes): the weights must agree with the direct computation using
       arb_inv, and the interpolating polynomial must not be finite. */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i, n, height, m, prec;
        arb_ptr xs, ys, w, w2, tmp;
        arb_ptr * tree;
        arb_poly_t S;
        int exact;

        n = 2 + n_randint(state, 8);
        prec = 2 + n_randint(state, 200);
        exact = n_randint(state, 2);

        xs = _arb_vec_init(n);
        ys = _arb_vec_init(n);
        w = _arb_vec_init(n);
        w2 = _arb_vec_init(n);
        tmp = _arb_vec_init(n + 1);
        arb_poly_init(S);

        for (i = 0; i < n; i++)
        {
            arb_set_si(xs + i, (slong) n_randint(state, 100) - 50);
            arb_set_si(ys + i, (slong) n_randint(state, 100) - 50);
        }

        /* force a repeated node */
        i = n_randint(state, n - 1);
        arb_set(xs + i + 1 + n_randint(state, n - 1 - i), xs + i);

        if (!exact)
            for (i = 0; i < n; i++)
                mag_set_ui_2exp_si(arb_radref(xs + i), 1, -4);

        tree = _arb_poly_tree_alloc(n);
        _arb_poly_tree_build(tree, xs, n, prec);

        _arb_poly_interpolation_weights(w, tree, n, prec);

        /* reference: the former implementation */
        height = FLINT_CLOG2(n);
        m = WORD(1) << (height - 1);
        _arb_poly_mul_monic(tmp, tree[height-1], m + 1,
                            tree[height-1] + (m + 1), (n - m + 1), prec);
        _arb_poly_derivative(tmp, tmp, n + 1, prec);
        _arb_poly_evaluate_vec_fast_precomp(w2, tmp, n, tree, n, prec);
        for (i = 0; i < n; i++)
            arb_inv(w2 + i, w2 + i, prec);

        if (!_arb_vec_equal(w, w2, n) || _arb_vec_is_finite(w, n))
        {
            flint_printf("FAIL (degenerate weights)\n");
            flint_printf("n = %wd, exact = %d\n", n, exact);
            flint_printf("xs = "); _arb_vec_printn(xs, n, 10, 0); flint_printf("\n");
            flint_printf("w = "); _arb_vec_printn(w, n, 10, 0); flint_printf("\n");
            flint_printf("w2 = "); _arb_vec_printn(w2, n, 10, 0); flint_printf("\n");
            flint_abort();
        }

        _arb_poly_tree_free(tree, n);

        arb_poly_interpolate_fast(S, xs, ys, n, prec);

        if (S->length == 0 || _arb_vec_is_finite(S->coeffs, S->length))
        {
            flint_printf("FAIL (degenerate interpolation)\n");
            flint_printf("n = %wd, exact = %d\n", n, exact);
            flint_printf("xs = "); _arb_vec_printn(xs, n, 10, 0); flint_printf("\n");
            flint_printf("S = "); arb_poly_printd(S, 10); flint_printf("\n");
            flint_abort();
        }

        _arb_vec_clear(xs, n);
        _arb_vec_clear(ys, n);
        _arb_vec_clear(w, n);
        _arb_vec_clear(w2, n);
        _arb_vec_clear(tmp, n + 1);
        arb_poly_clear(S);
    }

    TEST_FUNCTION_END(state);
}
