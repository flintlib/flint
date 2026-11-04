/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"
#include "../impl.h"

/* checks an enclosure x of x_k and w of its weight: P_n changes sign
   across x, x overlaps the root at low precision (by Newton iteration
   with evaluations of P_n), w overlaps the weight evaluated at x, and
   both are accurate */
static void
check_root(const char * method, const arb_t x, const arb_t w, ulong n, ulong k, slong prec)
{
    arb_t a, b, pa, pb, x0, t, u;
    slong wp = 2 * prec + 2 * n + 64;

    arb_init(a);
    arb_init(b);
    arb_init(pa);
    arb_init(pb);
    arb_init(x0);
    arb_init(t);
    arb_init(u);

    arb_get_lbound_arf(arb_midref(a), x, wp);
    arb_get_ubound_arf(arb_midref(b), x, wp);
    arb_hypgeom_legendre_p_ui(pa, NULL, n, a, wp);
    arb_hypgeom_legendre_p_ui(pb, NULL, n, b, wp);
    arb_mul(t, pa, pb, wp);

    /* at low precision, by Newton iteration with evaluations of P_n */
    arb_hypgeom_legendre_p_ui_root(x0, NULL, n, k, 64);

    arb_hypgeom_legendre_p_ui(NULL, u, n, x, wp);
    arb_mul(u, u, u, wp);

    {
        arb_t one_m_x2;
        arb_init(one_m_x2);
        arb_mul(one_m_x2, x, x, wp);
        arb_sub_ui(one_m_x2, one_m_x2, 1, wp);
        arb_neg(one_m_x2, one_m_x2);
        arb_mul(u, u, one_m_x2, wp);
        arb_ui_div(u, 2, u, wp);
        arb_clear(one_m_x2);
    }

    if (!arb_is_negative(t) || !arb_overlaps(x, x0) || !arb_overlaps(w, u)
        || arb_rel_accuracy_bits(x) < prec - 4 || arb_rel_accuracy_bits(w) < prec - 4)
    {
        flint_printf("FAIL (%s)\n\n", method);
        flint_printf("n = %wu, k = %wu, prec = %wd\n\n", n, k, prec);
        flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
        flint_printf("w = "); arb_printd(w, 30); flint_printf("\n\n");
        flint_printf("P(a) P(b) = "); arb_printd(t, 10); flint_printf("\n\n");
        flint_printf("x0 = "); arb_printd(x0, 30); flint_printf("\n\n");
        flint_printf("w(x) = "); arb_printd(u, 30); flint_printf("\n\n");
        flint_abort();
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(pa);
    arb_clear(pb);
    arb_clear(x0);
    arb_clear(t);
    arb_clear(u);
}

static void
report_failure(const char * method, ulong n, ulong k, slong prec)
{
    flint_printf("FAIL (%s returned 0)\n\n", method);
    flint_printf("n = %wu, k = %wu, prec = %wd\n\n", n, k, prec);
    flint_abort();
}

TEST_FUNCTION_START(arb_hypgeom_legendre_p_ui_root_methods, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        ulong n, k;
        slong prec;
        arb_t x, w;

        arb_init(x);
        arb_init(w);

        if (n_randint(state, 10) == 0)
            n = 2 + n_randint(state, 3000);
        else
            n = 2 + n_randint(state, 300);

        k = n_randint(state, n / 2);

        if (n_randint(state, 4) == 0)
            prec = 64 + n_randint(state, 6000);
        else
            prec = 64 + n_randint(state, 1000);

        /* the methods may fail (return 0) in principle, but should not
           at the precisions where they are used */
        if (_arb_hypgeom_legendre_p_ui_root_bitburst(x, w, n, k, prec))
            check_root("bitburst", x, w, n, k, prec);
        else if (prec >= 128)
            report_failure("bitburst", n, k, prec);

        /* the block variants of stage 0 and the transfers, and the FFT
           products, regardless of size */
        if (_arb_hypgeom_legendre_p_ui_root_bitburst_param(x, w, n, k, prec, 0, 0.0, 0))
            check_root("bitburst, blocks", x, w, n, k, prec);
        else if (prec >= 128)
            report_failure("bitburst, blocks", n, k, prec);

        if (_arb_hypgeom_legendre_p_ui_root_taylor(x, w, n, k, prec))
            check_root("taylor", x, w, n, k, prec);
        else if (prec >= 128)
            report_failure("taylor", n, k, prec);

        /* the dispatch, including the negative roots */
        k = n_randint(state, n);
        arb_hypgeom_legendre_p_ui_root(x, w, n, k, prec);
        if (n % 2 == 1 && k == n / 2)
        {
            if (!arb_is_zero(x))
            {
                flint_printf("FAIL (middle root)\n\n");
                flint_abort();
            }
        }
        else
        {
            if (k >= n / 2)
            {
                arb_neg(x, x);
                k = n - 1 - k;
            }
            check_root("dispatch", x, w, n, k, prec);
        }

        arb_clear(x);
        arb_clear(w);
    }

    TEST_FUNCTION_END(state);
}
