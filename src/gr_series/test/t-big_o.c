/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "gr.h"
#include "gr_series.h"

/* Parse s, print it, parse the printout, and print again. The string form is
   canonical, so the two printouts must agree: this verifies that "... + O(v^n)"
   error terms survive a parse/print round-trip even when the elements are
   inexact (and hence not provably equal via gr_equal). */
static int
check_fixed_point(gr_ctx_t ctx, const char * s)
{
    int status;
    gr_ptr a, b;
    char * s2 = NULL;
    char * s3 = NULL;
    int ok;

    GR_TMP_INIT2(a, b, ctx);

    status = gr_set_str(a, s, ctx);
    if (status == GR_SUCCESS)
        status |= gr_get_str(&s2, a, ctx);
    if (status == GR_SUCCESS)
        status |= gr_set_str(b, s2, ctx);
    if (status == GR_SUCCESS)
        status |= gr_get_str(&s3, b, ctx);

    ok = (status == GR_SUCCESS) && (s2 != NULL) && (s3 != NULL) && !strcmp(s2, s3);

    flint_free(s2);
    flint_free(s3);
    GR_TMP_CLEAR2(a, b, ctx);

    return ok;
}

TEST_FUNCTION_START(gr_series_big_o, state)
{
    gr_ctx_t QQ, Rx, Rxy, Rxpoly, Rxsermod;
    slong i;

    const char * nested_ok[] = {
        "(3 + x + O(x^5)) + (7 + O(x^2))*y + (5 + O(x^3))*y^2 + O(y^3)",
        "(x + O(x^5)) + O(y^5)",
        "x + O(x^4) + O(y^3)",
        "O(y^3)",
        "O(x^2)",
        NULL,
    };

    gr_ctx_init_fmpq(QQ);

    gr_ctx_init_gr_series(Rx, QQ, 30);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rx, "x"));

    gr_ctx_init_gr_series(Rxy, Rx, 30);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rxy, "y"));

    gr_ctx_init_gr_poly(Rxpoly, Rx);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rxpoly, "y"));

    gr_series_mod_ctx_init(Rxsermod, Rx, 30);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rxsermod, "y"));

    /* All of these are representable in (QQ[[x]])[[y]]. */
    for (i = 0; nested_ok[i] != NULL; i++)
        if (!check_fixed_point(Rxy, nested_ok[i]))
            TEST_FUNCTION_FAIL("expected a stable round-trip in QQ[[x]][[y]]:\n%s\n", nested_ok[i]);

    /* Univariate: the error of x + O(x^5) is exactly 5. */
    {
        gr_series_t a;
        gr_series_init(a, Rx);
        if (gr_set_str(a, "x + O(x^5)", Rx) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("failed to parse x + O(x^5)\n");
        if (GR_SERIES_ERROR(a) != 5)
            TEST_FUNCTION_FAIL("expected error 5, got %wd\n", GR_SERIES_ERROR(a));
        gr_series_clear(a, Rx);
    }

    /* Nested: O(x^k) lands in the y^0 coefficient only; O(y^m) sets the outer
       error. For (x + O(x^5)) + O(y^5): outer error 5, y^0 coefficient (= x)
       carries inner error 5, and there is no y^1 term. */
    {
        gr_series_t a;
        const gr_series_struct * c0;

        gr_series_init(a, Rxy);
        if (gr_set_str(a, "(x + O(x^5)) + O(y^5)", Rxy) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("failed to parse (x + O(x^5)) + O(y^5)\n");

        if (GR_SERIES_ERROR(a) != 5)
            TEST_FUNCTION_FAIL("expected outer error 5, got %wd\n", GR_SERIES_ERROR(a));

        if (GR_SERIES_POLY(a)->length != 1)
            TEST_FUNCTION_FAIL("expected a single y-coefficient, got length %wd\n",
                               GR_SERIES_POLY(a)->length);

        c0 = (const gr_series_struct *) GR_SERIES_POLY(a)->coeffs;
        if (GR_SERIES_ERROR(c0) != 5)
            TEST_FUNCTION_FAIL("expected inner error 5 on the y^0 coefficient, got %wd\n",
                               GR_SERIES_ERROR(c0));

        gr_series_clear(a, Rxy);
    }

    /* In the polynomial ring (QQ[[x]])[y], O on a constant (inner) variable is
       fine, but O on the polynomial generator y has no representation. */
    if (!check_fixed_point(Rxpoly,
            "(3 + x + O(x^5)) + (7 + O(x^2))*y + (5 + O(x^3))*y^2"))
        TEST_FUNCTION_FAIL("expected a stable round-trip in QQ[[x]][y]\n");

    {
        gr_ptr t;
        GR_TMP_INIT(t, Rxpoly);
        if (gr_set_str(t, "1 + O(y^3)", Rxpoly) == GR_SUCCESS)
            TEST_FUNCTION_FAIL("O(y^n) must be rejected in the exact polynomial ring QQ[[x]][y]\n");
        GR_TMP_CLEAR(t, Rxpoly);
    }

    {
        gr_ptr t;
        GR_TMP_INIT(t, Rxsermod);
        if (gr_set_str(t, "1 + O(y^3)", Rxsermod) == GR_SUCCESS)
            TEST_FUNCTION_FAIL("O(y^n) must be rejected in the series mod ring QQ[[x]][y] / y^m\n");
        GR_TMP_CLEAR(t, Rxsermod);
    }

    {
        gr_series_t a;
        gr_series_init(a, Rx);
        if (gr_set_str(a, "O(2^5)", Rx) == GR_SUCCESS)
            TEST_FUNCTION_FAIL("O(2^5) should be unrepresentable over QQ\n");
        gr_series_clear(a, Rx);
    }


    gr_ctx_clear(Rxsermod);
    gr_ctx_clear(Rxpoly);
    gr_ctx_clear(Rxy);
    gr_ctx_clear(Rx);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
