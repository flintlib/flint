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
#include "decimal.h"
#include "gr_poly.h"
#include "gr_mat.h"

/* parse s in ctx and check that the result prints as expected */
static void
_check_parse(const char * s, const char * expected, gr_ctx_t ctx)
{
    decball_t x;
    char * t;
    int status;

    decball_init(x, ctx);
    status = gr_set_str(x, s, ctx);
    t = (status == GR_SUCCESS) ? decball_get_str(x, ctx) : NULL;

    if (status != GR_SUCCESS || strcmp(t, expected) != 0)
    {
        flint_printf("FAIL: parse\n");
        gr_ctx_println(ctx);
        flint_printf("s = %s\n", s);
        flint_printf("expected = %s\n", expected);
        flint_printf("got = %s (status %d)\n", t ? t : "(null)", status);
        flint_abort();
    }

    flint_free(t);
    decball_clear(x, ctx);
}

TEST_FUNCTION_START(decball_str, state)
{
    slong iter;

    /* explicit cases */
    {
        gr_ctx_t ctx;

        gr_ctx_init_decball(ctx, 20, 0);
        decimal_ctx_set_rad_prec(ctx, 3);

        _check_parse("1.5", "1.5", ctx);
        _check_parse("[1.5]", "1.5", ctx);
        _check_parse("(1.5)", "1.5", ctx);
        _check_parse("[1.5 +/- 0.01]", "[1.5 +/- 0.01]", ctx);
        _check_parse("1.5 +/- 0.01", "[1.5 +/- 0.01]", ctx);
        _check_parse("[ 1.5 +/- 0.01 ]", "[1.5 +/- 0.01]", ctx);
        _check_parse("[1.5+/-0.01]", "[1.5 +/- 0.01]", ctx);
        _check_parse("+/- 0.5", "[0 +/- 0.5]", ctx);
        _check_parse("[+/- 0.5]", "[0 +/- 0.5]", ctx);
        _check_parse("[-1.5 +/- 0.01]", "[-1.5 +/- 0.01]", ctx);
        _check_parse("[1.5 +/- -0.01]", "[1.5 +/- 0.01]", ctx);
        _check_parse("[1.25e-3 +/- 3e-9]", "[0.00125 +/- 3e-9]", ctx);
        _check_parse("[1.25e30 +/- 3e15]", "[1.25e30 +/- 3e15]", ctx);
        _check_parse("[1.25e30 +/- 3e5]", "[1.25e30 +/- 300000]", ctx);
        _check_parse("[1.25e30 +/- 3e6]", "[1.25e30 +/- 3e6]", ctx);
        _check_parse("[1 +/- 0.001]", "[1 +/- 0.001]", ctx);
        _check_parse("[1 +/- 0.0001]", "[1 +/- 1e-4]", ctx);
        _check_parse("[1.25e30 +/- 3e25]", "[1.25e30 +/- 3e25]", ctx);
        _check_parse("[1 +/- 1e-10]*2", "[2 +/- 2e-10]", ctx);
        _check_parse("[1 +/- 1e-10] + [1 +/- 1e-10]", "[2 +/- 2e-10]", ctx);
        _check_parse("([1 +/- 1e-10])*([1 +/- 1e-10])", "[1 +/- 2.01e-10]", ctx);
        _check_parse("[1 +/- 0.001] +/- 0.001", "[1 +/- 0.002]", ctx);
        _check_parse("1/3", "[0.33333333333333333333 +/- 3.34e-21]", ctx);
        _check_parse("(1/3) +/- 1e-25", "[0.33333333333333333333 +/- 3.35e-21]", ctx);
        /* the radius gets rounded up to the radius precision */
        _check_parse("[1 +/- 0.123456]", "[1 +/- 0.124]", ctx);
        /* a radius given as a ball uses its upper bound */
        _check_parse("[1 +/- [0.1 +/- 0.01]]", "[1 +/- 0.11]", ctx);
        /* a midpoint with too many digits is rounded, enlarging the radius */
        _check_parse("1.2345678901234567890123456789", "[1.234567890123456789 +/- 1.24e-20]", ctx);

        {
            decball_t x;
            decball_init(x, ctx);
            if (gr_set_str(x, "[1 +/- x]", ctx) == GR_SUCCESS ||
                gr_set_str(x, "[1 +/-]", ctx) == GR_SUCCESS ||
                gr_set_str(x, "[1 +/- 1", ctx) == GR_SUCCESS)
            {
                flint_printf("FAIL: expected parse failure\n");
                flint_abort();
            }
            decball_clear(x, ctx);
        }

        gr_ctx_clear(ctx);
    }

    /* random roundtrips: x -> string -> y must satisfy x in y, and be
       equal when the printed digits fit the working precision and
       there are no exponent limits */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decball_t x, y;
        char * s;
        int status;
        slong prec;

        gr_ctx_init_decball_randtest(ctx, state, 60);
        prec = DECIMAL_CTX_PREC(ctx);

        /* a radius below 10^emin can only be parsed if underflow is
           allowed (it is then replaced by 10^emin) */
        if (DECIMAL_CTX_HAS_EXP_LIMITS(ctx))
            decimal_ctx_set_flags(ctx, DECIMAL_CTX_FLAGS(ctx) | DECIMAL_ALLOW_UNDERFLOW);
        decball_init(x, ctx);
        decball_init(y, ctx);

        GR_MUST_SUCCEED(decball_randtest(x, state, ctx));

        s = decball_get_str(x, ctx);
        status = gr_set_str(y, s, ctx);

        if (status != GR_SUCCESS || !_decball_contains(y, x, ctx))
        {
            flint_printf("FAIL: string roundtrip (containment)\n");
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_printf("s = %s\n", s);
            flint_printf("y = %{gr} (%d)\n", y, ctx, status);
            flint_abort();
        }

        if (decfloat_digits(&x->mid, ctx) <= prec && DECIMAL_CTX_RAD_PREC(ctx) <= prec &&
            DECIMAL_CTX(ctx)->emin == WORD_MIN && DECIMAL_CTX(ctx)->emax == WORD_MAX &&
            (decfloat_equal(&x->mid, &y->mid, ctx) != T_TRUE || _decmag_cmp(&x->rad, &y->rad, ctx) != 0))
        {
            flint_printf("FAIL: string roundtrip (equality)\n");
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_printf("s = %s\n", s);
            flint_printf("y = %{gr}\n", y, ctx);
            flint_abort();
        }

        flint_free(s);
        decball_clear(x, ctx);
        decball_clear(y, ctx);
        gr_ctx_clear(ctx);
    }

    /* polynomials and matrices over decballs */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, pctx, mctx;
        gr_poly_t f, g;
        gr_mat_t A, B;
        char * s;
        int status;
        slong i, j, n;

        gr_ctx_init_decball(ctx, 1 + n_randint(state, 40), 0);
        decimal_ctx_set_rad_prec(ctx, 1 + n_randint(state, 9));
        gr_ctx_init_gr_poly(pctx, ctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        status = gr_poly_randtest(f, state, 1 + n_randint(state, 6), ctx);
        status |= gr_get_str(&s, f, pctx);
        status |= gr_set_str(g, s, pctx);

        if (status != GR_SUCCESS || f->length != g->length)
        {
            flint_printf("FAIL: polynomial roundtrip\n");
            gr_ctx_println(ctx);
            flint_printf("f = %s\n", s);
            flint_printf("g = %{gr}  (%d)\n", g, pctx, status);
            flint_abort();
        }

        for (i = 0; i < f->length; i++)
        {
            if (!_decball_contains(gr_poly_coeff_ptr(g, i, ctx), gr_poly_coeff_ptr(f, i, ctx), ctx))
            {
                flint_printf("FAIL: polynomial roundtrip (containment)\n");
                gr_ctx_println(ctx);
                flint_printf("f = %s\n", s);
                flint_printf("g = %{gr}\n", g, pctx);
                flint_abort();
            }
        }

        flint_free(s);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_ctx_clear(pctx);

        n = n_randint(state, 4);
        gr_ctx_init_matrix_ring(mctx, ctx, n);
        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        status = gr_mat_randtest(A, state, ctx);
        status |= gr_get_str(&s, A, mctx);
        status |= gr_set_str(B, s, mctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: matrix roundtrip\n");
            gr_ctx_println(ctx);
            flint_printf("A = %s\n", s);
            flint_printf("B = %{gr}  (%d)\n", B, mctx, status);
            flint_abort();
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (!_decball_contains(gr_mat_entry_ptr(B, i, j, ctx), gr_mat_entry_ptr(A, i, j, ctx), ctx))
                {
                    flint_printf("FAIL: matrix roundtrip (containment)\n");
                    gr_ctx_println(ctx);
                    flint_printf("A = %s\n", s);
                    flint_printf("B = %{gr}\n", B, mctx);
                    flint_abort();
                }
            }
        }

        flint_free(s);
        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_ctx_clear(mctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
