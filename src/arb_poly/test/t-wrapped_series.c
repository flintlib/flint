/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "arb_poly.h"

/* Test that the arb_poly_XXX_series implementations agree with gr_poly_XXX_series.
 * which: 0=acosh, 1=asinh, 2=atanh, 3=cot, 4=coth, 5=tan_pi, 6=tanh */
void
test_series(flint_rand_t state, int which)
{
    slong prec, n, res_len, res_gr_len;
    arb_poly_t res, x;
    gr_poly_t res_gr, x_gr;
    gr_ctx_t ctx;
    int status = GR_SUCCESS;

    prec = 2 + n_randint(state, 200);

    gr_ctx_init_real_arb(ctx, prec);

    arb_poly_init(res);
    arb_poly_init(x);
    gr_poly_init(res_gr, ctx);
    gr_poly_init(x_gr, ctx);

    n = n_randint(state, 4);

    arb_poly_randtest(x, state, n, prec, 10);
    gr_poly_fit_length(x_gr, arb_poly_length(x), ctx);
    status |= _gr_vec_set(x_gr->coeffs, x->coeffs, arb_poly_length(x), ctx);
    _gr_poly_set_length(x_gr, arb_poly_length(x), ctx);
    _gr_poly_normalise(x_gr, ctx);

    switch (which)
    {
        case 0:
            arb_poly_acosh_series(res, x, n, prec);
            status |= gr_poly_acosh_series(res_gr, x_gr, n, ctx);
            break;
        case 1:
            arb_poly_asinh_series(res, x, n, prec);
            status |= gr_poly_asinh_series(res_gr, x_gr, n, ctx);
            break;
        case 2:
            arb_poly_atanh_series(res, x, n, prec);
            status |= gr_poly_atanh_series(res_gr, x_gr, n, ctx);
            break;
        case 3:
            arb_poly_cot_series(res, x, n, prec);
            status |= gr_poly_cot_series(res_gr, x_gr, n, ctx);
            break;
        case 4:
            arb_poly_coth_series(res, x, n, prec);
            status |= gr_poly_coth_series(res_gr, x_gr, n, ctx);
            break;
        case 5:
            arb_poly_tan_pi_series(res, x, n, prec);
            status |= gr_poly_tan_pi_series(res_gr, x_gr, n, ctx);
            break;
        case 6:
            arb_poly_tanh_series(res, x, n, prec);
            status |= gr_poly_tanh_series(res_gr, x_gr, n, ctx);
            break;
    }

    res_len = arb_poly_length(res);
    res_gr_len = gr_poly_length(res_gr, ctx);

    if (status == GR_SUCCESS)
    {
        if (res_len != res_gr_len)
        {
            flint_printf("FAIL\n");
            flint_printf("Different lengths\n");
            flint_printf("status = %d\n", status);
            flint_printf("which = %d\n", which);
            flint_printf("n = %wd\n", n);
            printf("x = "); arb_poly_printd(x, 5); printf("\n");
            printf("res = "); arb_poly_printd(res, 5); printf("\n");
            printf("res_gr = "); gr_poly_print(res_gr, ctx); printf("\n");
            flint_abort();
        }

        if (!_arb_poly_overlaps(res->coeffs, res_len, res_gr->coeffs, res_gr_len))
        {
            flint_printf("FAIL\n");
            flint_printf("which = %d\n", which);
            flint_printf("n = %wd\n", n);
            printf("x = "); arb_poly_printd(x, 5); printf("\n");
            printf("res = "); arb_poly_printd(res, 5); printf("\n");
            printf("res_gr = "); gr_poly_print(res_gr, ctx); printf("\n");
            flint_abort();
        }
    }
    else
    {
        slong k;

        for (k = 0; k < res_len; k++)
        {
            if (!arf_is_nan(arb_midref(res->coeffs + k)))
            {
                flint_printf("FAIL\n");
                flint_printf("which = %d\n", which);
                flint_printf("n = %wd\n", n);
                flint_printf("k = %wd\n", k);
                printf("x = "); arb_poly_printd(x, 5); printf("\n");
                printf("res = "); arb_poly_printd(res, 5); printf("\n");
                printf("res_gr = "); gr_poly_print(res_gr, ctx); printf("\n");
                flint_abort();
            }
        }
    }

    arb_poly_clear(res);
    arb_poly_clear(x);
    gr_poly_clear(res_gr, ctx);
    gr_poly_clear(x_gr, ctx);
    gr_ctx_clear(ctx);
}

TEST_FUNCTION_START(arb_poly_wrapped_series, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        test_series(state, 0);
        test_series(state, 1);
        test_series(state, 2);
        test_series(state, 3);
        test_series(state, 4);
        test_series(state, 5);
        test_series(state, 6);
    }

    TEST_FUNCTION_END(state);
}
