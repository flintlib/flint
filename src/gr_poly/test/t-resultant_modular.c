/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr_poly.h"

/* random bivariate polynomial of length at most leny in y whose coefficients
   have length at most lenx in x and entries of at most `bits` bits */
static int
_randtest_bivariate_zq(gr_poly_t f, flint_rand_t state, slong leny,
                            slong lenx, flint_bitcnt_t bits,
                            gr_ctx_t ctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    gr_poly_t c;
    slong i, k;

    gr_poly_init(c, cctx);
    status |= gr_poly_zero(f, ctx);

    for (i = 0; i < leny; i++)
    {
        status |= gr_poly_randtest(c, state, lenx, cctx);

        for (k = 0; k < c->length; k++)
        {
            fmpz_t v;
            fmpz_init(v);
            fmpz_randtest(v, state, bits);
            status |= gr_set_fmpz(GR_ENTRY(c->coeffs, k, cctx->sizeof_elem), v, cctx);
            fmpz_clear(v);
        }

        _gr_poly_normalise(c, cctx);
        status |= gr_poly_set_coeff_scalar(f, i, c, ctx);
    }

    gr_poly_clear(c, cctx);

    return status;
}

/* multiply every coefficient in y by a fixed polynomial in x, so that the
   input has a nontrivial content in x */
static int
_bivariate_mul_content_zq(gr_poly_t f, const gr_poly_t c,
                               gr_ctx_t ctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    gr_poly_t res, t;
    slong i;

    gr_poly_init(res, ctx);
    gr_poly_init(t, cctx);

    for (i = 0; i < f->length; i++)
    {
        status |= gr_poly_mul(t,
            (gr_poly_struct *) GR_ENTRY(f->coeffs, i, ctx->sizeof_elem), c, cctx);
        status |= gr_poly_set_coeff_scalar(res, i, t, ctx);
    }

    status |= gr_poly_set(f, res, ctx);

    gr_poly_clear(res, ctx);
    gr_poly_clear(t, cctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_resultant_modular, state)
{
    slong iter;

    /* Compare with the division-free Sylvester determinant over Z[x][y] and
       Q[x][y]. The inputs are also given common factors and contents, which
       exercise the content extraction and the early termination. */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, g;
        gr_ptr x, y;
        int rational = n_randint(state, 2);
        int mode = n_randint(state, 4);
        int status = GR_SUCCESS;
        int s1;

        /* the images at the several primes are split over the thread pool */
        flint_set_num_threads(1 + n_randint(state, 8));

        if (rational)
            gr_ctx_init_fmpq(cctx);
        else
            gr_ctx_init_fmpz(cctx);

        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _randtest_bivariate_zq(f, state, 1 + n_randint(state, 5),
                        1 + n_randint(state, 5), n_randint(state, 60), ctx, cctx);
        status |= _randtest_bivariate_zq(g, state, 1 + n_randint(state, 5),
                        1 + n_randint(state, 5), n_randint(state, 60), ctx, cctx);

        if (mode == 1)
        {
            /* a common factor in y makes the resultant zero */
            gr_poly_t h;
            gr_poly_init(h, ctx);
            status |= _randtest_bivariate_zq(h, state, 2 + n_randint(state, 2),
                            1 + n_randint(state, 3), n_randint(state, 20), ctx, cctx);
            status |= gr_poly_mul(f, f, h, ctx);
            status |= gr_poly_mul(g, g, h, ctx);
            gr_poly_clear(h, ctx);
        }
        else if (mode == 2)
        {
            /* a content in x */
            gr_poly_t c;
            gr_poly_init(c, cctx);
            status |= gr_poly_randtest(c, state, 1 + n_randint(state, 4), cctx);
            if (c->length != 0)
                status |= _bivariate_mul_content_zq(f, c, ctx, cctx);
            gr_poly_clear(c, cctx);
        }

        /* the two settings of proved must agree with each other and with
           the Sylvester determinant */
        s1 = gr_poly_resultant_modular(x, f, g, n_randint(state, 2), ctx);
        status |= gr_poly_resultant_sylvester(y, f, g, ctx);

        if (s1 == GR_SUCCESS && status == GR_SUCCESS && gr_equal(x, y, ctx) == T_FALSE)
        {
            flint_printf("FAIL (vs sylvester):\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, fh, g, h;
        gr_ptr x, y, z, yz;
        int status = GR_SUCCESS;

        flint_set_num_threads(1 + n_randint(state, 8));

        if (n_randint(state, 2))
            gr_ctx_init_fmpq(cctx);
        else
            gr_ctx_init_fmpz(cctx);

        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(fh, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(h, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        yz = gr_heap_init(ctx);

        status |= _randtest_bivariate_zq(f, state, 1 + n_randint(state, 4),
                        1 + n_randint(state, 4), n_randint(state, 40), ctx, cctx);
        status |= _randtest_bivariate_zq(g, state, 1 + n_randint(state, 4),
                        1 + n_randint(state, 4), n_randint(state, 40), ctx, cctx);
        status |= _randtest_bivariate_zq(h, state, 1 + n_randint(state, 4),
                        1 + n_randint(state, 4), n_randint(state, 40), ctx, cctx);

        status |= gr_poly_mul(fh, f, h, ctx);

        status |= gr_poly_resultant_modular(x, fh, g, n_randint(state, 2), ctx);
        status |= gr_poly_resultant_modular(y, f, g, n_randint(state, 2), ctx);
        status |= gr_poly_resultant_modular(z, h, g, n_randint(state, 2), ctx);
        status |= gr_mul(yz, y, z, ctx);

        if (status == GR_SUCCESS && gr_equal(x, yz, ctx) == T_FALSE)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("h = "); gr_poly_print(h, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("yz = "); gr_println(yz, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(fh, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(h, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_heap_clear(yz, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* Check that the generic resultant dispatches here and agrees with the
       subresultant algorithm at degrees above the cutoff */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, g;
        gr_ptr x, y;
        int status = GR_SUCCESS;

        flint_set_num_threads(1 + n_randint(state, 8));

        if (n_randint(state, 2))
            gr_ctx_init_fmpq(cctx);
        else
            gr_ctx_init_fmpz(cctx);

        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _randtest_bivariate_zq(f, state, 7 + n_randint(state, 4),
                        1 + n_randint(state, 6), 1 + n_randint(state, 80), ctx, cctx);
        status |= _randtest_bivariate_zq(g, state, 7 + n_randint(state, 4),
                        1 + n_randint(state, 6), 1 + n_randint(state, 80), ctx, cctx);

        status |= gr_poly_resultant(x, f, g, ctx);
        status |= gr_poly_resultant_subresultant(y, f, g, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: did not succeed over Z[x][y] or Q[x][y]\n\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (gr_equal(x, y, ctx) == T_FALSE)
        {
            flint_printf("FAIL (vs subresultant):\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
