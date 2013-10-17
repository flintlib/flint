/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("mulmod_preinv....");
    fflush(stdout);

    /* Aliasing res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, res, t, f, finv;
        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(finv, ctx);
        fq_poly_init(res, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        do {
            fq_poly_randtest(f, state, n_randint(state, 50), ctx);
        } while (fq_poly_is_zero(f, ctx));
        if (a->length >= f->length)
            fq_poly_rem (a, a, f, ctx);
        if (b->length >= f->length)
            fq_poly_rem (b, b, f, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_mulmod_preinv(res, a, b, f, finv, ctx);
        fq_poly_mulmod_preinv(a, a, b, f, finv, ctx);

        result = (fq_poly_equal(res, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(finv, ctx);
        fq_poly_clear(res, ctx);
        fq_poly_clear(t, ctx);
        fq_ctx_clear(ctx);
    }

    /* Aliasing res and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, res, t, f, finv;
        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(finv, ctx);
        fq_poly_init(res, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        do {
            fq_poly_randtest(f, state, n_randint(state, 50), ctx);
        } while (fq_poly_is_zero(f, ctx));
        if (a->length >= f->length)
            fq_poly_rem (a, a, f, ctx);
        if (b->length >= f->length)
            fq_poly_rem (b, b, f, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_mulmod_preinv(res, a, b, f, finv, ctx);
        fq_poly_mulmod_preinv(b, a, b, f, finv, ctx);

        result = (fq_poly_equal(res, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(finv, ctx);
        fq_poly_clear(res, ctx);
        fq_poly_clear(t, ctx);
        fq_ctx_clear(ctx);
    }

    /* Aliasing res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, res, t, f, finv;
        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(finv, ctx);
        fq_poly_init(res, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        do {
            fq_poly_randtest(f, state, n_randint(state, 50), ctx);
        } while (fq_poly_is_zero(f, ctx));
        if (a->length >= f->length)
            fq_poly_rem (a, a, f, ctx);
        if (b->length >= f->length)
            fq_poly_rem (b, b, f, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_mulmod_preinv(res, a, b, f, finv, ctx);
        fq_poly_mulmod_preinv(f, a, b, f, finv, ctx);

        result = (fq_poly_equal(res, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(finv, ctx);
        fq_poly_clear(res, ctx);
        fq_poly_clear(t, ctx);
        fq_ctx_clear(ctx);
    }

    /* Aliasing res and finv */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, res, t, f, finv;
        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(finv, ctx);
        fq_poly_init(res, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        do {
            fq_poly_randtest(f, state, n_randint(state, 50), ctx);
        } while (fq_poly_is_zero(f, ctx));
        if (a->length >= f->length)
            fq_poly_rem (a, a, f, ctx);
        if (b->length >= f->length)
            fq_poly_rem (b, b, f, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_mulmod_preinv(res, a, b, f, finv, ctx);
        fq_poly_mulmod_preinv(finv, a, b, f, finv, ctx);

        result = (fq_poly_equal(res, finv, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("finv:\n"); fq_poly_print(finv, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(finv, ctx);
        fq_poly_clear(res, ctx);
        fq_poly_clear(t, ctx);
        fq_ctx_clear(ctx);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_poly_t a, b, res1, res2, t, f, finv;

        fq_ctx_t ctx;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(a, ctx);
        fq_poly_init(b, ctx);
        fq_poly_init(f, ctx);
        fq_poly_init(finv, ctx);
        fq_poly_init(res1, ctx);
        fq_poly_init(res2, ctx);
        fq_poly_init(t, ctx);

        fq_poly_randtest(a, state, n_randint(state, 50), ctx);
        fq_poly_randtest(b, state, n_randint(state, 50), ctx);
        do {
            fq_poly_randtest(f, state, n_randint(state, 50), ctx);
        } while (fq_poly_is_zero(f, ctx));
        if (a->length >= f->length)
            fq_poly_rem (a, a, f, ctx);
        if (b->length >= f->length)
            fq_poly_rem (b, b, f, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        fq_poly_mulmod_preinv(res1, a, b, f, finv, ctx);

        fq_poly_mul(res2, a, b, ctx);
        fq_poly_divrem(t, res2, res2, f, ctx);

        result = (fq_poly_equal(res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fq_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b:\n"); fq_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fq_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fq_poly_print(res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n"); fq_poly_print(res2, ctx), flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(a, ctx);
        fq_poly_clear(b, ctx);
        fq_poly_clear(f, ctx);
        fq_poly_clear(finv, ctx);
        fq_poly_clear(res1, ctx);
        fq_poly_clear(res2, ctx);
        fq_poly_clear(t, ctx);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);

    flint_printf("PASS\n");
    return 0;
}
